# Samplers (internal) -----------------------------------------------------------
# Common interface: sampler(N, t, h) -> list(points = d x N, vol = scalar)
# Points are u = x - t (centred at t).
# All samplers return exactly N points.

#' Rejection sampler for an arbitrary analytic domain
#'
#' @param is_in_domain Vectorised function: n x d matrix -> logical vector.
#' @param n_pre Number of points used for the preliminary acceptance rate
#'   estimate.
#' @return Closure `function(N, t, h)`.
#' @importFrom stats runif
#' @keywords internal
sampler_rejection <- function(is_in_domain, n_pre = 2000L) {
  function(N, t, h) {
    d <- length(t)

    # Preliminary estimate of the acceptance rate and volume
    U_pre <- matrix(runif(n_pre * d, -h, h), n_pre, d)
    in_pre <- is_in_domain(sweep(U_pre, 2L, t, "+"))
    p_acc <- mean(in_pre)
    vol <- (2 * h)^d * p_acc

    if (p_acc < 1e-6) {
      stop(
        "The neighbourhood V(h) appears empty for h = ",
        h,
        ". Try reducing h or check the domain definition."
      )
    }

    # Rejection sampling until N points are accepted
    accepted <- matrix(NA_real_, d, 0L)
    while (ncol(accepted) < N) {
      need <- ceiling(1.5 * (N - ncol(accepted)) / p_acc)
      U <- matrix(runif(need * d, -h, h), need, d)
      ok <- is_in_domain(sweep(U, 2L, t, "+"))
      accepted <- cbind(accepted, t(U[ok, , drop = FALSE]))
    }
    list(points = accepted[, seq_len(N), drop = FALSE], vol = vol)
  }
}

#' Quasi-Monte Carlo sampler using a scrambled Sobol sequence
#'
#' Provides faster convergence than rejection sampling for smooth integrands:
#' \eqn{O((\log N)^d / N)} instead of \eqn{O(N^{-1/2})}.  Falls back to
#' standard Monte Carlo if the acceptance rate is too low.
#'
#' @param is_in_domain Vectorised function: n x d matrix -> logical vector.
#' @param scrambling Integer passed to [randtoolbox::sobol()] (default: 1).
#' @return Closure `function(N, t, h)`.
#' @importFrom randtoolbox sobol
#' @importFrom stats runif
#' @keywords internal
sampler_qmc <- function(is_in_domain, scrambling = 0L) {
  function(N, t, h) {
    d <- length(t)

    # Sobol sequence mapped to [-h, h]^d
    N_sob <- ceiling(3L * N)
    sob <- randtoolbox::sobol(N_sob, dim = d, scrambling = scrambling)
    U <- sweep(sob * 2 * h, 2L, h, "-") # N_sob x d

    in_V <- is_in_domain(sweep(U, 2L, t, "+"))
    p_acc <- mean(in_V)
    vol <- (2 * h)^d * p_acc
    pts <- t(U[in_V, , drop = FALSE]) # d x N_acc

    # Fall back to standard MC if acceptance rate is too low
    if (ncol(pts) < N) {
      warning(
        "QMC: low acceptance rate (",
        round(p_acc * 100, 1),
        "%), completing with standard MC."
      )
      while (ncol(pts) < N) {
        U2 <- matrix(runif(N * d, -h, h), N, d)
        ok <- is_in_domain(sweep(U2, 2L, t, "+"))
        pts <- cbind(pts, t(U2[ok, , drop = FALSE]))
      }
    }
    list(points = pts[, seq_len(N), drop = FALSE], vol = vol)
  }
}

#' Exact sampler for a polygonal domain (d = 2) using spatstat.geom
#'
#' Draws N uniform points in \eqn{V(h) = \mathcal{D} \cap [-h, h]^2} using
#' [spatstat.random::runifpoint()].  The area is computed exactly from the
#' polygon intersection.
#'
#' @param win An `owin` object (spatstat.geom) representing the domain
#'   \eqn{\mathcal{D}}.
#' @return Closure `function(N, t, h)`.
#' @importFrom spatstat.geom owin intersect.owin area.owin
#' @importFrom spatstat.random runifpoint
#' @keywords internal
sampler_spatstat <- function(win) {
  function(N, t, h) {
    box <- spatstat.geom::owin(
      c(t[1L] - h, t[1L] + h),
      c(t[2L] - h, t[2L] + h)
    )
    V_h <- spatstat.geom::intersect.owin(win, box)
    vol <- spatstat.geom::area.owin(V_h)

    if (vol < .Machine$double.eps) {
      stop(
        "V(h) has zero area for h = ",
        h,
        " at t = (",
        t[1L],
        ", ",
        t[2L],
        ")."
      )
    }

    pts <- spatstat.random::runifpoint(N, win = V_h)
    list(
      points = rbind(pts$x - t[1L], pts$y - t[2L]),
      vol = vol
    )
  }
}

#' Sampler for the polynomial sector D_k via rejection on the bounding box
#'
#' \eqn{D_k = \{(x, y) : x \in [0, 1],\, 0 \le y \le x^k\}}
#'
#' The volume of V(h) is computed exactly by 1D numerical integration.
#'
#' @param k Sector exponent (scalar > 0).
#' @return Closure `function(N, t, h)`.
#' @importFrom stats integrate runif
#' @keywords internal
sampler_sector <- function(k) {
  stopifnot(is.numeric(k), length(k) == 1L, k > 0)
  function(N, t, h) {
    tx <- t[1L]
    ty <- t[2L]
    x_lo <- max(0, tx - h)
    x_hi <- min(1, tx + h)

    if (x_lo >= x_hi) {
      stop(
        "V(h) is empty in the x direction for t = (",
        tx,
        ", ",
        ty,
        ") and h = ",
        h
      )
    }

    # Exact volume by 1D numerical integration
    vol <- integrate(
      function(x) pmax(0, pmin(x^k, ty + h) - pmax(0, ty - h)),
      x_lo,
      x_hi
    )$value

    if (vol < .Machine$double.eps) {
      stop("V(h) has zero volume for t = (", tx, ", ", ty, ") and h = ", h)
    }

    # Rejection sampling (efficient: domain is monotone in x)
    y_lo <- max(0, ty - h)
    y_hi <- ty + h
    accepted <- matrix(NA_real_, 2L, 0L)
    while (ncol(accepted) < N) {
      need <- ceiling(2 * (N - ncol(accepted)))
      x <- runif(need, x_lo, x_hi)
      y <- runif(need, y_lo, y_hi)
      ok <- (y >= 0) & (y <= x^k)
      accepted <- cbind(accepted, rbind(x[ok] - tx, y[ok] - ty))
    }
    list(points = accepted[, seq_len(N), drop = FALSE], vol = vol)
  }
}
