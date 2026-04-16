# Samplers (internal) -----------------------------------------------------------
# Common interface: sampler(N, t, h) -> list(points = d x N, vol = scalar)
# Points are u = x - t (centred at t).

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
