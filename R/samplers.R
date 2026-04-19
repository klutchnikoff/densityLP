# Samplers (internal) -----------------------------------------------------------
# Common interface: sampler(N_quad, t, h) -> list(points = d x N_in, n_total)
# Points are u = x - t (centred at t).
#
# $n_total : total draws in [-h,h]^d (exact, not estimated).
# $points  : all N_in accepted points — N_in is VARIABLE for MC samplers
#            (rejection) and fixed = N_quad for exact samplers.
#
# gram_matrix() uses weight = (2h)^d / (n_total * h^d).
# The domain information is carried by N_in (proportional to vol(V(h))):
#   E[ (2h)^d/n_total * sum_{k in V(h)} ... ] = B_gamma  (no vol needed).

#' Rejection sampler for an arbitrary analytic domain
#'
#' Draws `N_quad` points uniformly in \eqn{[-h,h]^d} and returns **all**
#' accepted points (those in \eqn{V(h)}).  The count `N_in` is therefore
#' variable and carries the domain-size information needed by `gram_matrix`.
#'
#' @param is_in_domain Vectorised function of `d` arguments (one per coordinate
#'   axis), returning a logical vector.
#' @return Closure `function(N_quad, t, h)`.
#' @importFrom stats runif
#' @keywords internal
sampler_rejection <- function(is_in_domain) {
  function(N_quad, t, h) {
    d <- length(t)
    U <- matrix(runif(N_quad * d, -h, h), N_quad, d)
    X <- sweep(U, 2L, t, "+")
    ok <- do.call(is_in_domain, lapply(seq_len(d), function(j) X[, j]))
    pts <- t(U[ok, , drop = FALSE])
    if (ncol(pts) == 0L) {
      stop(
        "V(h) is empty for h = ",
        h,
        ". Try increasing N_quad or reducing h."
      )
    }
    list(points = pts, n_total = N_quad)
  }
}

#' Exact sampler for a polygonal domain (d = 2) using spatstat.geom
#'
#' Draws N_quad uniform points in \eqn{V(h) = \mathcal{D} \cap [-h, h]^2}
#' using [spatstat.random::runifpoint()].  The area is computed exactly from
#' the polygon intersection; `n_total = round(N_quad * (2h)^2 / vol)`.
#'
#' @param win An `owin` object (spatstat.geom) representing the domain
#'   \eqn{\mathcal{D}}.
#' @return Closure `function(N_quad, t, h)`.
#' @importFrom spatstat.geom owin intersect.owin area.owin
#' @importFrom spatstat.random runifpoint
#' @keywords internal
sampler_owin <- function(win) {
  function(N_quad, t, h) {
    box <- spatstat.geom::owin(
      c(t[1L] - h, t[1L] + h),
      c(t[2L] - h, t[2L] + h)
    )
    V_h <- spatstat.geom::intersect.owin(win, box)
    vol <- spatstat.geom::area.owin(V_h)

    if (vol / (2 * h)^2 < 1e-4) {
      stop(sprintf(
        "V(h) has negligible area at t = (%.4g, %.4g) for h = %.4g: acceptance rate below 0.01%%, Monte Carlo estimates would be unreliable.",
        t[1L],
        t[2L],
        h
      ))
    }

    pts <- spatstat.random::runifpoint(N_quad, win = V_h)
    list(
      points = rbind(pts$x - t[1L], pts$y - t[2L]),
      n_total = round(N_quad * (2 * h)^2 / vol)
    )
  }
}
