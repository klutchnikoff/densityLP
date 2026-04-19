# Polynomial sector domain -------------------------------------------------------

#' Polynomial sector \eqn{D_k = \{(x, y) : x \in [0, 1],\, 0 \le y \le x^k\}}
#'
#' @param k Sector exponent (scalar > 0).
#' @return An `"domain_lp"` object (d = 2).
#' @export
domain_sector <- function(k) {
  stopifnot(is.numeric(k), length(k) == 1L, k > 0)
  new_domain_lp(
    2L,
    function() sampler_sector(k),
    label = paste0("polynomial sector D_", k)
  )
}

#' @importFrom stats integrate runif
#' @keywords internal
sampler_sector <- function(k) {
  stopifnot(is.numeric(k), length(k) == 1L, k > 0)
  function(N_quad, t, h) {
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
    while (ncol(accepted) < N_quad) {
      need <- ceiling(2 * (N_quad - ncol(accepted)))
      x <- runif(need, x_lo, x_hi)
      y <- runif(need, y_lo, y_hi)
      ok <- (y >= 0) & (y <= x^k)
      accepted <- cbind(accepted, rbind(x[ok] - tx, y[ok] - ty))
    }
    list(
      points = accepted[, seq_len(N_quad), drop = FALSE],
      n_total = round(N_quad * (2 * h)^2 / vol)
    )
  }
}
