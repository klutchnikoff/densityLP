# Grid estimator (exported) -----------------------------------------------------

#' Local polynomial density estimation on a grid of target points
#'
#' @param X Numeric matrix `n x d` of observations. Must be a matrix; no
#'   implicit conversion from data frames or tibbles is performed.
#' @param t_grid Numeric matrix `p x d` of target points.
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0).
#' @param domain An `"lp_domain"` object describing the support of the density
#'   (see [domain_Rd()], [domain_func()], [domain_sector()]).
#' @param N_quad Number of Monte Carlo quadrature points per target point
#'   (default: 500).
#'
#' @return An S3 object of class `"density_lp"` containing:
#'   - `$estimate`: numeric vector of length `p` (estimated values),
#'   - `$t_grid`:   the matrix of target points,
#'   - `$h`, `$m`, `$domain`, `$N_quad`: parameters used.
#'
#' @export
density_lp <- function(X, t_grid, h, m, domain, N_quad = 500L) {
  check_X(X)
  check_t_grid(t_grid, ncol(X))
  check_domain(domain, ncol(X))
  stopifnot(is.numeric(h), length(h) == 1L, h > 0)
  stopifnot(is.numeric(m), length(m) == 1L, m >= 0, m == floor(m))
  stopifnot(is.numeric(N_quad), length(N_quad) == 1L, N_quad >= 1L)

  p <- nrow(t_grid)
  estimate <- numeric(p)

  for (i in seq_len(p)) {
    estimate[i] <- density_lp_point(
      X = X,
      t = t_grid[i, ],
      h = h,
      m = m,
      domain = domain,
      N_quad = N_quad
    )
  }

  structure(
    list(
      estimate = estimate,
      t_grid = t_grid,
      h = h,
      m = m,
      domain = domain,
      N_quad = N_quad,
      call = match.call()
    ),
    class = "density_lp"
  )
}
