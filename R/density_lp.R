# Grid estimator (exported) -----------------------------------------------------

#' Local polynomial density estimation on a grid of target points
#'
#' @param X Numeric matrix `n x d` of observations. Must be a matrix; no
#'   implicit conversion from data frames or tibbles is performed.
#' @param t_grid Numeric matrix `p x d` of target points.
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0).
#' @param domain An `"domain_lp"` object describing the support of the density
#'   (see [domain_Rd()], [domain_from_indicator()], [domain_sector()]).
#' @param N_quad Number of Monte Carlo quadrature points per target point
#'   (default: 500).
#'
#' @return An S3 object of class `"density_lp"` containing:
#'   - `$estimate`: numeric vector of length `p` (estimated values),
#'   - `$variance`: numeric vector of length `p` (variance estimates),
#'   - `$t_grid`:   the matrix of target points,
#'   - `$h`, `$m`, `$domain`, `$N_quad`: parameters used.
#'
#' @export
density_lp <- function(X, t_grid, h, m = 0L, domain, N_quad = 500L) {
  check_X(X)
  check_t_grid(t_grid, ncol(X))
  check_domain(domain, ncol(X))
  check_h(h)
  check_m(m)
  check_N_quad(N_quad)

  p <- nrow(t_grid)
  estimate <- numeric(p)
  variance <- numeric(p)

  for (i in seq_len(p)) {
    res <- density_lp_point(
      X = X,
      t = t_grid[i, ],
      h = h,
      m = m,
      domain = domain,
      N_quad = N_quad
    )
    estimate[i] <- res["estimate"]
    variance[i] <- res["variance"]
  }

  structure(
    list(
      estimate = estimate,
      variance = variance,
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
