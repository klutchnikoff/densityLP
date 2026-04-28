# Grid estimator (exported) -----------------------------------------------------

#' Internal constructor for density_lp objects
#' @keywords internal
new_density_lp <- function(
  estimate,
  variance,
  t_grid,
  X,
  h,
  m,
  domain,
  N_quad,
  n_fail,
  call
) {
  structure(
    list(
      estimate = estimate,
      variance = variance,
      t_grid = t_grid,
      X = X,
      params = list(
        h = h,
        m = m,
        domain = domain,
        N_quad = N_quad
      ),
      stats = list(
        n_fail = n_fail,
        n_obs = nrow(X),
        p = nrow(t_grid),
        d = ncol(X)
      ),
      call = call
    ),
    class = "density_lp"
  )
}

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
#' @return An S3 object of class `"density_lp"`.
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
  n_fail <- 0L

  for (i in seq_len(p)) {
    res <- tryCatch(
      density_lp_point(
        X = X,
        t = t_grid[i, ],
        h = h,
        m = m,
        domain = domain,
        N_quad = N_quad
      ),
      error = function(e) NULL
    )
    if (is.null(res)) {
      estimate[i] <- NA_real_
      variance[i] <- NA_real_
      n_fail <- n_fail + 1L
    } else {
      estimate[i] <- res["estimate"]
      variance[i] <- res["variance"]
    }
  }

  if (n_fail > 0L) {
    warning(
      sprintf(
        "%d/%d grid point(s) failed (V(h) empty or Gram matrix singular). ",
        n_fail,
        p
      ),
      "Corresponding estimates are NA. Consider increasing N_quad or reducing h.",
      call. = FALSE
    )
  }

  new_density_lp(
    estimate = estimate,
    variance = variance,
    t_grid = t_grid,
    X = X,
    h = h,
    m = m,
    domain = domain,
    N_quad = N_quad,
    n_fail = n_fail,
    call = match.call()
  )
}
