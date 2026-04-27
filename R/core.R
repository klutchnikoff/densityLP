# Point estimator (internal) ----------------------------------------------------

#' Estimate the density at a single point t using local polynomials
#'
#' @param X Numeric matrix n x d of observations.
#' @param t Numeric vector of length d (target point).
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0).
#' @param domain An `"domain_lp"` object.
#' @param N_quad Number of Monte Carlo quadrature points.
#' @return Named numeric vector of length 2: `estimate` and `variance`.
#' @keywords internal
density_lp_point <- function(X, t, h, m, domain, N_quad) {
  d <- ncol(X)
  n <- nrow(X)
  alphas <- build_alphas(m, d)

  # 1. Draw quadrature points in V(h) via the domain sampler
  sampler <- domain$sampler_factory()
  s <- sampler(N_quad, t, h)

  if (ncol(s$points) < 2L) {
    stop(
      "V(h) is too small or empty for h = ",
      h,
      " at t = (",
      paste(round(t, 4), collapse = ", "),
      ")."
    )
  }

  # 2–7. Gram matrix, Cholesky, observation filtering, and final summation in C++
  # U_obs: d x n matrix of centred observations (X_i - t), columns = observations
  U_obs <- t(X) - t
  res <- lp_estimator_loo_cpp(s$points, s$n_total, U_obs, h, alphas)

  # res[1] is f_hat(t)
  # res[2] is norm_H0_sq = ||H_0||^2
  # res[3] is v_hat(t) = (1/n^2) * sum (Z_i^2)

  c(estimate = res[1L], norm_H0_sq = res[2L], variance = res[3L])
}
