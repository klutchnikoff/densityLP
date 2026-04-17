# Point estimator (internal) ----------------------------------------------------

#' Estimate the density at a single point t using local polynomials
#'
#' @param X Numeric matrix n x d of observations.
#' @param t Numeric vector of length d (target point).
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0).
#' @param domain An `"lp_domain"` object.
#' @param N_quad Number of Monte Carlo quadrature points.
#' @return Scalar estimate of f(t).
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

  # 2. Gram matrix and Cholesky decomposition
  B <- gram_matrix(s, h, alphas)

  L <- tryCatch(
    chol(B),
    error = function(e) {
      stop(
        "The Gram matrix is not positive definite. ",
        "Try increasing N_quad or check the domain. Error: ",
        e$message
      )
    }
  )
  # B = t(L) %*% L  (L upper triangular)
  # B^{-1} = L^{-1} t(L^{-1})
  # e_1^T B^{-1} phi = (t(L^{-1}) e_1)^T (t(L^{-1}) phi)
  # Solve t(L) x = v  <=>  x = t(L^{-1}) v  (t(L) is lower triangular)

  # 3. H_gamma(0) = t(L^{-1}) %*% e_1
  # Phi_gamma(0) = e_1 because only the constant monomial equals 1 at u = 0
  e1 <- numeric(nrow(alphas))
  e1[1L] <- 1.0
  H_0 <- forwardsolve(t(L), e1) # D_m vector

  # 4. Observations in V(h): U_i = X_i - t, keep |U_i|_inf <= h
  U <- sweep(X, 2L, t, "-") # n x d
  in_V <- apply(abs(U) <= h, 1L, all)
  U_V <- t(U[in_V, , drop = FALSE]) # d x N_V

  if (ncol(U_V) == 0L) {
    return(0)
  }

  # 5. Estimator: (h^{-d} / n) * sum_i H_0^T H(U_i)
  Phi_V <- build_Phi(U_V, h, alphas) # D_m x N_V
  H_V <- forwardsolve(t(L), Phi_V) # D_m x N_V

  (h^(-d) / n) * sum(H_0 * H_V)
}
