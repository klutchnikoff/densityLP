# Gram matrix -------------------------------------------------------------------

#' Compute the Gram matrix B_gamma by Monte Carlo approximation
#'
#' \eqn{B_\gamma = \int_{V(h)} \Phi \Phi^\top w_h \, du}
#'        ~= (vol(V(h)) / (N * h^d)) * sum_k phi(u_k) phi(u_k)^T
#'
#' @param s Output of a sampler: list with `$points` (d x N matrix of points
#'   u_k centred at t) and `$vol` (scalar, volume of V(h)).
#' @param h Bandwidth (scalar > 0).
#' @param alphas Integer matrix D_m x d of multi-indices.
#' @return Symmetric matrix of size D_m x D_m (positive definite when V(h) is
#'   non-degenerate).
#' @keywords internal
gram_matrix <- function(s, h, alphas) {
  stopifnot(is.list(s), !is.null(s$points), !is.null(s$vol))
  stopifnot(is.matrix(s$points), ncol(s$points) >= 1L)
  stopifnot(s$vol > 0)

  d <- nrow(s$points)
  N <- ncol(s$points)

  Phi <- build_Phi(s$points, h, alphas) # D_m x N
  weight <- s$vol / (N * h^d) # scalar
  weight * tcrossprod(Phi) # D_m x D_m
}
