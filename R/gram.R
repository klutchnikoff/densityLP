# Gram matrix -------------------------------------------------------------------

#' Compute the Gram matrix B_gamma
#'
#' \eqn{B_\gamma = \int_{[-h,h]^d} \Phi(u) \Phi^\top(u)\, w_h(u) \, du
#'              \approx \frac{(2h)^d}{N_{total} \cdot h^d}
#'                      \sum_{k \in V(h)} \Phi(u_k)\Phi^\top(u_k)}
#'
#' `s$n_total` is the total number of draws in \eqn{[-h,h]^d} that produced
#' the `N` returned points. The box volume \eqn{(2h)^d} is always exact.
#'
#' @param s Output of a sampler: list with `$points` (d x N matrix of points
#'   u_k centred at t) and `$n_total` (integer, total draws in the box).
#' @param h Bandwidth (scalar > 0).
#' @param alphas Integer matrix D_m x d of multi-indices.
#' @return Symmetric matrix of size D_m x D_m (positive definite when V(h) is
#'   non-degenerate).
#' @keywords internal
gram_matrix <- function(s, h, alphas) {
  stopifnot(is.list(s), !is.null(s$points), !is.null(s$n_total))
  stopifnot(is.matrix(s$points), ncol(s$points) >= 1L)
  stopifnot(is.numeric(s$n_total), s$n_total > 0)

  gram_matrix_cpp(s$points, h, alphas, s$n_total)
}
