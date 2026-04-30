# Multi-indices and rescaled monomials ------------------------------------------

#' Generate all multi-indices alpha in N_0^d with |alpha| <= m
#'
#' Multi-indices are ordered by increasing total degree, then lexicographically.
#' The number of rows equals D_m = choose(m + d, d).
#'
#' @param m Maximum degree (integer >= 0).
#' @param d Space dimension (integer >= 1).
#' @return Integer matrix of size D_m x d.
#' @keywords internal
build_alphas <- function(m, d) {
  stopifnot(is.numeric(m), length(m) == 1L, m >= 0, m == floor(m))
  stopifnot(is.numeric(d), length(d) == 1L, d >= 1, d == floor(d))

  grid <- do.call(expand.grid, rep(list(0:m), d))
  valid <- rowSums(grid) <= m
  alphas <- matrix(as.integer(as.matrix(grid[valid, , drop = FALSE])), ncol = d)

  # Sort by increasing total degree, then lexicographically column by column
  deg <- rowSums(alphas)
  ord <- do.call(order, as.data.frame(cbind(deg, alphas)))
  alphas <- alphas[ord, , drop = FALSE]

  rownames(alphas) <- NULL
  colnames(alphas) <- paste0("x", seq_len(d))
  alphas
}
