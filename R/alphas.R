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

  grid <- do.call(expand.grid, rep(list(seq.int(0L, m)), d))
  valid <- rowSums(grid) <= m
  alphas <- as.matrix(grid[valid, , drop = FALSE])
  mode(alphas) <- "integer"

  # Sort by increasing total degree, then lexicographically column by column
  deg <- rowSums(alphas)
  ord <- do.call(
    order,
    c(list(deg), lapply(seq_len(d), function(j) alphas[, j]))
  )
  alphas <- alphas[ord, , drop = FALSE]

  rownames(alphas) <- NULL
  colnames(alphas) <- paste0("x", seq_len(d))
  alphas
}

#' Evaluate rescaled monomials Phi_gamma(u) for a collection of points
#'
#' For each point u_k (column of U), computes the vector
#' \eqn{\varphi_\alpha(u_k / h) = \prod_j (u_{kj} / h)^{\alpha_j}} for all multi-indices
#' alpha in `alphas`.
#'
#' @param U Numeric matrix of size d x N (N points of dimension d).
#' @param h Bandwidth (scalar > 0).
#' @param alphas Integer matrix D_m x d of multi-indices (output of
#'   `build_alphas`).
#' @return Numeric matrix of size D_m x N.
#' @keywords internal
build_Phi <- function(U, h, alphas) {
  stopifnot(is.matrix(U), is.numeric(U))
  stopifnot(is.numeric(h), length(h) == 1L, h > 0)
  stopifnot(is.matrix(alphas), nrow(alphas) >= 1L)
  stopifnot(nrow(U) == ncol(alphas))

  d <- nrow(U)
  N <- ncol(U)
  Dm <- nrow(alphas)

  U_sc <- U / h # d x N: rescaled coordinates

  # Phi[j, k] = prod_{l=1}^d (U_sc[l, k])^{alphas[j, l]}
  Phi <- matrix(NA_real_, Dm, N)
  for (j in seq_len(Dm)) {
    powered <- U_sc^alphas[j, ] # d x N, element-wise
    Phi[j, ] <- col_prods(powered)
  }
  Phi
}

# Column-wise product of a matrix (equivalent to apply(M, 2, prod) but avoids
# overhead for small d).
col_prods <- function(M) {
  if (nrow(M) == 1L) {
    return(M[1L, ])
  }
  out <- M[1L, ]
  for (i in seq.int(2L, nrow(M))) {
    out <- out * M[i, ]
  }
  out
}
