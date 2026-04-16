# Domain objects "lp_domain" ----------------------------------------------------

new_lp_domain <- function(d, sampler_factory, label) {
  structure(
    list(d = d, sampler_factory = sampler_factory, label = label),
    class = "lp_domain"
  )
}

#' @export
print.lp_domain <- function(x, ...) {
  cat("lp_domain object:", x$label, "(d =", x$d, ")\n")
  invisible(x)
}

# ── Exported constructors -------------------------------------------------------

#' Domain R^d (no constraint: neighbourhood = full \eqn{[-h, h]^d})
#'
#' @param d Space dimension (integer >= 1).
#' @return An `"lp_domain"` object.
#' @export
domain_Rd <- function(d) {
  stopifnot(is.numeric(d), length(d) == 1L, d >= 1, d == floor(d))
  d <- as.integer(d)

  # is_in_domain always returns TRUE: the full box [-h, h]^d is valid
  is_in <- function(X_mat) rep(TRUE, nrow(X_mat))
  new_lp_domain(d, function() sampler_rejection(is_in), label = "R^d")
}

#' Domain defined by an indicator function
#'
#' @param is_in_domain Vectorised function `f(X_mat)` where `X_mat` is an
#'   `n x d` numeric matrix of points in R^d and the return value is a logical
#'   vector of length `n`.
#' @param d Space dimension (integer >= 1).
#' @return An `"lp_domain"` object.
#' @export
domain_func <- function(is_in_domain, d) {
  stopifnot(is.function(is_in_domain))
  stopifnot(is.numeric(d), length(d) == 1L, d >= 1, d == floor(d))
  d <- as.integer(d)
  new_lp_domain(
    d,
    function() sampler_rejection(is_in_domain),
    label = "analytic domain"
  )
}

#' Polynomial sector \eqn{D_k = \{(x, y) : x \in [0, 1],\, 0 \le y \le x^k\}}
#'
#' @param k Sector exponent (scalar > 0).
#' @return An `"lp_domain"` object (d = 2).
#' @export
domain_sector <- function(k) {
  stopifnot(is.numeric(k), length(k) == 1L, k > 0)
  new_lp_domain(
    2L,
    function() sampler_sector(k),
    label = paste0("polynomial sector D_", k)
  )
}
