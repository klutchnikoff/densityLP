# Domain objects "domain_lp" ----------------------------------------------------

#' @keywords internal
new_domain_lp <- function(d, sampler_factory, label) {
  structure(
    list(d = d, sampler_factory = sampler_factory, label = label),
    class = "domain_lp"
  )
}

#' Print a domain_lp object
#' @param x A `"domain_lp"` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.domain_lp <- function(x, ...) {
  cat("domain_lp object:", x$label, "(d =", x$d, ")\n")
  invisible(x)
}

# ── Exported constructors -------------------------------------------------------

#' Domain R^d (no constraint: neighbourhood = full \eqn{[-h, h]^d})
#'
#' @param d Space dimension (integer >= 1).
#' @return An `"domain_lp"` object.
#' @export
domain_Rd <- function(d) {
  stopifnot(is.numeric(d), length(d) == 1L, d >= 1, d == floor(d))
  d <- as.integer(d)

  is_in <- function(...) rep(TRUE, length(..1))
  new_domain_lp(d, function() sampler_rejection(is_in), label = "R^d")
}

#' Domain defined by an indicator function
#'
#' @param is_in_domain Vectorised indicator function of `d` arguments (one per
#'   coordinate axis), returning a logical vector of length `n`.  The dimension
#'   `d` is inferred from `length(formals(is_in_domain))`.
#' @return An `"domain_lp"` object.
#' @export
domain_from_indicator <- function(is_in_domain) {
  stopifnot(is.function(is_in_domain))
  d <- as.integer(length(formals(is_in_domain)))
  stopifnot(d >= 1L)
  new_domain_lp(
    d,
    function() sampler_rejection(is_in_domain),
    label = "analytic domain"
  )
}

#' Domain defined by a spatstat window (d = 2)
#'
#' Wraps a [spatstat.geom::owin()] object as an `"domain_lp"`.  The quadrature
#' sampler uses [spatstat.random::runifpoint()] on the exact intersection
#' \eqn{V(h) = \mathcal{D} \cap [-h, h]^2}, giving exact area and unbiased
#' Monte Carlo weights.
#'
#' @param win An `owin` object (spatstat.geom).
#' @return An `"domain_lp"` object (d = 2).
#' @export
domain_from_owin <- function(win) {
  stopifnot(inherits(win, "owin"))
  new_domain_lp(
    2L,
    function() sampler_owin(win),
    label = paste0("spatstat owin (", win$type, ")")
  )
}
