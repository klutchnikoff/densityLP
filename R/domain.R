# Domain objects "domain_lp" ----------------------------------------------------

#' @keywords internal
new_domain_lp <- function(
  d,
  sampler_factory,
  is_in_domain,
  label,
  call,
  bbox = NULL
) {
  structure(
    list(
      d = d,
      sampler_factory = sampler_factory,
      is_in_domain = is_in_domain,
      label = label,
      call = call,
      bbox = bbox
    ),
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
#' @examples
#' dom <- domain_Rd(2)
#' print(dom)
#' @export
domain_Rd <- function(d) {
  stopifnot(is.numeric(d), length(d) == 1L, is.finite(d), d >= 1, d == floor(d))
  d <- as.integer(d)

  is_in_domain <- function(...) rep(TRUE, length(..1))
  new_domain_lp(
    d,
    sampler_factory = function() sampler_rejection(is_in_domain),
    is_in_domain = is_in_domain,
    label = "R^d",
    call = match.call()
  )
}

#' Domain defined by an indicator function
#'
#' @param is_in_domain Vectorised indicator function of `d` arguments (one per
#'   coordinate axis), returning a logical vector of length `n`.  The dimension
#'   `d` is inferred from `length(formals(is_in_domain))`.
#' @return An `"domain_lp"` object.
#' @examples
#' is_in_domain <- function(x, y) x^2 + y^2 <= 1
#' dom <- domain_from_indicator(is_in_domain)
#' print(dom)
#' @export
domain_from_indicator <- function(is_in_domain) {
  stopifnot(is.function(is_in_domain))
  d <- as.integer(length(formals(is_in_domain)))
  stopifnot(d >= 1L)
  new_domain_lp(
    d,
    sampler_factory = function() sampler_rejection(is_in_domain),
    is_in_domain = is_in_domain,
    label = "analytic domain",
    call = match.call()
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
#' @examples
#' win <- spatstat.geom::owin(c(0, 1), c(0, 1))
#' dom <- domain_from_owin(win)
#' print(dom)
#' @export
domain_from_owin <- function(win) {
  stopifnot(inherits(win, "owin"))
  bbox <- matrix(c(win$xrange, win$yrange), nrow = 2L, byrow = FALSE)
  new_domain_lp(
    2L,
    sampler_factory = function() sampler_owin(win),
    is_in_domain = function(x, y) spatstat.geom::inside.owin(x, y, win),
    label = paste0("spatstat owin (", win$type, ")"),
    call = match.call(),
    bbox = bbox
  )
}
