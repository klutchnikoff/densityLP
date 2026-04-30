#' Unified Point Pattern for Local Polynomial Estimation
#'
#' Creates a unified object containing observations, a domain, and a bounding box.
#' This object serves as the standard input for estimation functions in the
#' \code{densityLP} package.
#'
#' @param X A matrix of observations (n x d).
#' @param domain A \code{domain_lp} object.
#' @param bbox A bounding box (matrix 2 x d, rows = min/max). If \code{NULL},
#'   it is inferred from the domain or the data.
#' @return An object of class \code{"pp_lp"} with fields \code{$X}, \code{$domain},
#'   \code{$bbox}, \code{$d}, and \code{$n}.
#' @examples
#' X <- matrix(runif(40), 20, 2)
#' dom <- domain_Rd(2L)
#' pp <- pp_lp(X, dom)
#' print(pp)
#' @export
pp_lp <- function(X, domain, bbox = NULL) {
  check_X(X)
  stopifnot(inherits(domain, "domain_lp"))

  d <- ncol(X)
  if (domain$d != d) {
    stop(sprintf(
      "Dimension mismatch: X has %d columns but domain is for d = %d",
      d,
      domain$d
    ))
  }

  if (is.null(bbox)) {
    bbox <- .get_bbox(X, domain)
  }

  structure(
    list(
      X = X,
      domain = domain,
      bbox = bbox,
      d = d,
      n = nrow(X)
    ),
    class = "pp_lp"
  )
}

#' @rdname pp_lp
#' @param x A \code{"pp_lp"} object.
#' @param ... Ignored.
#' @export
print.pp_lp <- function(x, ...) {
  cat(sprintf("Unified Point Pattern (d = %d)\n", x$d))
  cat(sprintf("  Observations : n = %d\n", x$n))
  cat(sprintf("  Domain       : %s\n", x$domain$label))
  cat("  Bounding Box :\n")
  print(x$bbox)
  invisible(x)
}

#' Convert an object to a pp_lp object
#'
#' @param x An object to convert. Can be a \code{matrix}, a \code{pp_lp}, or a
#'   \code{spatstat.geom::ppp}.
#' @param ... Additional arguments (e.g. \code{domain} for the matrix method).
#' @return An object of class \code{"pp_lp"}.
#' @examples
#' X <- matrix(runif(40), 20, 2)
#' dom <- domain_Rd(2L)
#' pp <- as_pp_lp(X, domain = dom)
#' print(pp)
#' @export
as_pp_lp <- function(x, ...) {
  UseMethod("as_pp_lp")
}

#' @rdname as_pp_lp
#' @export
as_pp_lp.pp_lp <- function(x, ...) {
  x
}

#' @rdname as_pp_lp
#' @param domain A \code{domain_lp} object (required for the matrix method).
#' @param bbox A bounding box matrix (optional).
#' @export
as_pp_lp.matrix <- function(x, domain, bbox = NULL, ...) {
  pp_lp(x, domain, bbox)
}

#' @rdname as_pp_lp
#' @importFrom spatstat.geom Window
#' @export
as_pp_lp.ppp <- function(x, ...) {
  X <- matrix(c(x$x, x$y), ncol = 2L)
  win <- spatstat.geom::Window(x)
  domain <- domain_from_owin(win)
  bbox <- matrix(c(win$xrange, win$yrange), nrow = 2L, byrow = FALSE)
  pp_lp(X, domain, bbox)
}

.get_bbox <- function(X, domain) {
  if (!is.null(domain$bbox)) {
    return(domain$bbox)
  }

  ranges <- apply(X, 2, range)
  extents <- ranges[2, ] - ranges[1, ]
  margin <- 0.05
  bbox <- rbind(
    ranges[1, ] - margin * extents,
    ranges[2, ] + margin * extents
  )
  colnames(bbox) <- colnames(X)
  rownames(bbox) <- c("min", "max")
  bbox
}
