# S3 methods for "density_lp_ppp" -----------------------------------------------

#' Print a density_lp_ppp object
#' @param x A `"density_lp_ppp"` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.density_lp_ppp <- function(x, ...) {
  cat("Local polynomial density estimate (spatstat)\n")
  cat(sprintf("  h      : %g\n", x$h))
  cat(sprintf("  m      : %d\n", x$m))
  cat(sprintf("  N_quad : %d\n", x$N_quad))
  cat(sprintf("  Grid   : %d x %d pixels\n", length(x$xcol), length(x$yrow)))
  rng <- range(x$v, na.rm = TRUE)
  cat(sprintf("  Range  : [%.5g, %.5g]\n", rng[1L], rng[2L]))
  invisible(x)
}

#' Plot a density_lp_ppp object
#' @param x A `"density_lp_ppp"` object.
#' @param ... Passed to [spatstat.geom::plot.im()].
#' @return `x` invisibly.
#' @export
plot.density_lp_ppp <- function(x, ...) {
  NextMethod()
  invisible(x)
}

# S3 methods for "density_lp" ---------------------------------------------------

#' Print a density_lp object
#' @param x A `"density_lp"` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.density_lp <- function(x, ...) {
  cat("Local polynomial density estimate\n")
  cat(sprintf("  Domain : %s\n", x$domain$label))
  cat(sprintf("  d      : %d\n", x$domain$d))
  cat(sprintf("  h      : %g\n", x$h))
  cat(sprintf("  m      : %d\n", x$m))
  cat(sprintf("  N_quad : %d\n", x$N_quad))
  cat(sprintf("  Points : %d\n", nrow(x$t_grid)))
  cat(sprintf(
    "  Range  : [%.5g, %.5g]\n",
    min(x$estimate, na.rm = TRUE),
    max(x$estimate, na.rm = TRUE)
  ))
  invisible(x)
}

#' Plot a density_lp object
#' @param x A `"density_lp"` object.
#' @param main Plot title (default: bandwidth and degree).
#' @param ... Passed to the underlying plot function.
#' @return `x` invisibly.
#' @importFrom grDevices grey
#' @importFrom graphics image
#' @export
plot.density_lp <- function(
  x,
  main = paste0("densityLP  (h = ", x$h, ", m = ", x$m, ")"),
  ...
) {
  d <- x$domain$d
  if (d == 1L) {
    plot(
      x$t_grid[, 1L],
      x$estimate,
      type = "l",
      xlab = "t",
      ylab = "Estimated density",
      main = main,
      ...
    )
  } else if (d == 2L) {
    xs <- sort(unique(x$t_grid[, 1L]))
    ys <- sort(unique(x$t_grid[, 2L]))
    if (length(xs) * length(ys) == nrow(x$t_grid)) {
      Z <- matrix(x$estimate, nrow = length(xs), ncol = length(ys))
      image(xs, ys, Z, xlab = "x", ylab = "y", main = main, ...)
    } else {
      warning("plot.density_lp: irregular grid, falling back to scatter plot.")
      plot(
        x$t_grid[, 1L],
        x$t_grid[, 2L],
        col = grey(1 - x$estimate / max(x$estimate, na.rm = TRUE)),
        pch = 15,
        cex = 0.5,
        xlab = "x",
        ylab = "y",
        main = main,
        ...
      )
    }
  } else {
    warning("plot.density_lp: no graphical method for d = ", d, ".")
  }
  invisible(x)
}
