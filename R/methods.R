# S3 methods for "density_lp" ---------------------------------------------------

#' Print a density_lp object
#' @param x A `"density_lp"` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @importFrom grDevices grey
#' @importFrom graphics image
#' @export
print.density_lp <- function(x, ...) {
  cat("Local polynomial density estimate\n")
  cat("  Domain :", x$domain$label, "\n")
  cat("  d      :", x$domain$d, "\n")
  cat("  h      :", x$h, "\n")
  cat("  m      :", x$m, "\n")
  cat("  N_quad :", x$N_quad, "\n")
  cat("  Points :", nrow(x$t_grid), "\n")
  cat(
    "  Range  : [",
    round(min(x$estimate), 5),
    ",",
    round(max(x$estimate), 5),
    "]\n"
  )
  invisible(x)
}

#' Plot a density_lp object
#' @param x A `"density_lp"` object.
#' @param main Plot title (default: bandwidth and degree).
#' @param ... Passed to the underlying plot function.
#' @return `x` invisibly.
#' @export
plot.density_lp <- function(x,
                            main = paste0("densityLP  (h = ", x$h,
                                          ", m = ", x$m, ")"),
                            ...) {
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
        pch = 15, cex = 0.5,
        xlab = "x", ylab = "y", main = main,
        ...
      )
    }
  } else {
    warning("plot.density_lp: no graphical method for d = ", d, ".")
  }
  invisible(x)
}
