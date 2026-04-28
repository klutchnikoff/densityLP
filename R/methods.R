# S3 methods for "density_lp_ppp" -----------------------------------------------

#' Print a density_lp_ppp object
#' @param x A `"density_lp_ppp"` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.density_lp_ppp <- function(x, ...) {
  cat("Local polynomial density estimate (spatstat)\n")
  cat(sprintf("  Points : %d\n", x$stats$n_obs))
  cat(sprintf("  h      : %g\n", x$params$h))
  cat(sprintf("  m      : %d\n", x$params$m))
  cat(sprintf("  N_quad : %d\n", x$params$N_quad))
  cat(sprintf("  Grid   : %d x %d pixels\n", length(x$xcol), length(x$yrow)))
  rng <- range(x$v, na.rm = TRUE)
  cat(sprintf("  Range  : [%.5g, %.5g]\n", rng[1L], rng[2L]))
  if (!is.null(x$stats$n_fail) && x$stats$n_fail > 0L) {
    cat(sprintf(
      "  Note   : %d failure(s) during calculation\n",
      x$stats$n_fail
    ))
  }
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
  cat(sprintf("  Points : %d\n", x$stats$p))
  cat(sprintf("  h      : %g\n", x$params$h))
  cat(sprintf("  m      : %d\n", x$params$m))
  rng <- range(x$estimate, na.rm = TRUE)
  cat(sprintf("  Range  : [%.5g, %.5g]\n", rng[1L], rng[2L]))
  if (x$stats$n_fail > 0L) {
    cat(sprintf(
      "  Note   : %d failure(s) during calculation\n",
      x$stats$n_fail
    ))
  }
  invisible(x)
}

#' Summary of a density_lp object
#' @param object A `"density_lp"` object.
#' @param ... Ignored.
#' @return `object` invisibly.
#' @export
summary.density_lp <- function(object, ...) {
  cat("Local Polynomial Density Estimation Summary\n")
  cat("-------------------------------------------\n")
  cat(sprintf(
    "Call         : %s\n",
    paste(deparse(object$call), collapse = "\n")
  ))
  cat(sprintf(
    "Domain       : %s (d = %d)\n",
    object$params$domain$label,
    object$stats$d
  ))
  cat(sprintf("Observations : n = %d\n", object$stats$n_obs))
  cat(sprintf(
    "Grid points  : p = %d (%d success, %d failure)\n",
    object$stats$p,
    object$stats$p - object$stats$n_fail,
    object$stats$n_fail
  ))
  cat("\nParameters:\n")
  cat(sprintf("  Bandwidth (h) : %g\n", object$params$h))
  cat(sprintf("  Degree (m)    : %d\n", object$params$m))
  cat(sprintf("  Quadrature    : %d MC points\n", object$params$N_quad))
  cat("\nResults (estimate):\n")
  print(summary(object$estimate))
  invisible(object)
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
  main = paste0("densityLP  (h = ", x$params$h, ", m = ", x$params$m, ")"),
  ...
) {
  d <- x$stats$d
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

# S3 methods for "cv_density_lp" ------------------------------------------------

#' Print a cv_density_lp object
#' @param x A \code{"cv_density_lp"} object.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.cv_density_lp <- function(x, ...) {
  score <- x$scores[paste0("m=", x$m_hat), as.character(round(x$h_hat, 6L))]
  cat("LOO cross-validation -- local polynomial density\n")
  cat(sprintf("  Selected m : %d\n", x$m_hat))
  cat(sprintf("  Selected h : %.6g\n", x$h_hat))
  cat(sprintf("  CV score   : %.5g\n", score))
  cat(sprintf(
    "  Grid       : %d degree(s) x %d bandwidth(s)\n",
    length(x$grids$m),
    length(x$grids$h)
  ))
  invisible(x)
}
