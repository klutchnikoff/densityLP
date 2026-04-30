# S3 methods for "density_lp" ---------------------------------------------------

#' @rdname density_lp
#' @export
print.density_lp <- function(x, ...) {
  if (inherits(x, "im")) {
    cat("Local polynomial density estimate (2D, spatstat-compatible)\n")
    cat(sprintf("  Points : %d\n", x$stats$n_obs))
    cat(sprintf("  h      : %g\n", x$params$h))
    cat(sprintf("  m      : %d\n", x$params$m))
    cat(sprintf("  Grid   : %d x %d pixels\n", length(x$xcol), length(x$yrow)))
    rng <- range(x$v, na.rm = TRUE)
    cat(sprintf("  Range  : [%.5g, %.5g]\n", rng[1L], rng[2L]))
  } else {
    cat("Local polynomial density estimate\n")
    cat(sprintf("  d      : %d\n", x$stats$d))
    cat(sprintf("  Points : %d\n", x$stats$p))
    cat(sprintf("  h      : %g\n", x$params$h))
    cat(sprintf("  m      : %d\n", x$params$m))
    rng <- range(x$estimate, na.rm = TRUE)
    cat(sprintf("  Range  : [%.5g, %.5g]\n", rng[1L], rng[2L]))
  }
  if (!is.null(x$stats$n_fail) && x$stats$n_fail > 0L) {
    cat(sprintf(
      "  Note   : %d failure(s) during calculation\n",
      x$stats$n_fail
    ))
  }
  invisible(x)
}

#' @rdname density_lp
#' @param object A \code{"density_lp"} object (for \code{summary}).
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
  if (inherits(object, "im")) {
    cat(sprintf(
      "Grid         : %d x %d pixels (%d evaluated, %d failure)\n",
      length(object$xcol),
      length(object$yrow),
      object$stats$p - object$stats$n_fail,
      object$stats$n_fail
    ))
  } else {
    cat(sprintf(
      "Grid points  : p = %d (%d success, %d failure)\n",
      object$stats$p,
      object$stats$p - object$stats$n_fail,
      object$stats$n_fail
    ))
  }
  cat("\nParameters:\n")
  cat(sprintf("  Bandwidth (h) : %g\n", object$params$h))
  cat(sprintf("  Degree (m)    : %d\n", object$params$m))
  cat(sprintf("  Quadrature    : %d MC points\n", object$params$N_quad))
  cat("\nResults (estimate):\n")
  if (inherits(object, "im")) {
    print(summary(as.vector(object$v), na.rm = TRUE))
  } else {
    print(summary(object$estimate, na.rm = TRUE))
  }
  invisible(object)
}

#' @rdname density_lp
#' @export
plot.density_lp <- function(x, ...) {
  if (inherits(x, "im")) {
    # Dispatch to spatstat.geom::plot.im
    NextMethod()
  } else {
    d <- x$stats$d
    if (d == 1L) {
      plot(
        x$t_grid[, 1L],
        x$estimate,
        type = "l",
        xlab = "t",
        ylab = "Estimated density",
        ...
      )
    } else {
      warning(
        "plot.density_lp: no default graphical method for d = ",
        d,
        ". Using scatter plot."
      )
      plot(x$t_grid, col = "blue", pch = 16, ...)
    }
  }
  invisible(x)
}

# Alias methods for backward compatibility --------------------------------------

#' @export
print.density_lp_ppp <- function(x, ...) print.density_lp(x, ...)

#' @export
summary.density_lp_ppp <- function(object, ...) summary.density_lp(object, ...)

#' @export
plot.density_lp_ppp <- function(x, ...) plot.density_lp(x, ...)


# S3 methods for "cv_density_lp" ------------------------------------------------

#' @rdname cv_density_lp
#' @param object A \code{"cv_density_lp"} object (for \code{summary}).
#' @export
summary.cv_density_lp <- function(object, ...) {
  row_idx <- which(object$grids$m == object$m_hat)
  col_idx <- which(object$grids$h == object$h_hat)
  score <- object$scores[row_idx, col_idx]
  cat("LOO Cross-Validation Summary -- local polynomial density\n")
  cat("---------------------------------------------------------\n")
  cat(sprintf(
    "Call      : %s\n",
    paste(deparse(object$call), collapse = "\n")
  ))
  cat(sprintf("m grid    : %s\n", paste(object$grids$m, collapse = ", ")))
  cat(sprintf(
    "h grid    : %s\n",
    paste(round(object$grids$h, 4L), collapse = ", ")
  ))
  cat("\nCV score matrix (lower is better):\n")
  print(round(object$scores, 4L))
  cat(sprintf("\nSelected  : m = %d, h = %.6g\n", object$m_hat, object$h_hat))
  cat(sprintf("CV score  : %.5g\n", score))
  invisible(object)
}

#' @rdname cv_density_lp
#' @export
print.cv_density_lp <- function(x, ...) {
  row_idx <- which(x$grids$m == x$m_hat)
  col_idx <- which(x$grids$h == x$h_hat)
  score <- x$scores[row_idx, col_idx]
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
