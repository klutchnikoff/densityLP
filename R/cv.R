# LOO cross-validation for bandwidth and degree selection -----------------------

#' LOO cross-validation for the local polynomial density estimator
#'
#' For each pair \eqn{(m, h)} in the grid, computes the leave-one-out
#' log-likelihood score
#' \deqn{\mathrm{CV}(m, h) = -\frac{1}{n}\sum_{i=1}^n
#'       \log \hat{f}_{\gamma,-i}(X_i)}
#' using the algebraic LOO correction
#' \deqn{\hat{f}_{\gamma,-i}(X_i) =
#'       \frac{n\,\hat{f}_\gamma(X_i) - h^{-d}\|H_0^{(X_i)}\|^2}{n-1},}
#' where \eqn{\|H_0^{(X_i)}\|^2 = (B_{X_i}^{-1})_{11}} is the self-influence
#' scalar returned as a byproduct of the C++ pipeline.
#'
#' @param X Numeric matrix `n x d` of observations.
#' @param h_grid Numeric vector of bandwidth values (all > 0).
#' @param m_grid Integer vector of polynomial degrees (default `0:3`).
#' @param domain An `"domain_lp"` object.
#' @param N_quad Number of quadrature points per evaluation (default 500).
#' @return An S3 object of class `"cv_density_lp"` containing:
#'   - `$m_hat`, `$h_hat`: selected degree and bandwidth,
#'   - `$scores`: matrix `|m_grid| x |h_grid|` of CV scores,
#'   - `$m_grid`, `$h_grid`: the grids used.
#' @seealso [density_lp()], [cv_density_lp_ppp()]
#' @export
cv_density_lp <- function(X, h_grid, m_grid = 0:3, domain, N_quad = 500L) {
  check_X(X)
  check_domain(domain, ncol(X))
  check_h_grid(h_grid)
  check_m_grid(m_grid)
  check_N_quad(N_quad)

  m_grid <- as.integer(m_grid)
  n <- nrow(X)
  d <- ncol(X)

  scores <- matrix(
    NA_real_,
    length(m_grid),
    length(h_grid),
    dimnames = list(paste0("m=", m_grid), as.character(round(h_grid, 6L)))
  )
  n_fail <- matrix(0L, length(m_grid), length(h_grid))

  sampler <- domain$sampler_factory()
  # We transpose X once for C++ (d x n)
  Xt <- t(X)

  for (ih in seq_along(h_grid)) {
    h <- h_grid[ih]

    # Draw quadrature points for each observation i (only once per h)
    s_list <- vector("list", n)
    n_total_vec <- integer(n)
    for (i in seq_len(n)) {
      t_i <- X[i, , drop = TRUE]
      s <- tryCatch(sampler(N_quad, t_i, h), error = function(e) NULL)
      if (is.null(s)) {
        s_list[[i]] <- matrix(0, d, 0)
        n_total_vec[i] <- 0L
      } else {
        s_list[[i]] <- s$points
        n_total_vec[i] <- s$n_total
      }
    }

    for (im in seq_along(m_grid)) {
      m <- m_grid[im]
      alphas <- build_alphas(m, d)

      # Call the optimized C++ backend for the full observation loop
      log_loo <- cv_lp_fixed_h_m_cpp(Xt, s_list, n_total_vec, h, alphas)

      # Track failures (infinite scores)
      n_fail[im, ih] <- sum(!is.finite(log_loo))
      scores[im, ih] <- -mean(log_loo)
    }
  }

  idx <- which(scores == min(scores, na.rm = TRUE), arr.ind = TRUE)[1L, ]

  fail_rate <- n_fail[idx[1L], idx[2L]] / n
  if (fail_rate >= 0.25) {
    warning(sprintf(
      paste0(
        "%d/%d LOO evaluations failed for the selected parameters",
        " (m = %d, h = %g).",
        " The CV score is unreliable: consider a larger h or a coarser grid."
      ),
      n_fail[idx[1L], idx[2L]],
      n,
      m_grid[idx[1L]],
      h_grid[idx[2L]]
    ))
  } else if (fail_rate >= 0.10) {
    warning(sprintf(
      paste0(
        "%d/%d LOO evaluations failed for the selected parameters",
        " (m = %d, h = %g).",
        " Results should be interpreted with caution."
      ),
      n_fail[idx[1L], idx[2L]],
      n,
      m_grid[idx[1L]],
      h_grid[idx[2L]]
    ))
  }
  structure(
    list(
      m_hat = m_grid[idx[1L]],
      h_hat = h_grid[idx[2L]],
      scores = scores,
      m_grid = m_grid,
      h_grid = h_grid,
      call = match.call()
    ),
    class = "cv_density_lp"
  )
}

#' LOO cross-validation for a spatial point pattern
#'
#' Convenience wrapper around [cv_density_lp()] for `ppp` objects. The domain
#' is derived automatically from `Window(pp)`.
#'
#' @param pp A `ppp` object (spatstat.geom).
#' @param h_grid Numeric vector of bandwidth values (all > 0).
#' @param m_grid Integer vector of polynomial degrees (default `0:3`).
#' @param N_quad Number of quadrature points per evaluation (default 500).
#' @return An S3 object of class `"cv_density_lp"` (see [cv_density_lp()]).
#' @seealso [cv_density_lp()], [density_lp_ppp()]
#' @importFrom spatstat.geom Window
#' @export
cv_density_lp_ppp <- function(pp, h_grid, m_grid = 0:3, N_quad = 500L) {
  stopifnot(inherits(pp, "ppp"))
  X <- cbind(pp$x, pp$y)
  dom <- domain_from_owin(spatstat.geom::Window(pp))
  cv_density_lp(
    X,
    h_grid = h_grid,
    m_grid = m_grid,
    domain = dom,
    N_quad = N_quad
  )
}

#' Print a cv_density_lp object
#' @param x A `"cv_density_lp"` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.cv_density_lp <- function(x, ...) {
  score <- x$scores[paste0("m=", x$m_hat), as.character(round(x$h_hat, 6L))]
  cat("LOO cross-validation -- local polynomial density\n")
  cat(sprintf("  Selected m : %d\n", x$m_hat))
  cat(sprintf("  Selected h : %.6g\n", x$h_hat))
  cat(sprintf("  CV score   : %.5g\n", score))
  cat(sprintf(
    "  Grid       : %d degree(s) x %d bandwidth(s)\n",
    length(x$m_grid),
    length(x$h_grid)
  ))
  invisible(x)
}
