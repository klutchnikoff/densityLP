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
#' @param domain An `"lp_domain"` object.
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
  stopifnot(is.numeric(h_grid), length(h_grid) >= 1L, all(h_grid > 0))
  stopifnot(
    is.numeric(m_grid),
    length(m_grid) >= 1L,
    all(m_grid >= 0),
    all(m_grid == floor(m_grid))
  )
  stopifnot(is.numeric(N_quad), length(N_quad) == 1L, N_quad >= 1L)

  m_grid <- as.integer(m_grid)
  n <- nrow(X)
  d <- ncol(X)

  scores <- matrix(
    NA_real_,
    length(m_grid),
    length(h_grid),
    dimnames = list(paste0("m=", m_grid), as.character(round(h_grid, 6L)))
  )

  sampler <- domain$sampler_factory()

  for (im in seq_along(m_grid)) {
    m <- m_grid[im]
    alphas <- build_alphas(m, d)

    for (ih in seq_along(h_grid)) {
      h <- h_grid[ih]
      log_loo <- numeric(n)

      for (i in seq_len(n)) {
        t_i <- X[i, , drop = TRUE]
        s <- tryCatch(sampler(N_quad, t_i, h), error = function(e) NULL)

        if (is.null(s) || ncol(s$points) < 2L) {
          log_loo[i] <- -Inf
          next
        }

        U_obs <- t(sweep(X, 2L, t_i, "-"))
        res <- lp_estimator_loo_cpp(s$points, s$n_total, U_obs, h, alphas)

        f_hat_i <- res[1L]
        norm_H0_sq <- res[2L]
        f_loo_i <- (n * f_hat_i - norm_H0_sq / h^d) / (n - 1L)

        log_loo[i] <- if (f_loo_i > 0) log(f_loo_i) else -Inf
      }

      scores[im, ih] <- -mean(log_loo)
    }
  }

  idx <- which(scores == min(scores, na.rm = TRUE), arr.ind = TRUE)[1L, ]
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
  cat("LOO cross-validation -- local polynomial density\n")
  cat("  Selected m :", x$m_hat, "\n")
  cat("  Selected h :", round(x$h_hat, 6L), "\n")
  cat(
    "  CV score   :",
    round(
      x$scores[paste0("m=", x$m_hat), as.character(round(x$h_hat, 6L))],
      5L
    ),
    "\n"
  )
  cat(
    "  Grid       :",
    length(x$m_grid),
    "degree(s) x",
    length(x$h_grid),
    "bandwidth(s)\n"
  )
  invisible(x)
}
