# LOO cross-validation for bandwidth and degree selection -----------------------

#' Internal constructor for cv_density_lp objects
#' @keywords internal
new_cv_density_lp <- function(m_hat, h_hat, scores, grids, call) {
  structure(
    list(
      m_hat = m_hat,
      h_hat = h_hat,
      scores = scores,
      grids = grids,
      call = call
    ),
    class = "cv_density_lp"
  )
}

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
#' @param x Observations. Can be a numeric \code{matrix}, a \code{spatstat.geom::ppp}
#'   object, or a \code{pp_lp} object.
#' @param h_grid Numeric vector of bandwidth values (all > 0).
#' @param m_grid Integer vector of polynomial degrees (default \code{0:3}).
#' @param N_quad Number of quadrature (Monte Carlo) points used to approximate
#'   the normalising integral (default 500).
#' @param ... Additional arguments passed to methods.
#'
#' @return An S3 object of class \code{"cv_density_lp"} with fields:
#'   \describe{
#'     \item{m_hat}{Selected polynomial degree.}
#'     \item{h_hat}{Selected bandwidth.}
#'     \item{scores}{Matrix of CV scores (rows = degrees, columns = bandwidths).}
#'     \item{grids}{List with \code{$m} and \code{$h} grids used.}
#'   }
#' @seealso [density_lp()]
#' @examples
#' set.seed(1)
#' X <- matrix(runif(30 * 2), ncol = 2)
#' dom <- domain_Rd(2L)
#' cv <- cv_density_lp(X, h_grid = c(0.2, 0.4), m_grid = 0:1, domain = dom, N_quad = 100L)
#' print(cv)
#' summary(cv)
#' @export
cv_density_lp <- function(x, h_grid, m_grid = 0:3, N_quad = 500L, ...) {
  UseMethod("cv_density_lp")
}

#' @export
cv_density_lp.matrix <- function(
  x,
  h_grid,
  m_grid = 0:3,
  N_quad = 500L,
  domain,
  ...
) {
  pp <- pp_lp(x, domain)
  cv_density_lp(pp, h_grid = h_grid, m_grid = m_grid, N_quad = N_quad, ...)
}

#' @export
cv_density_lp.data.frame <- function(
  x,
  h_grid,
  m_grid = 0:3,
  N_quad = 500L,
  domain,
  ...
) {
  pp <- pp_lp(as.matrix(x), domain)
  cv_density_lp(pp, h_grid = h_grid, m_grid = m_grid, N_quad = N_quad, ...)
}

#' @export
cv_density_lp.ppp <- function(x, h_grid, m_grid = 0:3, N_quad = 500L, ...) {
  pp <- as_pp_lp(x)
  cv_density_lp(pp, h_grid = h_grid, m_grid = m_grid, N_quad = N_quad, ...)
}

#' @export
cv_density_lp.pp_lp <- function(x, h_grid, m_grid = 0:3, N_quad = 500L, ...) {
  check_h_grid(h_grid)
  check_m_grid(m_grid)
  check_N_quad(N_quad)

  m_grid <- as.integer(m_grid)
  n <- x$n
  d <- x$d

  scores <- matrix(
    NA_real_,
    length(m_grid),
    length(h_grid),
    dimnames = list(paste0("m=", m_grid), as.character(round(h_grid, 6L)))
  )
  n_fail <- matrix(0L, length(m_grid), length(h_grid))

  sampler <- x$domain$sampler_factory()
  Xt <- t(x$X)

  for (ih in seq_along(h_grid)) {
    h <- h_grid[ih]

    # Draw quadrature points for each observation i (only once per h)
    s_list <- vector("list", n)
    n_total_vec <- integer(n)
    for (i in seq_len(n)) {
      t_i <- x$X[i, , drop = TRUE]
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

      log_loo <- cv_lp_fixed_h_m_cpp(Xt, s_list, n_total_vec, h, alphas)

      n_fail[im, ih] <- sum(!is.finite(log_loo))
      scores[im, ih] <- -mean(log_loo)
    }
  }

  idx <- arrayInd(which.min(scores), dim(scores))

  fail_rate <- n_fail[idx[1L], idx[2L]] / n
  if (fail_rate >= 0.25) {
    warning(sprintf(
      paste0(
        "%d/%d LOO evaluations failed for the selected parameters",
        " (m = %d, h = %g): V(h) empty or Gram matrix singular.",
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
        " (m = %d, h = %g): V(h) empty or Gram matrix singular.",
        " Results should be interpreted with caution."
      ),
      n_fail[idx[1L], idx[2L]],
      n,
      m_grid[idx[1L]],
      h_grid[idx[2L]]
    ))
  }
  new_cv_density_lp(
    m_hat = m_grid[idx[1L]],
    h_hat = h_grid[idx[2L]],
    scores = scores,
    grids = list(
      m = m_grid,
      h = h_grid
    ),
    call = match.call()
  )
}

#' LOO cross-validation for a spatial point pattern
#'
#' This function is now an alias for [cv_density_lp()] when called on a \code{ppp}
#' object.
#'
#' @param pp A `ppp` object (spatstat.geom).
#' @param h_grid Numeric vector of bandwidth values (all > 0).
#' @param m_grid Integer vector of polynomial degrees (default `0:3`).
#' @param N_quad Number of quadrature points per evaluation (default 500).
#' @param ... Additional arguments passed to [cv_density_lp()].
#' @return An S3 object of class \code{"cv_density_lp"} (see [cv_density_lp()]).
#' @examples
#' \donttest{
#' set.seed(1)
#' pp <- spatstat.geom::ppp(runif(50), runif(50), window = spatstat.geom::owin())
#' cv <- cv_density_lp_ppp(pp, h_grid = c(0.2, 0.4), m_grid = 0:1, N_quad = 100L)
#' print(cv)
#' }
#' @export
cv_density_lp_ppp <- function(pp, h_grid, m_grid = 0:3, N_quad = 500L, ...) {
  if (!inherits(pp, "ppp")) {
    stop("cv_density_lp_ppp requires a 'ppp' object from spatstat.geom")
  }
  cv_density_lp(pp, h_grid = h_grid, m_grid = m_grid, N_quad = N_quad, ...)
}
