# Grid estimator (exported) -----------------------------------------------------

#' Internal constructor for density_lp objects
#' @keywords internal
new_density_lp <- function(
  estimate,
  variance,
  t_grid,
  X,
  h,
  m,
  domain,
  N_quad,
  n_fail,
  call
) {
  structure(
    list(
      estimate = estimate,
      variance = variance,
      t_grid = t_grid,
      X = X,
      params = list(
        h = h,
        m = m,
        domain = domain,
        N_quad = N_quad
      ),
      stats = list(
        n_fail = n_fail,
        n_obs = nrow(X),
        p = nrow(t_grid),
        d = ncol(X)
      ),
      call = call
    ),
    class = "density_lp"
  )
}

#' Local polynomial density estimation
#'
#' @param x Observations. Can be a numeric \code{matrix}, a \code{spatstat.geom::ppp}
#'   object, or a \code{pp_lp} object.
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0, default 0).
#' @param t_grid Numeric matrix \code{p x d} of target points. If \code{NULL},
#'   it is automatically generated from the bounding box of \code{x}.
#' @param N_quad Number of quadrature (Monte Carlo) points used to approximate
#'   the normalising integral (default 500).
#' @param ... Additional arguments passed to methods or for grid generation
#'   (\code{nx}, \code{ny} for 2D).
#'
#' @return An S3 object of class \code{"density_lp"}. In 2D, it also inherits
#'   from \code{spatstat.geom::im}.
#'
#' @examples
#' # 2D estimation on the unit square
#' set.seed(1)
#' X <- matrix(runif(50 * 2), ncol = 2)
#' dom <- domain_Rd(2L)
#' fit <- density_lp(X, h = 0.3, m = 1L, domain = dom, nx = 10L, ny = 10L, N_quad = 100L)
#' print(fit)
#' summary(fit)
#'
#' # 1D estimation with plot
#' set.seed(2)
#' X1 <- matrix(runif(50), ncol = 1)
#' t1 <- matrix(seq(0.1, 0.9, length.out = 20), ncol = 1)
#' fit1 <- density_lp(X1, h = 0.3, domain = domain_Rd(1L), t_grid = t1, N_quad = 100L)
#' plot(fit1)
#'
#' @export
density_lp <- function(x, h, m = 0L, t_grid = NULL, N_quad = 500L, ...) {
  UseMethod("density_lp")
}

#' @export
density_lp.matrix <- function(
  x,
  h,
  m = 0L,
  t_grid = NULL,
  N_quad = 500L,
  domain,
  ...
) {
  pp <- pp_lp(x, domain)
  obj <- density_lp(pp, h = h, m = m, t_grid = t_grid, N_quad = N_quad, ...)

  # If t_grid was provided manually but is a regular grid, we could try to promote to im.
  # For now, let's ensure that if it's already a density_lp from pp_lp, we don't break it.
  obj
}

#' @export
density_lp.data.frame <- function(
  x,
  h,
  m = 0L,
  t_grid = NULL,
  N_quad = 500L,
  domain,
  ...
) {
  pp <- pp_lp(as.matrix(x), domain)
  density_lp(pp, h = h, m = m, t_grid = t_grid, N_quad = N_quad, ...)
}

#' @importFrom spatstat.geom im
#' @export
density_lp.ppp <- function(x, h, m = 0L, t_grid = NULL, N_quad = 500L, ...) {
  pp <- as_pp_lp(x)
  density_lp(pp, h = h, m = m, t_grid = t_grid, N_quad = N_quad, ...)
}

#' @export
density_lp.pp_lp <- function(x, h, m = 0L, t_grid = NULL, N_quad = 500L, ...) {
  check_h(h)
  check_m(m)
  check_N_quad(N_quad)

  if (is.null(t_grid)) {
    t_grid <- .make_auto_grid(x, ...)
  }
  check_t_grid(t_grid, x$d)

  p <- nrow(t_grid)
  estimate <- numeric(p)
  variance <- numeric(p)
  n_fail <- 0L

  for (i in seq_len(p)) {
    res <- tryCatch(
      density_lp_point(
        X = x$X,
        t = t_grid[i, ],
        h = h,
        m = m,
        domain = x$domain,
        N_quad = N_quad
      ),
      error = function(e) NULL
    )
    if (is.null(res)) {
      estimate[i] <- NA_real_
      variance[i] <- NA_real_
      n_fail <- n_fail + 1L
    } else {
      estimate[i] <- res["estimate"]
      variance[i] <- res["variance"]
    }
  }

  if (n_fail > 0L) {
    warning(
      sprintf(
        "%d/%d grid point(s) failed (V(h) empty or Gram matrix singular). ",
        n_fail,
        p
      ),
      "Corresponding estimates are NA. Consider increasing N_quad or reducing h.",
      call. = FALSE
    )
  }

  obj <- new_density_lp(
    estimate = estimate,
    variance = variance,
    t_grid = t_grid,
    X = x$X,
    h = h,
    m = as.integer(m),
    domain = x$domain,
    N_quad = as.integer(N_quad),
    n_fail = n_fail,
    call = match.call()
  )

  # Harmonization for 2D (spatstat compatibility)
  if (x$d == 2L && !is.null(attr(t_grid, "grid_info"))) {
    info <- attr(t_grid, "grid_info")
    values <- rep(NA_real_, info$nx * info$ny)
    values[info$inside] <- obj$estimate
    im_obj <- spatstat.geom::im(
      t(matrix(values, nrow = info$nx, ncol = info$ny)),
      xcol = info$xcol,
      yrow = info$yrow
    )
    obj <- structure(
      c(unclass(obj), unclass(im_obj)),
      class = c("density_lp", "im")
    )
  }

  obj
}

.make_auto_grid <- function(x, nx = 128L, ny = 128L, ...) {
  if (x$d == 2L) {
    xcol <- seq(x$bbox[1, 1], x$bbox[2, 1], length.out = nx)
    yrow <- seq(x$bbox[1, 2], x$bbox[2, 2], length.out = ny)
    grid <- as.matrix(expand.grid(x = xcol, y = yrow))
    inside <- x$domain$is_in_domain(grid[, 1], grid[, 2])
    t_grid <- grid[inside, , drop = FALSE]
    attr(t_grid, "grid_info") <- list(
      nx = nx,
      ny = ny,
      xcol = xcol,
      yrow = yrow,
      inside = inside
    )
    return(t_grid)
  } else {
    # For d > 2, we could implement a default, but it's risky.
    # For now, require t_grid if not 2D, or use a coarse grid.
    stop(
      "Automatic grid generation is currently only supported for d = 2. Please provide t_grid."
    )
  }
}
