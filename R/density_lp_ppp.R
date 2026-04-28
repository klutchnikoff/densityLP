# spatstat wrapper ---------------------------------------------------------------

#' Internal constructor for density_lp_ppp objects
#' @keywords internal
new_density_lp_ppp <- function(im_obj, params, stats, call) {
  structure(
    c(
      unclass(im_obj),
      list(
        params = params,
        stats = stats,
        call = call
      )
    ),
    class = c("density_lp_ppp", "im")
  )
}

#' Local polynomial density estimation for a spatial point pattern
#'
#' Convenience wrapper around [density_lp()] for `ppp` objects (spatstat).
#' The evaluation grid and the domain are derived automatically from
#' `Window(pp)`.  The result is a `spatstat.geom::im` pixel image, fully
#' compatible with the spatstat plotting and analysis ecosystem.
#'
#' @param pp A `ppp` object (spatstat.geom).
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0, default 0).
#' @param N_quad Number of quadrature points per target point (default 500).
#' @param nx Number of pixel columns (default 128).
#' @param ny Number of pixel rows (default 128).
#' @return An S3 object of class `"density_lp_ppp"`, which also inherits from
#'   `spatstat.geom::im`. Pixels outside `Window(pp)` are `NA`.
#' @seealso [density_lp()], [domain_from_owin()]
#' @examples
#' set.seed(1)
#' pp <- spatstat.geom::ppp(runif(50), runif(50), window = spatstat.geom::owin())
#' fit <- density_lp_ppp(pp, h = 0.3, m = 1L, N_quad = 100L, nx = 20L, ny = 20L)
#' print(fit)
#'
#' @importFrom spatstat.geom Window inside.owin im unitname
#' @export
density_lp_ppp <- function(pp, h, m = 0L, N_quad = 500L, nx = 128L, ny = 128L) {
  stopifnot(inherits(pp, "ppp"))
  check_h(h)
  check_m(m)
  check_N_quad(N_quad)
  stopifnot(is.numeric(nx), nx >= 2L)
  stopifnot(is.numeric(ny), ny >= 2L)

  win <- spatstat.geom::Window(pp)
  dom <- domain_from_owin(win)

  # Regular grid over the bounding box
  grid_x <- seq(win$xrange[1L], win$xrange[2L], length.out = nx)
  grid_y <- seq(win$yrange[1L], win$yrange[2L], length.out = ny)
  t_all <- as.matrix(expand.grid(x = grid_x, y = grid_y))

  # Keep only grid points inside the window
  inside <- spatstat.geom::inside.owin(t_all[, 1L], t_all[, 2L], w = win)
  t_grid <- t_all[inside, , drop = FALSE]

  if (nrow(t_grid) == 0L) {
    stop("No grid points fall inside Window(pp). Try increasing nx or ny.")
  }

  mc <- match.call()

  X <- cbind(pp$x, pp$y)
  fit <- density_lp(
    X,
    t_grid = t_grid,
    h = h,
    m = m,
    domain = dom,
    N_quad = N_quad
  )

  # Fill full nx*ny grid (NA outside window), then reshape for im
  # expand.grid fills x fastest: point (grid_x[j], grid_y[i]) has index (i-1)*nx+j
  # im(mat, xcol, yrow) expects mat[i,j] = value at (xcol[j], yrow[i])
  # => mat_im = t(matrix(values, nrow = nx, ncol = ny))
  values <- rep(NA_real_, nx * ny)
  values[inside] <- fit$estimate

  im_obj <- spatstat.geom::im(
    t(matrix(values, nrow = nx, ncol = ny)),
    xcol = grid_x,
    yrow = grid_y,
    unitname = spatstat.geom::unitname(pp)
  )

  new_density_lp_ppp(
    im_obj = im_obj,
    params = list(
      h = h,
      m = as.integer(m),
      N_quad = as.integer(N_quad),
      domain = dom
    ),
    stats = fit$stats,
    call = mc
  )
}
