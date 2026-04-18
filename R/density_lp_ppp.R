# spatstat wrapper ---------------------------------------------------------------

#' Local polynomial density estimation for a spatial point pattern
#'
#' Convenience wrapper around [density_lp()] for `ppp` objects (spatstat).
#' The evaluation grid and the domain are derived automatically from
#' `Window(pp)`.  The result is a `spatstat.geom::im` pixel image, fully
#' compatible with the spatstat plotting and analysis ecosystem.
#'
#' @param pp A `ppp` object (spatstat.geom).
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0, default 1).
#' @param N_quad Number of quadrature points per target point (default 500).
#' @param nx Number of pixel columns (default 128).
#' @param ny Number of pixel rows (default 128).
#' @return An `im` object (spatstat.geom). Pixels outside `Window(pp)` are
#'   `NA`.
#' @seealso [density_lp()], [domain_from_owin()]
#' @importFrom spatstat.geom Window inside.owin im unitname
#' @export
density_lp_ppp <- function(pp, h, m = 1L, N_quad = 500L, nx = 128L, ny = 128L) {
  stopifnot(inherits(pp, "ppp"))
  stopifnot(is.numeric(h), length(h) == 1L, h > 0)
  stopifnot(is.numeric(m), length(m) == 1L, m >= 0, m == floor(m))
  stopifnot(is.numeric(N_quad), length(N_quad) == 1L, N_quad >= 1L)
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

  spatstat.geom::im(
    t(matrix(values, nrow = nx, ncol = ny)),
    xcol = grid_x,
    yrow = grid_y,
    unitname = spatstat.geom::unitname(pp)
  )
}
