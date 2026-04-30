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
#' This function is now an alias for [density_lp()] when called on a `ppp`
#' object.
#'
#' @param pp A `ppp` object (spatstat.geom).
#' @param h Bandwidth (scalar > 0).
#' @param m Polynomial degree (integer >= 0, default 0).
#' @param N_quad Number of quadrature points per target point (default 500).
#' @param nx Number of pixel columns (default 128).
#' @param ny Number of pixel rows (default 128).
#' @param ... Additional arguments passed to [density_lp()].
#' @return An S3 object of class \code{"density_lp"}, which also inherits from
#'   \code{spatstat.geom::im}. Pixels outside \code{Window(pp)} are \code{NA}.
#' @seealso [density_lp()]
#' @examples
#' \donttest{
#' set.seed(1)
#' pp <- spatstat.geom::ppp(runif(50), runif(50), window = spatstat.geom::owin())
#' fit <- density_lp_ppp(pp, h = 0.3, m = 1L, N_quad = 100L, nx = 20L, ny = 20L)
#' print(fit)
#' }
#' @export
density_lp_ppp <- function(
  pp,
  h,
  m = 0L,
  N_quad = 500L,
  nx = 128L,
  ny = 128L,
  ...
) {
  if (!inherits(pp, "ppp")) {
    stop("density_lp_ppp requires a 'ppp' object from spatstat.geom")
  }
  density_lp(pp, h = h, m = m, N_quad = N_quad, nx = nx, ny = ny, ...)
}
