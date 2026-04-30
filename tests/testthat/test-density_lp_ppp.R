# ── density_lp.ppp ────────────────────────────────────────────────────────────

test_that("density_lp.ppp: returns density_lp and im classes", {
  skip_if_not_installed("spatstat.geom")
  set.seed(1L)
  r <- sqrt(runif(80L))
  theta <- runif(80L) * 2 * pi
  pp <- spatstat.geom::ppp(
    r * cos(theta),
    r * sin(theta),
    window = spatstat.geom::disc()
  )
  fit <- density_lp(pp, h = 0.4, m = 1L, nx = 32L, ny = 32L)
  expect_s3_class(fit, "density_lp")
  expect_s3_class(fit, "im")
})

test_that("density_lp.ppp: print method output matches 2D", {
  skip_if_not_installed("spatstat.geom")
  set.seed(2L)
  pp <- spatstat.geom::ppp(
    runif(100L),
    runif(100L),
    window = spatstat.geom::owin()
  )
  fit <- density_lp(pp, h = 0.5, nx = 32L, ny = 32L)
  expect_output(print(fit), "spatstat-compatible")
  expect_identical(withVisible(print(fit))$visible, FALSE)
})

test_that("density_lp.ppp: respects nx and ny", {
  skip_if_not_installed("spatstat.geom")
  set.seed(3L)
  pp <- spatstat.geom::ppp(
    runif(50L),
    runif(50L),
    window = spatstat.geom::owin()
  )
  nx <- 20L
  ny <- 25L
  fit <- density_lp(pp, h = 0.4, nx = nx, ny = ny)
  expect_equal(length(fit$xcol), nx)
  expect_equal(length(fit$yrow), ny)
  expect_equal(dim(fit$v), c(ny, nx))
})

test_that("density_lp.ppp: error on non-positive bandwidth", {
  skip_if_not_installed("spatstat.geom")
  pp <- spatstat.geom::ppp(0.5, 0.5, window = spatstat.geom::owin())
  expect_error(density_lp(pp, h = -0.1), "positive")
})

test_that("density_lp.ppp: m=0 behaves correctly", {
  skip_if_not_installed("spatstat.geom")
  set.seed(4L)
  pp <- spatstat.geom::ppp(
    runif(100L),
    runif(100L),
    window = spatstat.geom::owin()
  )
  fit <- density_lp(pp, h = 0.3, m = 0L, nx = 20L, ny = 20L)
  expect_true(all(fit$v >= 0, na.rm = TRUE))
})

test_that("density_lp.ppp: NA outside the window", {
  skip_if_not_installed("spatstat.geom")
  # Circular window
  win <- spatstat.geom::disc(radius = 0.5, centre = c(0.5, 0.5))
  pp <- spatstat.geom::ppp(0.5, 0.5, window = win)
  fit <- density_lp(pp, h = 0.2, nx = 32L, ny = 32L)
  # Corner (0,0) is outside the disc
  expect_true(is.na(fit$v[1L, 1L]))
})

test_that("density_lp_ppp alias: error if input is not a ppp", {
  expect_error(density_lp_ppp(matrix(1, 10, 2), h = 0.3), "ppp")
})

test_that("density_lp.ppp: integration on unit square close to 1", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(101L)
  pp <- spatstat.random::runifpoint(1000L, win = spatstat.geom::owin())
  fit <- density_lp(pp, h = 0.3, nx = 50L, ny = 50L)
  # Area of pixel
  area <- (fit$xcol[2L] - fit$xcol[1L]) * (fit$yrow[2L] - fit$yrow[1L])
  int_val <- sum(fit$v, na.rm = TRUE) * area
  expect_gt(int_val, 0.5)
  expect_lt(int_val, 1.5)
})
