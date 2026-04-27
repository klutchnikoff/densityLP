# ── density_lp_ppp ────────────────────────────────────────────────────────────

test_that("density_lp_ppp: returns a density_lp_ppp object inheriting from im", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(1L)
  pp <- spatstat.random::runifpoint(100L, win = spatstat.geom::disc())
  fit <- density_lp_ppp(pp, h = 0.5, m = 0L, N_quad = 100L, nx = 32L, ny = 32L)
  expect_s3_class(fit, "density_lp_ppp")
  expect_s3_class(fit, "im")
  expect_equal(fit$h, 0.5)
  expect_equal(fit$m, 0L)
})

test_that("print.density_lp_ppp: returns x invisibly and produces output", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(7L)
  pp <- spatstat.random::runifpoint(100L, win = spatstat.geom::disc())
  fit <- density_lp_ppp(pp, h = 0.5, m = 0L, N_quad = 100L, nx = 32L, ny = 32L)
  expect_output(print(fit), "spatstat")
  expect_identical(withVisible(print(fit))$visible, FALSE)
})

test_that("density_lp_ppp: pixel matrix dimensions match nx and ny", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(2L)
  nx <- 20L
  ny <- 25L
  pp <- spatstat.random::runifpoint(100L, win = spatstat.geom::disc())
  fit <- density_lp_ppp(pp, h = 0.5, m = 0L, N_quad = 100L, nx = nx, ny = ny)
  expect_equal(dim(fit$v), c(ny, nx))
})

test_that("density_lp_ppp: pixels outside window are NA", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(3L)
  pp <- spatstat.random::runifpoint(100L, win = spatstat.geom::disc())
  fit <- density_lp_ppp(pp, h = 0.5, m = 0L, N_quad = 100L, nx = 32L, ny = 32L)
  # Disc bounding box [-1,1]^2: corner pixels are outside the disc
  expect_true(is.na(fit$v[1L, 1L]))
  expect_true(any(is.na(fit$v)))
})

test_that("density_lp_ppp: non-NA estimates are non-negative (m = 0)", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(4L)
  pp <- spatstat.random::runifpoint(150L, win = spatstat.geom::disc())
  fit <- density_lp_ppp(pp, h = 0.5, m = 0L, N_quad = 150L, nx = 32L, ny = 32L)
  expect_true(all(fit$v[!is.na(fit$v)] >= 0))
})

test_that("density_lp_ppp: integral is approximately 1 for uniform data", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(5L)
  pp <- spatstat.random::runifpoint(300L, win = spatstat.geom::disc())
  fit <- density_lp_ppp(pp, h = 0.5, m = 0L, N_quad = 200L, nx = 64L, ny = 64L)
  int_val <- sum(fit$v, na.rm = TRUE) * fit$xstep * fit$ystep
  expect_gt(int_val, 0.5)
  expect_lt(int_val, 1.8)
})

test_that("density_lp_ppp: error if input is not a ppp", {
  expect_error(density_lp_ppp(matrix(1, 10, 2), h = 0.3), "ppp")
})

test_that("density_lp_ppp: error if h <= 0", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set.seed(6L)
  pp <- spatstat.random::runifpoint(50L, win = spatstat.geom::disc())
  expect_error(density_lp_ppp(pp, h = -0.1))
})
