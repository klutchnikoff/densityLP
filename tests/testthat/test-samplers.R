# ── sampler_qmc ────────────────────────────────────────────────────────────────

test_that("sampler_qmc: returns N points of correct dimension", {
  set.seed(1L)
  is_in <- function(X) rep(TRUE, nrow(X))
  s <- sampler_qmc(is_in)(100L, t = c(0.5, 0.5), h = 0.3)
  expect_equal(nrow(s$points), 2L)
  expect_equal(ncol(s$points), 100L)
  expect_true(s$vol > 0)
})

test_that("sampler_qmc: all points lie in [-h, h]^d", {
  set.seed(2L)
  h <- 0.4
  is_in <- function(X) rep(TRUE, nrow(X))
  s <- sampler_qmc(is_in)(200L, t = c(0.5, 0.5), h = h)
  expect_true(all(abs(s$points) <= h + 1e-9))
})

test_that("sampler_qmc: volume estimate close to (2h)^d for full box", {
  set.seed(3L)
  d <- 2L
  h <- 0.3
  is_in <- function(X) rep(TRUE, nrow(X))
  s <- sampler_qmc(is_in)(500L, t = rep(0.5, d), h = h)
  expect_equal(s$vol, (2 * h)^d, tolerance = 1e-10)
})

test_that("sampler_qmc: points respect domain constraint", {
  set.seed(4L)
  is_in <- function(X) rowSums(X^2) <= 1
  s <- sampler_qmc(is_in)(200L, t = c(0, 0), h = 0.5)
  pts_global <- sweep(t(s$points), 2L, c(0, 0), "+")
  expect_true(all(rowSums(pts_global^2) <= 1 + 1e-9))
})

test_that("sampler_qmc: lower variance than rejection for smooth domain", {
  # Estimate vol(disk of radius 0.3) ~ pi * 0.3^2 = 0.2827
  # QMC should have smaller variance over 20 replications
  is_in <- function(X) rowSums(X^2) <= 0.09
  N <- 300L
  t_pt <- c(0, 0)
  h <- 0.4
  R <- 20L

  set.seed(5L)
  vols_rej <- replicate(R, sampler_rejection(is_in)(N, t_pt, h)$vol)
  set.seed(6L)
  vols_qmc <- replicate(R, sampler_qmc(is_in)(N, t_pt, h)$vol)

  expect_lt(var(vols_qmc), var(vols_rej))
})

# ── domain_func with method = "qmc" ───────────────────────────────────────────

test_that("domain_func(method='qmc'): correct label", {
  dom <- domain_func(function(X) rep(TRUE, nrow(X)), d = 2L, method = "qmc")
  expect_match(dom$label, "qmc")
})

test_that("domain_func(method='qmc'): produces valid samples", {
  set.seed(10L)
  dom <- domain_func(
    function(X) rowSums(X^2) <= 1,
    d = 2L,
    method = "qmc"
  )
  s <- dom$sampler_factory()(100L, t = c(0, 0), h = 0.5)
  expect_equal(nrow(s$points), 2L)
  expect_equal(ncol(s$points), 100L)
  expect_true(s$vol > 0)
})

# ── sampler_spatstat ───────────────────────────────────────────────────────────

test_that("sampler_spatstat: returns N points of dimension 2", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  s <- sampler_spatstat(win)(100L, t = c(0.5, 0.5), h = 0.3)
  expect_equal(nrow(s$points), 2L)
  expect_equal(ncol(s$points), 100L)
})

test_that("sampler_spatstat: volume is exact (area of intersection)", {
  skip_if_not_installed("spatstat.geom")
  set.seed(21L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  h <- 0.3
  t_pt <- c(0.5, 0.5)
  s <- sampler_spatstat(win)(200L, t = t_pt, h = h)
  # Full box [0.2, 0.8]^2 is inside [0,1]^2 -> vol = (2h)^2
  expect_equal(s$vol, (2 * h)^2, tolerance = 1e-10)
})

test_that("sampler_spatstat: volume correct at boundary (clipped box)", {
  skip_if_not_installed("spatstat.geom")
  set.seed(22L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  h <- 0.3
  t_pt <- c(0.1, 0.5) # box clips at x=0
  s <- sampler_spatstat(win)(200L, t = t_pt, h = h)
  # x range: [max(0, -0.2), 0.4] = [0, 0.4], width = 0.4 < 2h = 0.6
  expect_lt(s$vol, (2 * h)^2)
  expect_equal(s$vol, 0.4 * 0.6, tolerance = 1e-10)
})

test_that("sampler_spatstat: points lie within the owin", {
  skip_if_not_installed("spatstat.geom")
  set.seed(23L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  t_pt <- c(0.5, 0.5)
  h <- 0.4
  s <- sampler_spatstat(win)(150L, t = t_pt, h = h)
  pts_x <- s$points[1L, ] + t_pt[1L]
  pts_y <- s$points[2L, ] + t_pt[2L]
  expect_true(all(pts_x >= 0 - 1e-9 & pts_x <= 1 + 1e-9))
  expect_true(all(pts_y >= 0 - 1e-9 & pts_y <= 1 + 1e-9))
})
