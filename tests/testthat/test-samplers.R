# ── sampler_rejection ──────────────────────────────────────────────────────────

test_that("sampler_rejection: points have correct dimension", {
  set.seed(1L)
  is_in <- function(X) rep(TRUE, nrow(X))
  s <- sampler_rejection(is_in)(100L, t = c(0.5, 0.5), h = 0.3)
  expect_equal(nrow(s$points), 2L)
  expect_true(ncol(s$points) > 0L)
})

test_that("sampler_rejection: n_total equals N_quad", {
  set.seed(2L)
  is_in <- function(X) rep(TRUE, nrow(X))
  s <- sampler_rejection(is_in)(200L, t = c(0.5, 0.5), h = 0.3)
  expect_equal(s$n_total, 200L)
})

test_that("sampler_rejection: all points lie in [-h, h]^d", {
  set.seed(3L)
  h <- 0.4
  is_in <- function(X) rep(TRUE, nrow(X))
  s <- sampler_rejection(is_in)(300L, t = c(0.5, 0.5), h = h)
  expect_true(all(abs(s$points) <= h + 1e-9))
})

test_that("sampler_rejection: points respect domain constraint", {
  set.seed(4L)
  is_in <- function(X) rowSums(X^2) <= 1
  s <- sampler_rejection(is_in)(300L, t = c(0, 0), h = 0.5)
  pts_global <- sweep(t(s$points), 2L, c(0, 0), "+")
  expect_true(all(rowSums(pts_global^2) <= 1 + 1e-9))
})

test_that("sampler_rejection: error when V(h) is empty", {
  is_in <- function(X) rep(FALSE, nrow(X))
  expect_error(sampler_rejection(is_in)(100L, t = c(0.5, 0.5), h = 0.3))
})

# ── sampler_sector ─────────────────────────────────────────────────────────────

test_that("sampler_sector: returns exactly N_quad points of dimension 2", {
  set.seed(10L)
  s <- sampler_sector(2)(300L, t = c(0.5, 0.1), h = 0.2)
  expect_equal(nrow(s$points), 2L)
  expect_equal(ncol(s$points), 300L)
})

test_that("sampler_sector: points lie within D_k", {
  set.seed(11L)
  k <- 2
  t_pt <- c(0.5, 0.1)
  s <- sampler_sector(k)(300L, t = t_pt, h = 0.2)
  pts <- sweep(t(s$points), 2L, t_pt, "+")
  x <- pts[, 1L]
  y <- pts[, 2L]
  expect_true(all(x >= 0 - 1e-9 & x <= 1 + 1e-9))
  expect_true(all(y >= 0 - 1e-9 & y <= x^k + 1e-9))
})

test_that("sampler_sector: n_total > N_quad near apex (small volume)", {
  set.seed(12L)
  N <- 200L
  s <- sampler_sector(2)(N, t = c(0.05, 0.01), h = 0.1)
  expect_gt(s$n_total, N)
})

test_that("sampler_sector: error when x range is empty", {
  expect_error(sampler_sector(2)(100L, t = c(-0.5, 0), h = 0.3))
})

# ── sampler_owin ───────────────────────────────────────────────────────────

test_that("sampler_owin: returns N points of dimension 2", {
  skip_if_not_installed("spatstat.geom")
  set.seed(20L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  s <- sampler_owin(win)(100L, t = c(0.5, 0.5), h = 0.3)
  expect_equal(nrow(s$points), 2L)
  expect_equal(ncol(s$points), 100L)
})

test_that("sampler_owin: n_total equals N_quad for full box", {
  skip_if_not_installed("spatstat.geom")
  set.seed(21L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  h <- 0.3
  t_pt <- c(0.5, 0.5)
  s <- sampler_owin(win)(200L, t = t_pt, h = h)
  # Full box [0.2, 0.8]^2 is inside [0,1]^2 -> vol = (2h)^2 -> n_total = N
  expect_equal(s$n_total, 200L)
})

test_that("sampler_owin: n_total > N_quad at boundary (clipped box)", {
  skip_if_not_installed("spatstat.geom")
  set.seed(22L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  h <- 0.3
  t_pt <- c(0.1, 0.5) # box clips at x=0
  N <- 200L
  s <- sampler_owin(win)(N, t = t_pt, h = h)
  # x range: [0, 0.4], width = 0.4 < 2h = 0.6 -> vol < (2h)^2
  # n_total = round(N * (2h)^2 / vol) > N
  expect_gt(s$n_total, N)
})

test_that("sampler_owin: points lie within the owin", {
  skip_if_not_installed("spatstat.geom")
  set.seed(23L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  t_pt <- c(0.5, 0.5)
  h <- 0.4
  s <- sampler_owin(win)(150L, t = t_pt, h = h)
  pts_x <- s$points[1L, ] + t_pt[1L]
  pts_y <- s$points[2L, ] + t_pt[2L]
  expect_true(all(pts_x >= 0 - 1e-9 & pts_x <= 1 + 1e-9))
  expect_true(all(pts_y >= 0 - 1e-9 & pts_y <= 1 + 1e-9))
})
