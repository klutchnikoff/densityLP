test_that("domain_Rd: correct class and dimension", {
  dom <- domain_Rd(2L)
  expect_s3_class(dom, "lp_domain")
  expect_equal(dom$d, 2L)
  expect_true(is.function(dom$sampler_factory))
})

test_that("domain_Rd: sampler returns points of correct dimension", {
  set.seed(10L)
  dom <- domain_Rd(3L)
  s <- dom$sampler_factory()(50L, t = c(0.5, 0.5, 0.5), h = 0.2)
  expect_equal(nrow(s$points), 3L)
  expect_true(s$n_total > 0L)
})

test_that("domain_Rd: all sampled points lie in [-h, h]^d", {
  set.seed(11L)
  dom <- domain_Rd(2L)
  h <- 0.3
  s <- dom$sampler_factory()(200L, t = c(0.4, 0.6), h = h)
  expect_true(all(abs(s$points) <= h))
})

test_that("domain_Rd: all N_quad points accepted (no rejection)", {
  set.seed(12L)
  d <- 2L
  h <- 0.4
  N_quad <- 2000L
  dom <- domain_Rd(d)
  s <- dom$sampler_factory()(N_quad, t = rep(0.5, d), h = h)
  # domain_Rd accepts everything: N_in = N_quad and n_total = N_quad
  expect_equal(ncol(s$points), N_quad)
  expect_equal(s$n_total, N_quad)
})

test_that("domain_from_indicator: rejects points outside unit disk", {
  set.seed(20L)
  is_in <- function(X) rowSums(X^2) <= 1
  dom <- domain_from_indicator(is_in, d = 2L)
  expect_s3_class(dom, "lp_domain")
  s <- dom$sampler_factory()(100L, t = c(0, 0), h = 0.5)
  pts_global <- sweep(t(s$points), 2L, c(0, 0), "+")
  expect_true(all(rowSums(pts_global^2) <= 1 + 1e-9))
})

test_that("domain_from_indicator: n_total is positive and points are accepted", {
  set.seed(21L)
  is_in <- function(X) rowSums(X^2) <= 1
  dom <- domain_from_indicator(is_in, d = 2L)
  s <- dom$sampler_factory()(500L, t = c(0, 0), h = 0.5)
  expect_true(s$n_total > 0L)
  expect_true(ncol(s$points) > 0L)
})

test_that("domain_sector: correct class and dimension", {
  dom <- domain_sector(2L)
  expect_s3_class(dom, "lp_domain")
  expect_equal(dom$d, 2L)
})

test_that("domain_sector: sampled points lie in D_k", {
  set.seed(30L)
  k <- 2
  dom <- domain_sector(k)
  t_pt <- c(0.5, 0.1)
  h <- 0.2
  s <- dom$sampler_factory()(200L, t = t_pt, h = h)
  pts <- sweep(t(s$points), 2L, t_pt, "+")
  x <- pts[, 1L]
  y <- pts[, 2L]
  expect_true(all(x >= 0 - 1e-9))
  expect_true(all(x <= 1 + 1e-9))
  expect_true(all(y >= 0 - 1e-9))
  expect_true(all(y <= x^k + 1e-9))
})

test_that("domain_sector: n_total is positive and finite", {
  set.seed(31L)
  dom <- domain_sector(1)
  s <- dom$sampler_factory()(300L, t = c(0.5, 0.2), h = 0.2)
  expect_true(is.finite(s$n_total))
  expect_true(s$n_total > 0L)
})

test_that("check_domain: error on wrong class", {
  expect_error(check_domain(list(d = 2L), d = 2L), "lp_domain")
})

test_that("check_domain: error on dimension mismatch", {
  dom <- domain_Rd(2L)
  expect_error(check_domain(dom, d = 3L), "dimension")
})

# ── check_X / check_t_grid ─────────────────────────────────────────────────────

test_that("check_X: error if not a matrix", {
  expect_error(check_X(as.data.frame(matrix(1, 3, 2))), "matrix")
})

test_that("check_X: error if not numeric", {
  expect_error(check_X(matrix("a", 2, 2)), "numeric")
})

test_that("check_X: error if contains NA", {
  X <- matrix(1:4, 2, 2)
  X[1, 1] <- NA_real_
  expect_error(check_X(X), "missing")
})

test_that("check_X: error if zero rows", {
  expect_error(check_X(matrix(numeric(0), 0, 2)), "one observation")
})

test_that("check_t_grid: error if not a matrix", {
  expect_error(check_t_grid(c(0.5, 0.5), d = 2L), "matrix")
})

test_that("check_t_grid: error on dimension mismatch", {
  expect_error(check_t_grid(matrix(0.5, 1, 3), d = 2L), "column")
})

test_that("check_t_grid: error if contains NA", {
  t <- matrix(c(0.5, NA_real_), 1, 2)
  expect_error(check_t_grid(t, d = 2L), "missing")
})

# ── print.lp_domain ────────────────────────────────────────────────────────────

test_that("print.lp_domain: returns x invisibly", {
  dom <- domain_Rd(2L)
  expect_identical(withVisible(print(dom))$visible, FALSE)
})

test_that("print.lp_domain: produces output", {
  dom <- domain_Rd(2L)
  expect_output(print(dom), "lp_domain")
})

# ── domain_from_owin ───────────────────────────────────────────────────────────

test_that("domain_from_owin: returns lp_domain with d = 2", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  dom <- domain_from_owin(win)
  expect_s3_class(dom, "lp_domain")
  expect_equal(dom$d, 2L)
})

test_that("domain_from_owin: label contains 'owin'", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::disc()
  dom <- domain_from_owin(win)
  expect_match(dom$label, "owin")
})

test_that("domain_from_owin: error on non-owin input", {
  expect_error(domain_from_owin(list(x = 1:3)), "owin")
})

test_that("domain_from_owin: sampled points lie inside the window", {
  skip_if_not_installed("spatstat.geom")
  set.seed(40L)
  win <- spatstat.geom::owin(c(0, 1), c(0, 1))
  dom <- domain_from_owin(win)
  t_pt <- c(0.5, 0.5)
  h <- 0.3
  s <- dom$sampler_factory()(100L, t = t_pt, h = h)
  pts_x <- s$points[1L, ] + t_pt[1L]
  pts_y <- s$points[2L, ] + t_pt[2L]
  expect_true(all(pts_x >= 0 - 1e-9 & pts_x <= 1 + 1e-9))
  expect_true(all(pts_y >= 0 - 1e-9 & pts_y <= 1 + 1e-9))
})
