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
