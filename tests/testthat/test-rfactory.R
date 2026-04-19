# ── rfactory ──────────────────────────────────────────────────────────────────

f_disc <- function(x, y) x^2 + abs(y)^3
in_disc <- function(x, y) x^2 + y^2 <= 1
box_disc <- list(lower = c(-1, -1), upper = c(1, 1))

test_that("rfactory: returns a function", {
  r <- rfactory(f_disc, in_disc, box_disc)
  expect_true(is.function(r))
})

test_that("rfactory: output is a matrix with correct dimensions", {
  set.seed(1L)
  r <- rfactory(f_disc, in_disc, box_disc)
  pts <- r(100L)
  expect_true(is.matrix(pts))
  expect_equal(dim(pts), c(100L, 2L))
})

test_that("rfactory: column names match formals of f", {
  set.seed(2L)
  r <- rfactory(f_disc, in_disc, box_disc)
  pts <- r(50L)
  expect_equal(colnames(pts), c("x", "y"))
})

test_that("rfactory: all points satisfy the domain indicator", {
  set.seed(3L)
  r <- rfactory(f_disc, in_disc, box_disc)
  pts <- r(200L)
  expect_true(all(in_disc(pts[, "x"], pts[, "y"])))
})

test_that("rfactory: non-uniform density shifts the distribution", {
  set.seed(4L)
  r <- rfactory(function(x, y) x^2, in_disc, box_disc)
  pts <- r(500L)
  expect_gt(mean(abs(pts[, "x"])), 0.6)
})

test_that("rfactory: simplified form (masking in f) works", {
  set.seed(5L)
  g <- function(x, y) (x^2 + abs(y)^3) * (x^2 + y^2 <= 1)
  r <- rfactory(g, box = box_disc)
  pts <- r(200L)
  expect_equal(dim(pts), c(200L, 2L))
  expect_true(all(pts[, "x"]^2 + pts[, "y"]^2 <= 1 + 1e-12))
})

test_that("rfactory: works in d = 1", {
  set.seed(6L)
  r <- rfactory(
    function(x) x^2,
    function(x) x >= 0 & x <= 1,
    box = list(lower = 0, upper = 1)
  )
  pts <- r(100L)
  expect_equal(dim(pts), c(100L, 1L))
  expect_equal(colnames(pts), "x")
  expect_true(all(pts >= 0 & pts <= 1))
})

test_that("rfactory: user-supplied M is respected", {
  set.seed(7L)
  r <- rfactory(f_disc, in_disc, box_disc, M = 2.0)
  pts <- r(100L)
  expect_equal(dim(pts), c(100L, 2L))
})

test_that("rfactory: error if lower >= upper", {
  expect_error(
    rfactory(f_disc, in_disc, list(lower = c(1, -1), upper = c(-1, 1)))
  )
})

test_that("rfactory: error if all pilot values are zero", {
  expect_error(
    rfactory(function(x, y) 0, in_disc, box_disc, N_pilot = 200L),
    "zero"
  )
})
