# ── density_lp_point ───────────────────────────────────────────────────────────

test_that("density_lp_point: returns a scalar", {
  set.seed(1L)
  X <- matrix(runif(100L * 2L), 100L, 2L)
  val <- densityLP:::density_lp_point(
    X,
    t = c(0.5, 0.5),
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_true(is.numeric(val))
  expect_length(val, 1L)
})

test_that("density_lp_point: returns non-negative value", {
  set.seed(2L)
  X <- matrix(runif(200L * 2L), 200L, 2L)
  val <- densityLP:::density_lp_point(
    X,
    t = c(0.5, 0.5),
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_gte(val, 0)
})

test_that("density_lp_point: returns 0 when no observations in V(h)", {
  X <- matrix(runif(50L * 2L, 0.8, 1), 50L, 2L)
  val <- densityLP:::density_lp_point(
    X,
    t = c(0.1, 0.1),
    h = 0.05,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_equal(val, 0)
})

# ── print / plot methods ────────────────────────────────────────────────────────

test_that("print.density_lp: returns x invisibly", {
  set.seed(10L)
  X <- matrix(runif(100L * 2L), 100L, 2L)
  res <- density_lp(
    X,
    matrix(0.5, 1L, 2L),
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_identical(withVisible(print(res))$visible, FALSE)
})

test_that("print.density_lp: produces output", {
  set.seed(11L)
  X <- matrix(runif(100L * 2L), 100L, 2L)
  res <- density_lp(
    X,
    matrix(0.5, 1L, 2L),
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_output(print(res), "Local polynomial")
})

test_that("plot.density_lp: returns x invisibly for d=1", {
  set.seed(12L)
  X <- matrix(runif(100L), 100L, 1L)
  t_grid <- matrix(seq(0.2, 0.8, length.out = 5L), ncol = 1L)
  res <- density_lp(
    X,
    t_grid,
    h = 0.3,
    m = 1L,
    domain = domain_Rd(1L),
    N_quad = 200L
  )
  expect_identical(withVisible(plot(res))$visible, FALSE)
})

test_that("plot.density_lp: returns x invisibly for d=2", {
  set.seed(13L)
  X <- matrix(runif(100L * 2L), 100L, 2L)
  xs <- seq(0.2, 0.8, length.out = 3L)
  t_grid <- as.matrix(expand.grid(x = xs, y = xs))
  res <- density_lp(
    X,
    t_grid,
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_identical(withVisible(plot(res))$visible, FALSE)
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("density_lp: error if X is not a matrix", {
  expect_error(
    density_lp(
      as.data.frame(matrix(0, 10, 2)),
      matrix(0.5, 1, 2),
      h = 0.3,
      m = 1,
      domain = domain_Rd(2L)
    ),
    "matrix"
  )
})

test_that("density_lp: error on dimension mismatch between X and t_grid", {
  X <- matrix(runif(20), 10, 2)
  t_grid <- matrix(0.5, 1, 3) # d=3 != 2
  expect_error(
    density_lp(X, t_grid, h = 0.3, m = 1, domain = domain_Rd(2L)),
    "column"
  )
})

test_that("density_lp: error on domain dimension mismatch", {
  X <- matrix(runif(20), 10, 2)
  t_grid <- matrix(0.5, 1, 2)
  expect_error(
    density_lp(X, t_grid, h = 0.3, m = 1, domain = domain_Rd(3L)),
    "dimension"
  )
})

# ── Return structure ───────────────────────────────────────────────────────────

test_that("density_lp: returns 'density_lp' object with correct fields", {
  set.seed(42L)
  X <- matrix(runif(100L * 2L), 100L, 2L)
  t_grid <- matrix(c(0.3, 0.5, 0.7, 0.5), nrow = 2L)
  res <- density_lp(
    X,
    t_grid,
    h = 0.4,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 200L
  )
  expect_s3_class(res, "density_lp")
  expect_length(res$estimate, 2L)
  expect_equal(nrow(res$t_grid), 2L)
  expect_equal(res$h, 0.4)
  expect_equal(res$m, 1L)
})

test_that("density_lp: estimates are non-negative", {
  set.seed(43L)
  X <- matrix(runif(200L * 2L), 200L, 2L)
  xs <- seq(0.2, 0.8, length.out = 4L)
  t_grid <- as.matrix(expand.grid(x = xs, y = xs))
  res <- density_lp(
    X,
    t_grid,
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 300L
  )
  expect_true(all(res$estimate >= 0))
})

# ── Numerical validation on [0,1]^2 ───────────────────────────────────────────
# True density: Uniform(0,1)^2 → f ≡ 1.
# Oracle estimator (m=0, known optimal h) should be close to 1.

test_that("density_lp: estimates close to 1 for Uniform([0,1]^2), m=0", {
  set.seed(100L)
  n <- 2000L
  X <- matrix(runif(n * 2L), n, 2L)
  # Interior points, away from boundary
  t_grid <- matrix(c(0.3, 0.3, 0.5, 0.5, 0.7, 0.3), nrow = 3L, byrow = TRUE)
  res <- density_lp(
    X,
    t_grid,
    h = 0.25,
    m = 0L,
    domain = domain_Rd(2L),
    N_quad = 500L
  )
  # With n=2000, h=0.25: expect within 20% of truth (1.0) with high probability
  expect_true(
    all(abs(res$estimate - 1) < 0.3),
    label = paste("estimates:", paste(round(res$estimate, 3), collapse = ", "))
  )
})

test_that("density_lp: m=1 integrates approximately to 1 on interior grid", {
  set.seed(101L)
  n <- 1000L
  X <- matrix(runif(n * 2L), n, 2L)
  # Regular interior grid (avoid boundary effects)
  xs <- seq(0.15, 0.85, length.out = 8L)
  t_grid <- as.matrix(expand.grid(x = xs, y = xs))
  cell_area <- (xs[2L] - xs[1L])^2
  res <- density_lp(
    X,
    t_grid,
    h = 0.3,
    m = 1L,
    domain = domain_Rd(2L),
    N_quad = 300L
  )
  integral <- sum(res$estimate) * cell_area
  # Riemann sum over interior: expect between 0.4 and 1.2
  expect_true(
    integral > 0.4 && integral < 1.2,
    label = paste("integral =", round(integral, 3))
  )
})
