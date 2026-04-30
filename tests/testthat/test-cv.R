# ── cv_density_lp ─────────────────────────────────────────────────────────────

set_seed <- function() set.seed(42L)

# Shared fixture: uniform on [0,1]^2, small n for speed
make_square_cv <- function(
  n = 80L,
  h_grid = c(0.2, 0.3, 0.4),
  m_grid = 0:1,
  N_quad = 150L
) {
  set_seed()
  X <- matrix(runif(n * 2L), n, 2L)
  is_sq <- function(x, y) x >= 0 & x <= 1 & y >= 0 & y <= 1
  domain <- domain_from_indicator(is_sq)
  cv_density_lp(
    X,
    h_grid = h_grid,
    m_grid = m_grid,
    domain = domain,
    N_quad = N_quad
  )
}

test_that("cv_density_lp: returns cv_density_lp class", {
  cv <- make_square_cv()
  expect_s3_class(cv, "cv_density_lp")
})

test_that("cv_density_lp: score matrix has correct dimensions", {
  h_grid <- c(0.2, 0.3, 0.4)
  m_grid <- 0:2
  set_seed()
  X <- matrix(runif(60L * 2L), 60L, 2L)
  is_sq <- function(x, y) x >= 0 & x <= 1 & y >= 0 & y <= 1
  cv <- cv_density_lp(
    X,
    h_grid = h_grid,
    m_grid = m_grid,
    domain = domain_from_indicator(is_sq),
    N_quad = 100L
  )
  expect_equal(dim(cv$scores), c(length(m_grid), length(h_grid)))
})

test_that("cv_density_lp: m_hat is in m_grid", {
  cv <- make_square_cv()
  expect_true(cv$m_hat %in% cv$grids$m)
})

test_that("cv_density_lp: h_hat is in h_grid", {
  cv <- make_square_cv()
  expect_true(cv$h_hat %in% cv$grids$h)
})

test_that("cv_density_lp: selected score is the minimum finite score", {
  cv <- make_square_cv()
  row_idx <- which(cv$grids$m == cv$m_hat)
  col_idx <- which(cv$grids$h == cv$h_hat)
  best <- cv$scores[row_idx, col_idx]
  expect_equal(best, min(cv$scores, na.rm = TRUE))
})

test_that("cv_density_lp: grids stored in output match inputs", {
  h_grid <- c(0.25, 0.35)
  m_grid <- c(0L, 1L)
  set_seed()
  X <- matrix(runif(60L * 2L), 60L, 2L)
  is_sq <- function(x, y) x >= 0 & x <= 1 & y >= 0 & y <= 1
  cv <- cv_density_lp(
    X,
    h_grid = h_grid,
    m_grid = m_grid,
    domain = domain_from_indicator(is_sq),
    N_quad = 100L
  )
  expect_equal(cv$grids$h, h_grid)
  expect_equal(cv$grids$m, m_grid)
})

test_that("cv_density_lp: allows data.frame and error if input is a list", {
  X_df <- as.data.frame(matrix(runif(20), 10, 2))
  dom <- domain_Rd(2L)
  expect_no_error(cv_density_lp(X_df, h_grid = 0.3, domain = dom))
  expect_error(cv_density_lp(as.list(X_df), h_grid = 0.3, domain = dom))
})

test_that("cv_density_lp: error if h_grid contains non-positive value", {
  set_seed()
  X <- matrix(runif(40L * 2L), 40L, 2L)
  expect_error(
    cv_density_lp(X, h_grid = c(0.2, -0.1), domain = domain_Rd(2L))
  )
})

# ── print.cv_density_lp ───────────────────────────────────────────────────────

test_that("print.cv_density_lp: returns x invisibly", {
  cv <- make_square_cv()
  expect_identical(withVisible(print(cv))$visible, FALSE)
})

test_that("print.cv_density_lp: output mentions selected m and h", {
  cv <- make_square_cv()
  expect_output(print(cv), "Selected m")
  expect_output(print(cv), "Selected h")
})

# ── cv_density_lp_ppp ─────────────────────────────────────────────────────────

test_that("cv_density_lp_ppp: returns cv_density_lp class", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set_seed()
  pp <- spatstat.random::runifpoint(80L, win = spatstat.geom::disc())
  cv <- cv_density_lp_ppp(pp, h_grid = c(0.3, 0.5), m_grid = 0:1, N_quad = 150L)
  expect_s3_class(cv, "cv_density_lp")
})

test_that("cv_density_lp_ppp: m_hat and h_hat are in the grids", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")
  set_seed()
  pp <- spatstat.random::runifpoint(80L, win = spatstat.geom::disc())
  h_grid <- c(0.3, 0.5)
  m_grid <- 0:1
  cv <- cv_density_lp_ppp(pp, h_grid = h_grid, m_grid = m_grid, N_quad = 150L)
  expect_true(cv$m_hat %in% m_grid)
  expect_true(cv$h_hat %in% h_grid)
})

test_that("cv_density_lp_ppp: error if input is not a ppp", {
  expect_error(cv_density_lp_ppp(matrix(1, 10, 2), h_grid = 0.3), "ppp")
})

# ── Integer storage ────────────────────────────────────────────────────────────

test_that("cv_density_lp: m_grid stored as integer in grids$m", {
  cv <- make_square_cv(m_grid = c(0, 1)) # passed as double
  expect_type(cv$grids$m, "integer")
})

# ── print CV score not NA ──────────────────────────────────────────────────────

test_that("print.cv_density_lp: CV score is not NA", {
  cv <- make_square_cv()
  row_idx <- which(cv$grids$m == cv$m_hat)
  col_idx <- which(cv$grids$h == cv$h_hat)
  score <- cv$scores[row_idx, col_idx]
  expect_false(is.na(score))
  expect_true(is.finite(score))
})

# ── summary.cv_density_lp ─────────────────────────────────────────────────────

test_that("summary.cv_density_lp: produces output and returns invisibly", {
  cv <- make_square_cv()
  expect_output(summary(cv), "Selected")
  expect_output(summary(cv), "CV score matrix")
  expect_identical(withVisible(summary(cv))$visible, FALSE)
})

# ── Warning on high failure rate ───────────────────────────────────────────────

test_that("cv_density_lp: warns when >25% LOO evaluations fail", {
  set_seed()
  n <- 30L
  X <- matrix(runif(n * 2L), n, 2L)
  dom <- domain_Rd(2L)
  # h very small -> most V(h) neighbourhoods empty -> high failure rate
  expect_warning(
    cv_density_lp(X, h_grid = c(0.01), m_grid = 0L, domain = dom, N_quad = 50L),
    "failed"
  )
})
