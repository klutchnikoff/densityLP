test_that("gram_matrix: output is square D_m x D_m", {
  set.seed(1L)
  for (m in 0:2) {
    for (d in 1:3) {
      alphas <- build_alphas(m, d)
      Dm <- nrow(alphas)
      N <- 200L
      pts <- matrix(runif(d * N, -0.5, 0.5), d, N)
      s <- list(points = pts, n_total = N) # vol = (2*0.5)^d, n_total = N
      B <- gram_matrix(s, h = 0.5, alphas = alphas)
      expect_equal(dim(B), c(Dm, Dm))
    }
  }
})

test_that("gram_matrix: output is symmetric", {
  set.seed(2L)
  alphas <- build_alphas(2L, 2L)
  N <- 300L
  pts <- matrix(runif(2L * N, -0.3, 0.3), 2L, N)
  s <- list(points = pts, n_total = N) # vol = (0.6)^2 = 0.36, n_total = N
  B <- gram_matrix(s, h = 0.3, alphas = alphas)
  expect_equal(B, t(B))
})

test_that("gram_matrix: output is positive definite (all eigenvalues > 0)", {
  set.seed(3L)
  alphas <- build_alphas(2L, 2L)
  N <- 500L
  pts <- matrix(runif(2L * N, -0.4, 0.4), 2L, N)
  s <- list(points = pts, n_total = N) # vol = (0.8)^2 = 0.64, n_total = N
  B <- gram_matrix(s, h = 0.4, alphas = alphas)
  ev <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("gram_matrix: m=0 gives 1x1 matrix equal to (2h)^d / h^d = 2^d", {
  # B = (2h)^d / (n_total * h^d) * sum_k phi_0(u_k)^2
  #   = (2h)^d / (n_total * h^d) * N
  # With n_total = N (full box), B[1,1] = (2h)^d / h^d = 2^d.
  set.seed(4L)
  d <- 2L
  h <- 0.5
  N <- 100L
  alphas <- build_alphas(0L, d)
  pts <- matrix(runif(d * N, -h, h), d, N)
  s <- list(points = pts, n_total = N)
  B <- gram_matrix(s, h = h, alphas = alphas)
  expect_equal(dim(B), c(1L, 1L))
  expect_equal(B[1L, 1L], (2 * h)^d / h^d, tolerance = 1e-12)
})

test_that("gram_matrix: Cholesky succeeds (implies positive definiteness)", {
  set.seed(5L)
  alphas <- build_alphas(3L, 2L)
  N <- 1000L
  pts <- matrix(runif(2L * N, -0.5, 0.5), 2L, N)
  s <- list(points = pts, n_total = N)
  B <- gram_matrix(s, h = 0.5, alphas = alphas)
  expect_no_error(chol(B))
})
