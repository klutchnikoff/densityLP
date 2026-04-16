test_that("build_alphas : cardinalite D_m = choose(m+d, d)", {
  expect_equal(nrow(build_alphas(0L, 1L)), choose(0 + 1, 1))
  expect_equal(nrow(build_alphas(1L, 1L)), choose(1 + 1, 1))
  expect_equal(nrow(build_alphas(2L, 2L)), choose(2 + 2, 2))
  expect_equal(nrow(build_alphas(3L, 3L)), choose(3 + 3, 3))
})

test_that("build_alphas : dimension des colonnes = d", {
  expect_equal(ncol(build_alphas(2L, 3L)), 3L)
  expect_equal(ncol(build_alphas(0L, 5L)), 5L)
})

test_that("build_alphas : tous les |alpha| <= m", {
  for (m in 0:3) {
    for (d in 1:3) {
      A <- build_alphas(m, d)
      expect_true(all(rowSums(A) <= m))
    }
  }
})

test_that("build_alphas : ordre correct (degre croissant, puis lex)", {
  A <- build_alphas(2L, 2L)
  deg <- rowSums(A)
  # Degrés non décroissants
  expect_true(all(diff(deg) >= 0))
  # Dans chaque bloc de même degré, ordre lexicographique colonne 1
  for (d_val in unique(deg)) {
    bloc <- A[deg == d_val, , drop = FALSE]
    expect_equal(
      bloc,
      bloc[
        do.call(order, lapply(seq_len(ncol(bloc)), function(j) bloc[, j])),
        ,
        drop = FALSE
      ]
    )
  }
})

test_that("build_alphas : valeurs exactes pour m=1, d=2", {
  A <- build_alphas(1L, 2L)
  # Ordre : degre 0 -> (0,0) ; degre 1, lex -> (0,1) < (1,0)
  expected <- matrix(c(0L, 0L, 0L, 1L, 1L, 0L), nrow = 3L, byrow = TRUE)
  colnames(expected) <- c("x1", "x2")
  expect_equal(A, expected)
})

test_that("build_alphas : valeurs exactes pour m=2, d=1", {
  A <- build_alphas(2L, 1L)
  expected <- matrix(0L:2L, nrow = 3L, ncol = 1L)
  colnames(expected) <- "x1"
  expect_equal(A, expected)
})

test_that("build_alphas : rejette les entrees invalides", {
  expect_error(build_alphas(-1L, 2L))
  expect_error(build_alphas(1L, 0L))
  expect_error(build_alphas(1.5, 2L))
})

# ── build_Phi ------------------------------------------------------------------

test_that("build_Phi : dimensions correctes", {
  alphas <- build_alphas(2L, 2L)
  U <- matrix(runif(2L * 10L), 2L, 10L)
  Phi <- build_Phi(U, h = 1, alphas = alphas)
  expect_equal(dim(Phi), c(nrow(alphas), 10L))
})

test_that("build_Phi : monome constant = 1 en u = 0 (alpha = (0,...,0))", {
  for (d in 1:3) {
    alphas <- build_alphas(2L, d)
    U0 <- matrix(0, d, 1L)
    Phi0 <- build_Phi(U0, h = 1, alphas = alphas)
    # Premier monome (|alpha|=0) vaut toujours 1
    expect_equal(Phi0[1L, 1L], 1)
    # Les monomes de degre >= 1 valent 0 en u=0
    expect_true(all(Phi0[-1L, 1L] == 0))
  }
})

test_that("build_Phi : valeurs exactes pour m=1, d=2, h=1", {
  alphas <- build_alphas(1L, 2L) # (0,0), (0,1), (1,0)
  u <- c(0.3, 0.7)
  U <- matrix(u, 2L, 1L)
  Phi <- build_Phi(U, h = 1, alphas = alphas)
  # alpha=(0,0)->1, alpha=(0,1)->u[2]=0.7, alpha=(1,0)->u[1]=0.3
  expect_equal(as.numeric(Phi), c(1, u[2], u[1]))
})

test_that("build_Phi : rescaling par h correct", {
  alphas <- build_alphas(1L, 1L) # (0), (1)
  u <- 0.4
  h <- 2
  U <- matrix(u, 1L, 1L)
  Phi <- build_Phi(U, h = h, alphas = alphas)
  expect_equal(as.numeric(Phi), c(1, u / h))
})
