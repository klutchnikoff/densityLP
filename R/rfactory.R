# Random point generators -------------------------------------------------------

#' Factory for random point generators from an unnormalized density
#'
#' Returns a closure `function(n)` that draws `n` points from the probability
#' distribution proportional to `f * 1_{is_in_domain}` via rejection sampling.
#'
#' The dimension `d` is inferred from the number of formal arguments of `f`.
#' Both `f` and `is_in_domain` must be vectorised: each argument receives a
#' numeric vector of length `n_candidates` (one per coordinate axis).
#'
#' If `is_in_domain` is omitted, `f` itself must return `0` outside the domain.
#'
#' @param f Vectorised positive function of `d` arguments (one per coordinate).
#'   Returns a numeric vector of length equal to the argument length.
#' @param is_in_domain Optional vectorised indicator: same signature as `f`,
#'   returns a logical vector.  When `NULL`, masking is embedded in `f`.
#' @param box Named list with elements `lower` and `upper` (numeric vectors of
#'   length `d`) bounding the support of `f`.
#' @param M Scalar upper bound of `f` on the domain.  Estimated automatically
#'   from a pilot sample when `NULL`.
#' @param N_pilot Number of pilot draws used to estimate `M`.
#' @return A function `function(n)` returning an `n x d` numeric matrix of
#'   random points, with column names taken from the formals of `f`.
#' @examples
#' f            <- function(x, y) x^2 + abs(y)^3
#' is_in_domain <- function(x, y) x^2 + y^2 <= 1
#' box          <- list(lower = c(-1, -1), upper = c(1, 1))
#' rpts         <- rfactory(f, is_in_domain, box)
#' pts          <- rpts(200L)
#'
#' # Simplified: masking embedded in f
#' g    <- function(x, y) (x^2 + abs(y)^3) * (x^2 + y^2 <= 1)
#' rpts <- rfactory(g, box = list(lower = c(-1, -1), upper = c(1, 1)))
#' @importFrom stats runif
#' @export
rfactory <- function(f, is_in_domain = NULL, box, M = NULL, N_pilot = 2000L) {
  d <- length(formals(f))
  stopifnot(d >= 1L)

  lower <- box$lower
  upper <- box$upper
  stopifnot(
    is.numeric(lower),
    is.numeric(upper),
    length(lower) == d,
    length(upper) == d,
    all(is.finite(lower)),
    all(is.finite(upper)),
    all(upper > lower)
  )

  coord_names <- names(formals(f))

  .call_fn <- function(fn, mat) {
    args <- lapply(seq_len(d), function(j) mat[, j])
    do.call(fn, args)
  }

  .sample_box <- function(n) {
    mat <- matrix(runif(n * d), n, d)
    sweep(sweep(mat, 2L, upper - lower, "*"), 2L, lower, "+")
  }

  if (is.null(M)) {
    pilots <- .sample_box(N_pilot)
    fvals <- .call_fn(f, pilots)
    if (!is.null(is_in_domain)) {
      fvals <- fvals * .call_fn(is_in_domain, pilots)
    }
    if (max(fvals) <= 0) {
      stop("All pilot values are zero or negative. Check `box` and `f`.")
    }
    M <- max(fvals) * 1.1
  }
  stopifnot(is.numeric(M), length(M) == 1L, is.finite(M), M > 0)

  function(n) {
    n <- as.integer(n)
    stopifnot(n >= 1L)
    batch <- max(1000L, n * 10L)
    out <- matrix(NA_real_, n, d)
    filled <- 0L

    while (filled < n) {
      cands <- .sample_box(batch)
      if (!is.null(is_in_domain)) {
        cands <- cands[.call_fn(is_in_domain, cands), , drop = FALSE]
      }
      if (nrow(cands) == 0L) {
        next
      }

      fvals <- .call_fn(f, cands)
      if (any(fvals > M)) {
        warning(
          "some fvals exceed M: M is underestimated, sample may be biased."
        )
      }
      keep <- runif(nrow(cands)) < fvals / M
      pts <- cands[keep, , drop = FALSE]

      take <- min(nrow(pts), n - filled)
      if (take > 0L) {
        out[seq(filled + 1L, filled + take), ] <- pts[seq_len(take), ]
        filled <- filled + take
      }
    }

    colnames(out) <- coord_names
    out
  }
}
