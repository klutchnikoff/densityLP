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
#' @param box Numeric vector of length `2d` encoding the bounding box as
#'   `c(x1_min, x1_max, x2_min, x2_max, ...)`.
#' @param M Scalar upper bound of `f` on the domain.  Estimated automatically
#'   from a pilot sample when `NULL`.
#' @param n_pilot Number of pilot draws used to estimate `M`.
#' @return A function `function(n)` returning an `n x d` numeric matrix of
#'   random points, with column names taken from the formals of `f`.
#' @examples
#' f            <- function(x, y) x^2 + abs(y)^3
#' is_in_domain <- function(x, y) x^2 + y^2 <= 1
#' box          <- c(-1, 1, -1, 1)
#' rpts         <- rfactory(f, is_in_domain, box)
#' pts          <- rpts(200L)
#'
#' # Simplified: masking embedded in f
#' g    <- function(x, y) (x^2 + abs(y)^3) * (x^2 + y^2 <= 1)
#' rpts <- rfactory(g, box = c(-1, 1, -1, 1))
#' @importFrom stats runif
#' @export
rfactory <- function(f, is_in_domain = NULL, box, M = NULL, n_pilot = 2000L) {
  d <- length(formals(f))
  stopifnot(d >= 1L)
  stopifnot(is.numeric(box), length(box) == 2L * d)

  lower <- box[seq(1L, 2L * d - 1L, by = 2L)]
  upper <- box[seq(2L, 2L * d, by = 2L)]
  stopifnot(all(upper > lower))

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
    pilots <- .sample_box(n_pilot)
    fvals <- .call_fn(f, pilots)
    if (!is.null(is_in_domain)) {
      fvals <- fvals * .call_fn(is_in_domain, pilots)
    }
    if (max(fvals) <= 0) {
      stop("All pilot values are zero or negative. Check `box` and `f`.")
    }
    M <- max(fvals) * 1.1
  }
  stopifnot(is.numeric(M), length(M) == 1L, M > 0)

  function(n) {
    n <- as.integer(n)
    batch <- max(1000L, n * 10L)
    acc <- vector("list")
    total <- 0L

    while (total < n) {
      cands <- .sample_box(batch)

      if (!is.null(is_in_domain)) {
        cands <- cands[.call_fn(is_in_domain, cands), , drop = FALSE]
      }
      if (nrow(cands) == 0L) {
        next
      }

      fvals <- .call_fn(f, cands)
      keep <- runif(nrow(cands)) < fvals / M
      pts <- cands[keep, , drop = FALSE]

      if (nrow(pts) > 0L) {
        acc <- c(acc, list(pts))
        total <- total + nrow(pts)
      }
    }

    out <- do.call(rbind, acc)[seq_len(n), , drop = FALSE]
    colnames(out) <- coord_names
    out
  }
}
