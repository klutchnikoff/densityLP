# Input checks (internal) -------------------------------------------------------

#' @keywords internal
check_X <- function(X) {
  if (!is.matrix(X)) {
    stop("`X` must be a numeric matrix (n x d). Use `as.matrix(X)` if needed.")
  }
  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  if (anyNA(X)) {
    stop("`X` contains missing values (NA).")
  }
  if (nrow(X) < 1L) {
    stop("`X` must have at least one observation.")
  }
}

#' @keywords internal
check_t_grid <- function(t_grid, d) {
  if (!is.matrix(t_grid)) {
    stop("`t_grid` must be a numeric matrix (p x d).")
  }
  if (!is.numeric(t_grid)) {
    stop("`t_grid` must be a numeric matrix.")
  }
  if (ncol(t_grid) != d) {
    stop(
      "`t_grid` must have ",
      d,
      " column(s) (d = ",
      d,
      "), but has ",
      ncol(t_grid),
      "."
    )
  }
  if (anyNA(t_grid)) {
    stop("`t_grid` contains missing values (NA).")
  }
}

#' @keywords internal
check_domain <- function(domain, d) {
  if (!inherits(domain, "lp_domain")) {
    stop(
      "`domain` must be an 'lp_domain' object ",
      "(see domain_Rd(), domain_func(), domain_sector())."
    )
  }
  if (domain$d != d) {
    stop(
      "Domain dimension (d = ",
      domain$d,
      ") does not match ",
      "the dimension of X (d = ",
      d,
      ")."
    )
  }
}
