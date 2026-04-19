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
check_h <- function(h) {
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("`h` must be a finite positive scalar.")
  }
}

#' @keywords internal
check_m <- function(m) {
  if (
    !is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 0 || m != floor(m)
  ) {
    stop("`m` must be a finite non-negative integer.")
  }
}

#' @keywords internal
check_N_quad <- function(N_quad) {
  if (
    !is.numeric(N_quad) ||
      length(N_quad) != 1L ||
      !is.finite(N_quad) ||
      N_quad < 1L
  ) {
    stop("`N_quad` must be a finite positive integer.")
  }
}

#' @keywords internal
check_h_grid <- function(h_grid) {
  if (
    !is.numeric(h_grid) ||
      length(h_grid) < 1L ||
      !all(is.finite(h_grid)) ||
      !all(h_grid > 0)
  ) {
    stop("`h_grid` must be a non-empty vector of finite positive values.")
  }
}

#' @keywords internal
check_m_grid <- function(m_grid) {
  if (
    !is.numeric(m_grid) ||
      length(m_grid) < 1L ||
      !all(is.finite(m_grid)) ||
      !all(m_grid >= 0) ||
      !all(m_grid == floor(m_grid))
  ) {
    stop("`m_grid` must be a non-empty vector of finite non-negative integers.")
  }
}

#' @keywords internal
check_domain <- function(domain, d) {
  if (!inherits(domain, "domain_lp")) {
    stop(
      "`domain` must be an 'domain_lp' object ",
      "(see domain_Rd(), domain_from_indicator(), domain_sector())."
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
