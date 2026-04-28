#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Gram matrix weight: vol([-h,h]^d) / (n_total * h^d) = 2^d / n_total.
static inline double gram_weight(int d, int n_total) {
  return std::pow(2.0, d) / static_cast<double>(n_total);
}

// Integer power: x^n for small non-negative n.
// Faster than std::pow for the integer exponents arising from multi-indices.
static inline double ipow(double x, int n) {
  if (n == 0) return 1.0;
  if (n == 1) return x;
  double r = 1.0;
  while (n-- > 0) r *= x;
  return r;
}

// Build the Phi matrix (Dm x N): Phi(j,k) = prod_l (U_sc(l,k))^alphas(j,l)
// U_sc must already be rescaled (U / h), column-major d x N.
static arma::mat build_phi_mat(const arma::mat& U_sc,
                                const arma::Mat<int>& alphas) {
  int d  = U_sc.n_rows;
  int N  = U_sc.n_cols;
  int Dm = alphas.n_rows;

  arma::mat Phi(Dm, N, arma::fill::ones);
  for (int k = 0; k < N; k++) {
    for (int j = 0; j < Dm; j++) {
      double val = 1.0;
      for (int l = 0; l < d; l++) {
        int a = alphas(j, l);
        if (a != 0) val *= ipow(U_sc(l, k), a);
      }
      Phi(j, k) = val;
    }
  }
  return Phi;
}

//' Compute the Gram matrix B_gamma in C++
//'
//' @param U Numeric matrix d x N of quadrature points centred at t.
//' @param h Bandwidth (scalar > 0).
//' @param alphas Integer matrix D_m x d of multi-indices.
//' @param n_total Total draws in the box \eqn{[-h,h]^d} (from the sampler).
//' @return Symmetric matrix D_m x D_m.
//' @keywords internal
// [[Rcpp::export]]
arma::mat gram_matrix_cpp(const arma::mat&      U,
                          double                h,
                          const arma::Mat<int>& alphas,
                          int                   n_total) {
  int d = U.n_rows;

  arma::mat Phi = build_phi_mat(U / h, alphas);

  // B = weight * Phi * Phi^T  (BLAS DSYRK via Armadillo)
  return gram_weight(d, n_total) * (Phi * Phi.t());
}

//' Full LP estimator with LOO self-influence term and variance (C++ backend)
//'
//' @param U_quad Numeric matrix d x N of quadrature points centred at t.
//' @param n_total Integer, total draws in \eqn{[-h,h]^d}.
//' @param U_obs Numeric matrix d x n of all observations centred at t.
//' @param h Bandwidth (scalar > 0).
//' @param alphas Integer matrix D_m x d of multi-indices.
//' @return Numeric vector of length 3: \code{c(estimate, norm_H0_sq, variance)}.
//' @keywords internal
// [[Rcpp::export]]
arma::vec lp_estimator_loo_cpp(const arma::mat&      U_quad,
                               int                   n_total,
                               const arma::mat&      U_obs,
                               double                h,
                               const arma::Mat<int>& alphas) {
  int d  = U_quad.n_rows;
  int n  = U_obs.n_cols;
  int Dm = alphas.n_rows;

  // 1. Gram matrix
  arma::mat Phi_q = build_phi_mat(U_quad / h, alphas);
  arma::mat B = gram_weight(d, n_total) * (Phi_q * Phi_q.t());

  // 2. Cholesky: B = L_upper^T * L_upper
  arma::mat L_upper;
  if (!arma::chol(L_upper, B)) {
    Rcpp::stop("The Gram matrix is not positive definite. "
               "Try increasing N_quad or check the domain.");
  }
  arma::mat L_lower = L_upper.t();

  // 3. H_0 = L_lower^{-1} e_1
  arma::vec e1(Dm, arma::fill::zeros);
  e1(0) = 1.0;
  arma::vec H_0       = arma::solve(arma::trimatl(L_lower), e1);
  double norm_H0_sq   = arma::dot(H_0, H_0);  // ||H_0||^2 = (B^{-1})_{11}

  // 4. Filter observations in V(h) = [-h, h]^d
  arma::uvec idx(n);
  arma::uword N_V = 0;
  for (int i = 0; i < n; i++) {
    if (arma::all(arma::abs(U_obs.col(i)) <= h)) idx(N_V++) = i;
  }

  double estimate = 0.0;
  double variance = 0.0;
  if (N_V > 0) {
    // 5. Phi_V and H_V for observations in V(h)
    arma::mat Phi_V = build_phi_mat(U_obs.cols(idx.head(N_V)) / h, alphas);
    arma::mat H_V   = arma::solve(arma::trimatl(L_lower), Phi_V);

    // 6. Estimator and variance:
    // Z_i = h^-d * dot(H_0, H_V.col(i))
    // estimate = (1/n) * sum Z_i
    // variance = (1/n^2) * sum Z_i^2
    arma::vec Z = (H_V.t() * H_0) / std::pow(h, d);
    estimate = arma::sum(Z) / static_cast<double>(n);

    double sum_Z2 = arma::dot(Z, Z);
    variance = sum_Z2 / (static_cast<double>(n) * static_cast<double>(n));
  }

  return arma::vec{estimate, norm_H0_sq, variance};
}

//' Full local-polynomial density estimator at one point (C++ backend)
//'
//' @param U_quad Numeric matrix d x N of quadrature points centred at t
//'   (output of the domain sampler).
//' @param n_total Integer, total draws in \eqn{[-h,h]^d}.
//' @param U_obs Numeric matrix d x n of all observations centred at t
//'   (i.e. \eqn{X_i - t}, transposed so columns are observations).
//' @param h Bandwidth (scalar > 0).
//' @param alphas Integer matrix D_m x d of multi-indices.
//' @return Scalar density estimate \eqn{\hat f(t)}.
//' @keywords internal
// [[Rcpp::export]]
double lp_estimator_cpp(const arma::mat&      U_quad,
                        int                   n_total,
                        const arma::mat&      U_obs,
                        double                h,
                        const arma::Mat<int>& alphas) {
  return lp_estimator_loo_cpp(U_quad, n_total, U_obs, h, alphas)[0];
}

//' Optimized LOO CV scores for a fixed (m, h) grid point (C++ backend)
//'
//' @param X Numeric matrix d x n of all observations.
//' @param U_quad_list List of d x N_i matrices of quadrature points.
//' @param n_total_vec Integer vector of total draws in the box for each point.
//'   (length n).
//' @param h Bandwidth (scalar > 0).
//' @param alphas Integer matrix D_m x d of multi-indices.
//' @return Numeric vector of length n: log(f_loo(X_i)).
//' @keywords internal
// [[Rcpp::export]]
arma::vec cv_lp_fixed_h_m_cpp(const arma::mat&      X,
                              const Rcpp::List&    U_quad_list,
                              const arma::uvec&    n_total_vec,
                              double                h,
                              const arma::Mat<int>& alphas) {
  int n  = X.n_cols;
  int d  = X.n_rows;
  int Dm = alphas.n_rows;
  double h_pow_d = std::pow(h, d);
  arma::vec log_loo(n);

  for (int i = 0; i < n; i++) {
    arma::mat U_quad = Rcpp::as<arma::mat>(U_quad_list[i]);
    // Return -inf (not an error) so the loop continues for other observations.
    // The R side counts !is.finite(log_loo) as n_fail and warns accordingly.
    if (U_quad.n_cols < 2) {          // empty neighbourhood V(h)
      log_loo(i) = -arma::datum::inf;
      continue;
    }

    // 1. Gram matrix for observation i
    arma::mat Phi_q = build_phi_mat(U_quad / h, alphas);
    arma::mat B = gram_weight(d, n_total_vec(i)) * (Phi_q * Phi_q.t());

    // 2. Cholesky
    arma::mat L_upper;
    if (!arma::chol(L_upper, B)) {    // singular Gram matrix
      log_loo(i) = -arma::datum::inf;
      continue;
    }
    arma::mat L_lower = L_upper.t();

    // 3. H_0 = L_lower^{-1} e_1
    arma::vec e1(Dm, arma::fill::zeros);
    e1(0) = 1.0;
    arma::vec H_0 = arma::solve(arma::trimatl(L_lower), e1);
    double norm_H0_sq = arma::dot(H_0, H_0);

    // 4. Filter neighbors of X_i in the h-box directly from X
    // Instead of creating U_obs = X - X_i, we filter indices
    arma::uvec neighbors_idx(n);
    arma::uword N_V = 0;
    for (int j = 0; j < n; j++) {
      bool in_box = true;
      for (int l = 0; l < d; l++) {
        if (std::abs(X(l, j) - X(l, i)) > h) {
          in_box = false;
          break;
        }
      }
      if (in_box) neighbors_idx(N_V++) = j;
    }

    double estimate = 0.0;
    if (N_V > 0) {
      // 5. Phi_V and H_V for observations in V(h)
      // We must center them for the estimator: (X_j - X_i) / h
      arma::mat U_V(d, N_V);
      for (arma::uword k = 0; k < N_V; k++) {
        U_V.col(k) = (X.col(neighbors_idx(k)) - X.col(i)) / h;
      }
      arma::mat Phi_V = build_phi_mat(U_V, alphas);
      arma::mat H_V   = arma::solve(arma::trimatl(L_lower), Phi_V);

      estimate = arma::dot(H_0, arma::sum(H_V, 1)) / (static_cast<double>(n) * h_pow_d);
    }

    // 6. LOO correction
    double f_loo_i = (static_cast<double>(n) * estimate - norm_H0_sq / h_pow_d) / (static_cast<double>(n) - 1.0);
    log_loo(i) = (f_loo_i > 0) ? std::log(f_loo_i) : -arma::datum::inf;
  }

  return log_loo;
}
