#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
  double weight = std::pow(2.0 * h, d) /
                  (static_cast<double>(n_total) * std::pow(h, d));
  return weight * (Phi * Phi.t());
}

//' Full local-polynomial density estimator at one point (C++ backend)
//'
//' Performs steps 2–7 of the point estimator: Gram matrix, Cholesky,
//' forward solves, observation filtering, and final summation.
//' Step 1 (sampling) and alpha construction remain in R.
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
  int d  = U_quad.n_rows;
  int n  = U_obs.n_cols;
  int Dm = alphas.n_rows;

  // 1. Gram matrix
  arma::mat Phi_q = build_phi_mat(U_quad / h, alphas);
  double weight   = std::pow(2.0 * h, d) /
                    (static_cast<double>(n_total) * std::pow(h, d));
  arma::mat B = weight * (Phi_q * Phi_q.t());

  // 2. Cholesky: B = L_upper^T * L_upper
  arma::mat L_upper;
  if (!arma::chol(L_upper, B)) {
    Rcpp::stop("The Gram matrix is not positive definite. "
               "Try increasing N_quad or check the domain.");
  }
  arma::mat L_lower = L_upper.t();

  // 3. H_0 = L_lower^{-1} e_1   (e_1: first standard basis vector)
  arma::vec e1(Dm, arma::fill::zeros);
  e1(0) = 1.0;
  arma::vec H_0 = arma::solve(arma::trimatl(L_lower), e1);

  // 4. Filter observations in V(h) = [-h, h]^d
  arma::uvec idx(n);
  arma::uword N_V = 0;
  for (int i = 0; i < n; i++) {
    if (arma::all(arma::abs(U_obs.col(i)) <= h)) idx(N_V++) = i;
  }
  if (N_V == 0) return 0.0;

  // 5. Phi matrix for observations in V(h)  (Dm x N_V)
  arma::mat Phi_V = build_phi_mat(U_obs.cols(idx.head(N_V)) / h, alphas);

  // 6. H_V = L_lower^{-1} Phi_V   (Dm x N_V)
  arma::mat H_V = arma::solve(arma::trimatl(L_lower), Phi_V);

  // 7. Estimator: h^{-d}/n * H_0^T (H_V 1_{N_V})
  double result = arma::dot(H_0, arma::sum(H_V, 1)) /
                  (static_cast<double>(n) * std::pow(h, d));
  return result;
}
