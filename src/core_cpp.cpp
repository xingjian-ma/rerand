#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// sample n1 indices from {0,...,n-1} without replacement
inline arma::uvec sample_n1(const int n, const int n1) {
  return arma::randperm(n, n1);
}

// [[Rcpp::export]]
Rcpp::List rem_core_cpp(const arma::mat& X,
                        const arma::mat& Y,
                        const int n1,
                        const double a,
                        const int max_tries) {
  
  const int n = X.n_rows;
  const int K = X.n_cols;
  const int n0 = n - n1;
  
  // --- covariance + inverse (inside C++) ---
  arma::mat S = arma::cov(X);
  arma::mat S_inv;
  
  bool ok = arma::inv_sympd(S_inv, S);
  if (!ok) {
    S_inv = arma::pinv(S);
  }
  
  arma::vec Z(n, arma::fill::zeros);
  arma::vec xbar1(K), xbar0(K), diff(K);
  
  double M = arma::datum::inf;
  int tries = max_tries;
  bool accepted = false;
  
  // placeholders for output
  arma::vec Y_obs(n);
  
  for (int t = 0; t < max_tries; t++) {
    
    // draw assignment with exactly n1 treated
    Z.zeros();
    arma::uvec idx1 = sample_n1(n, n1);
    Z.elem(idx1).ones();
    
    arma::uvec idx0 = arma::find(Z == 0.0);
    
    // group means
    xbar1 = arma::mean(X.rows(idx1), 0).t();
    xbar0 = arma::mean(X.rows(idx0), 0).t();
    diff  = xbar1 - xbar0;
    
    // Mahalanobis distance
    M = arma::as_scalar(diff.t() * S_inv * diff) *
      ( (double)n1 * (double)n0 / (double)n );
    
    if (M <= a) {
      tries = t + 1;
      accepted = true;
      
      // observed outcomes under accepted Z
      // Y_obs = Y0*(1-Z) + Y1*Z
      Y_obs = Y.col(0) % (1.0 - Z) + Y.col(1) % Z;
      
      break;
    }
  }
  
  if (!accepted) {
    warning("Maximum tries exceeded without reaching threshold. Returning last assignment anyway.");
    // compute Y_obs for last Z
    Y_obs = Y.col(0) % (1.0 - Z) + Y.col(1) % Z;
  }
  
  return Rcpp::List::create(
    _["Z"] = Rcpp::NumericVector(Z.begin(), Z.end()),
    _["Y_obs"] = Rcpp::NumericVector(Y_obs.begin(), Y_obs.end()),
    _["tries"] = tries,
    _["M"] = M,
    _["accepted"] = accepted
  );
}
