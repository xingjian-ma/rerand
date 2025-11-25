#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// partial Fisher–Yates: shuffle first n1 positions only
inline void sample_Z(arma::vec& Z, arma::uvec& perm, int n, int n1) {
  // reset Z
  Z.zeros();
  
  // partial shuffle perm in-place
  for (int i = 0; i < n1; i++) {
    int j = i + (int)std::floor(R::runif(0.0, 1.0) * (n - i));
    std::swap(perm[i], perm[j]);
    Z[ perm[i] ] = 1.0;  // first n1 indices treated
  }
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
  arma::vec Z0(n);
  arma::vec xbar1(K), xbar0(K), diff(K);
  arma::uvec perm = arma::regspace<arma::uvec>(0, n-1);
  
  const double scale = (double)n1 * (double)n0 / (double)n;
  double M = arma::datum::inf;
  int tries = max_tries;
  bool accepted = false;
  
  // placeholders for output
  arma::vec Y_obs(n);
  
  
  for (int t = 0; t < max_tries; t++) {
    
    // draw assignment with exactly n1 treated
    sample_Z(Z, perm, n, n1);
    Z0 = 1.0 - Z;
    
    // group means
    xbar1 = (X.t() * Z) / (double)n1;
    xbar0 = (X.t() * Z0) / (double)n0;
    diff = xbar1 - xbar0;
    
    // Mahalanobis distance
    M = scale * arma::dot( diff, S_inv * diff );
    
    
    if (M <= a) {
      tries = t + 1;
      accepted = true;
      
      // observed outcomes under accepted Z
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
