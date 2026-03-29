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
Rcpp::List design_cpp(const arma::mat& X,
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
      break;
    }
  }


  if (!accepted) {
    warning("Maximum tries exceeded without reaching threshold. Returning last assignment anyway.");
  }


  return Rcpp::List::create(
    _["Z"] = Rcpp::NumericVector(Z.begin(), Z.end()),
    _["tries"] = tries,
    _["M"] = M,
    _["accepted"] = accepted
  );
}

// [[Rcpp::export]]
double get_quantile_cpp(double R2, int K, double p_a, double alpha, int n_sim) {

  double a = R::qchisq(p_a, (double)K, 1, 0);
  double target_prob = 2.0 * alpha - 1.0;

  arma::vec abs_dist(n_sim);
  int count = 0;

  while (count < n_sim) {
    double d1 = R::rnorm(0, 1);
    double w = (K > 1) ? R::rchisq((double)K - 1) : 0;

    if ((d1 * d1 + w) <= a) {
      double val = std::sqrt(R2) * d1 + std::sqrt(1.0 - R2) * R::rnorm(0, 1);
      abs_dist(count) = std::abs(val);
      count++;
    }
  }

  arma::vec sorted_abs_dist = arma::sort(abs_dist);

  int target_idx = (int)(target_prob * n_sim) - 1;
  if (target_idx < 0) target_idx = 0;

  return sorted_abs_dist(target_idx);
}
