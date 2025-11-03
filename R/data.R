#' Centered Covariate Generation Function
#'
#' Generates n observations of K covariates drawn from a multivariate normal distribution with zero mean and specified covariance structure.
#' @param n Interger scalar; number of observations. n > 2.
#' @param K Integer scalar; number of covariates. K > 0.
#' @param Sigma_X Numeric matrix; K × K covariance matrix for X and must be positive semi definite. Default is identity matrix (independent covariates).
#' @param seed Integer scalar; random seed  reproducibility. Default is NULL (no seed set).
#'
#' @return Numeric matrix; n × K covariate matrix.
#'
#' @examples
#' set.seed(123)
#' X_data <- gen_X(n = 100, K = 5)
#'
#' @export
gen_X <- function(n, K, Sigma_X = diag(K), seed = NULL) {

  # Check inputs
  checkmate::assert_count(n)
  checkmate::assert_count(K)
  checkmate::assert_matrix(Sigma_X, mode = "numeric", nrows = K, ncols = K, any.missing = FALSE)
  # still need to check psd
  checkmate::assert_count(seed)

  as.integer(n)
  as.integer(K)

  checkmate::assert_integer(n, lower = 2, any.missing = FALSE)
  checkmate::assert_integer(K, lower = 1, any.missing = FALSE)


  # Set parameters
  if (!is.null(seed)) set.seed(seed)


  # Generate X ~ N(0, Sigma)
  X <- MASS::mvrnorm(n = n, mu = rep(0, K), Sigma_X = Sigma_X)
  return(X)
}



#' Beta Coefficient Generation Function
#'
#' Constructs a beta coefficient vector for covariates X such that the variance explained by the linear predictor Xβ matches a specified target R² value.
#' @param Sigma_X Numeric matrix; K × K covariance matrix for X.
#' @param R2 Numeric scalar; target R² value in [0, 1).
#' @param Sigma_e Numeric scalar; standard deviation of the error term. Default is 1.
#' @param v Numeric vector; unit direction vector for beta. If NULL, a random direction is generated. Default is NULL.
#' @param seed Integer scalar; random seed for reproducibility. Default is NULL (no seed set).
#'
#' @return Numeric vector; beta coefficient vector of length K.
#'
#' @examples
#' set.seed(123)
#' beta <- gen_beta(Sigma_X = diag(5), R2 = 0.5)
#'
#' @export
gen_beta <- function(Sigma_X, R2, Sigma_e = 1, v = NULL, seed = NULL) {

  # Check inputs
  checkmate::assert_matrix(Sigma_X, mode = "numeric", any.missing = FALSE)
  checkmate::assert_numeric(R2, lower = 0, upper = 1, any.missing = FALSE, len = 1)
  checkmate::assert_true(R2 < 1)
  checkmate::assert_numeric(Sigma_e, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_true(Sigma_e > 0)
  checkmate::assert_nu(v, any.missing = FALSE, len = ncol(Sigma_X))
  checkmate::assert_count(seed)

  # Set parameters
  K <- ncol(Sigma_X)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(v)) v <- rnorm(K)
  v <- v / sqrt(sum(v^2))

  # Compute target variance scale
  c_target <- sqrt(R2 / (1 - R2)) * Sigma_e
  scale_fac <- c_target / sqrt(t(v) %*% Sigma_X %*% v)
  beta <- as.numeric(scale_fac * v)
  beta
}



#' Potential Outcomes Generation Function
#'
#' Generates potential outcomes Y(0) and Y(1) based on a linear model with covariates X, coefficients beta, treatment effect tau, and error term. If a treatment assignment vector Z is provided, it also computes the observed outcomes Y.
#' @param X Numeric matrix; n × K covariate matrix.
#' @param beta Numeric vector; coefficient vector of length K.
#' @param tau Numeric scalar; average treatment effect. Default is 1.
#' @param mu Numeric scalar; intercept term. Default is 0.
#' @param Sigma_e Numeric scalar; standard deviation of the error term. Default is 1.
#' @param Sigma_tau Numeric scalar; standard deviation of individual treatment effects. Default is 0 (constant treatment effect).
#' @param seed Integer scalar; random seed for reproducibility. Default is NULL (no seed set).
#'
#' @return Numeric matrix; n × 2 matrix of potential outcomes with first column Y(0) and second column Y(1).
#'
#' @examples
#' set.seed(123)
#' X <- gen_X(n = 100, K = 5)
#' beta <- gen_beta(Sigma_X = diag(5), R2 = 0.5)
#' Y <- gen_Y(X = X, beta = beta)
#'
#' @export
gen_Y <- function(X, beta, tau = 1, mu = 0, Sigma_e = 1, Sigma_tau = 0, seed = NULL) {

  # Check inputs
  checkmate::assert_matrix(X, mode = "numeric", any.missing = FALSE)
  checkmate::assert_numeric(beta, any.missing = FALSE, len = ncol(X))
  checkmate::assert_numeric(tau, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(mu, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(Sigma_e, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_true(Sigma_e > 0)
  checkmate::assert_numeric(Sigma_tau, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_count(seed)

  # Set parameters
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  lin <- as.numeric(X %*% beta)
  e   <- rnorm(n, 0, Sigma_e)

  if (Sigma_tau > 0) {
    tau <- rnorm(n, mean = tau, sd = Sigma_tau)
  }

  Y0 <- mu + lin + e
  Y1 <- mu + tau + lin + e
  Y <- matrix(cbind(Y0, Y1), ncol = 2)
  colnames(Y) <- c("Y0", "Y1")

  return(Y)
}






#
#
#
# X <- gen_X(n = 100, K = 5, seed = 123)
# beta <- gen_beta(Sigma = diag(5), R2 = 0.5, seed = 123)
# Y <- gen_Y(X = X$X, beta = beta)
#
#
# Z <- ReM(X = X$X, n_1 = 50, p_a = 0.1, seed = 123)$Z
#
# est_unadjusted(Y = Y$Y1 * Z + Y$Y0 * (1 - Z), Z = Z)
# est_adjusted(Y = Y$Y1 * Z + Y$Y0 * (1 - Z), Z = Z, X = X$X)




