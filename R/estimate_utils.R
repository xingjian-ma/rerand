#' Difference-in-Means Estimator under Rerandomization
#'
#' Computes the difference-in-means point estimate and conservative standard errors,
#' including rerandomization adjustment when covariates are provided.
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment with length n (0 = control, 1 = treatment).
#' @param X Numeric matrix; n x K covariate matrix. Optional. Default is NULL.
#' @param p_accept Numeric scalar; acceptance probability.
#'
#' @return A list containing:
#' \describe{
#' \item{tau_hat}{Numeric scalar; estimated ATE.}
#' \item{se_neyman}{Numeric scalar; standard error of the estimate using Neyman's conservative variance estimator.}
#' \item{se_ding}{Numeric scalar; standard error of the estimate using Ding's improved variance estimator.}
#' \item{R2_hat}{Numeric scalar; estimated variance explained by covariates, if X is provided.}
#' \item{V_tt_hat_1}{Numeric scalar; baseline conservative variance component.}
#' \item{V_tt_hat_2}{Numeric scalar; adjusted conservative variance component, if X is provided.}
#' }
#'
#' @examples
#' res <- est_dim(Y_obs = rnorm(100), Z = rbinom(100, 1, 0.5), p_accept = 1)
#'
#' @keywords internal
calc_sample_stats <- function(Y_obs, Z, X = NULL, p_accept = 1) {

  p_a <- p_accept

  Y1 <- Y_obs[Z == 1]
  Y0 <- Y_obs[Z == 0]
  n1 <- length(Y1)
  n0 <- length(Y0)
  n <- length(Y_obs)
  r1 <- n1 / n
  r0 <- n0 / n
  tau_dim <- mean(Y1) - mean(Y0)

  V_tt_hat_1 <- stats::var(Y1) / r1 + stats::var(Y0) / r0

  if (is.null(X)) {
    return(list(
      n = n,
      n1 = n1,
      n0 = n0,
      r1 = r1,
      r0 = r0,
      tau_dim = tau_dim,
      V_tt_hat_1 = as.numeric(V_tt_hat_1),
      V_tt_hat_2 = NULL,
      R2_hat = NULL,
      v_K_a = NULL,
      correction_factor = 1
    ))
  }

  # Rerandomization-adjusted conservative variance
  K <- ncol(X)
  a <- stats::qchisq(p = p_a, df = K)
  v_K_a <- stats::pchisq(q = a, df = K + 2) / p_a

  X1 <- X[Z == 1, , drop = FALSE]
  X0 <- X[Z == 0, , drop = FALSE]

  S_inv <- solve(stats::cov(X))

  S_Y1X <- matrix(stats::cov(Y1, X1), nrow = 1)
  S_Y0X <- matrix(stats::cov(Y0, X0), nrow = 1)

  S_tauX <- as.numeric((S_Y1X - S_Y0X) %*% S_inv %*% t(S_Y1X - S_Y0X))
  V_tt_hat_2 <- as.numeric(V_tt_hat_1 - S_tauX)

  S_Y1_given_X <- as.numeric(S_Y1X %*% S_inv %*% t(S_Y1X))
  S_Y0_given_X <- as.numeric(S_Y0X %*% S_inv %*% t(S_Y0X))

  R2_hat <- as.numeric((S_Y0_given_X / r0 + S_Y1_given_X / r1 - S_tauX) / V_tt_hat_2)
  correction_factor <- as.numeric(1 - (1 - v_K_a) * R2_hat)

  list(
    n = n,
    n1 = n1,
    n0 = n0,
    r1 = r1,
    r0 = r0,
    tau_dim = tau_dim,
    V_tt_hat_1 = as.numeric(V_tt_hat_1),
    V_tt_hat_2 = as.numeric(V_tt_hat_2),
    R2_hat = as.numeric(R2_hat),
    v_K_a = as.numeric(v_K_a),
    correction_factor = as.numeric(correction_factor)
  )
}


#' Difference-in-Means Estimator under Rerandomization
#'
#' Computes the difference-in-means point estimate and conservative standard errors,
#' including rerandomization adjustment when covariates are provided.
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment with length n (0 = control, 1 = treatment).
#' @param X Numeric matrix; n x K covariate matrix. Optional. Default is NULL.
#' @param p_accept Numeric scalar; acceptance probability.
#' @param ... Additional arguments for forward compatibility.
#'
#' @return A list containing:
#' \describe{
#' \item{tau_hat}{Numeric scalar; estimated ATE.}
#' \item{se_neyman}{Numeric scalar; standard error of the estimate using Neyman's conservative variance estimator.}
#' \item{se_ding}{Numeric scalar; standard error of the estimate using Ding's improved variance estimator.}
#' \item{R2_hat}{Numeric scalar; estimated variance explained by covariates, if X is provided.}
#' \item{V_tt_hat_1}{Numeric scalar; baseline conservative variance component.}
#' \item{V_tt_hat_2}{Numeric scalar; adjusted conservative variance component, if X is provided.}
#' }
#'
#' @examples
#' res <- est_dim(Y_obs = rnorm(100), Z = rbinom(100, 1, 0.5), p_accept = 1)
#'
#' @keywords internal
#' @export
est_dim <- function(Y_obs, Z, X = NULL, p_accept = 1, ...) {

  if (is.null(stats)) {
    sample_stats <- calc_sample_stats(Y_obs = Y_obs, Z = Z, X = X, p_accept = p_accept)
  } else {
    sample_stats <- stats
  }

  tau_hat <- sample_stats$tau_dim

  if (is.null(X)) {
    se_neyman <- sqrt(sample_stats$V_tt_hat_1) / sqrt(sample_stats$n)
    se_ding <- NULL
  } else {
    se_neyman <- sqrt(sample_stats$V_tt_hat_1 * sample_stats$correction_factor) / sqrt(sample_stats$n)
    se_ding <- sqrt(sample_stats$V_tt_hat_2 * sample_stats$correction_factor) / sqrt(sample_stats$n)
  }

  return(list(
    tau_hat = as.numeric(tau_hat),
    se_neyman = as.numeric(se_neyman),
    se_ding = if (is.null(se_ding)) NULL else as.numeric(se_ding),
    R2_hat = if (is.null(sample_stats$R2_hat)) NULL else as.numeric(sample_stats$R2_hat),
    V_tt_hat_1 = as.numeric(sample_stats$V_tt_hat_1),
    V_tt_hat_2 = if (is.null(sample_stats$V_tt_hat_2)) NULL else as.numeric(sample_stats$V_tt_hat_2)
  ))
}


#' Lin (2013) Regression-Adjusted Estimator
#'
#' Computes the Lin (2013) fully interacted regression-adjusted ATE estimator
#' with heteroskedasticity-consistent standard errors.
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment with length n (0 = control, 1 = treatment).
#' @param X Numeric matrix; n x K covariate matrix.
#' @param ... Additional arguments for forward compatibility.
#'
#' @return A list containing:
#' \describe{
#' \item{tau_hat}{Numeric scalar; estimated ATE.}
#' \item{se_ehw}{Numeric scalar; heteroskedasticity-consistent standard error.}
#' \item{fit}{lm object; fitted linear model.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(200), nrow = 100, ncol = 2)
#' res <- est_lin(Y_obs = rnorm(100), Z = rbinom(100, 1, 0.5), X = X)
#'
#' @keywords internal
#' @export
est_lin <- function(Y_obs, Z, X, ...) {

  if (is.null(stats)) {
    stats <- calc_sample_stats(Y_obs = Y_obs, Z = Z, X = X, p_accept = 1)
  }

  X_means <- colMeans(X)
  Xc <- sweep(X, 2, X_means, FUN = "-")

  # Build data.frame with explicit column names
  df <- data.frame(Y = as.numeric(Y_obs), Z = as.numeric(Z), Xc)
  colnames(df) <- c("Y", "Z", paste0("X.", seq_len(ncol(Xc))))

  # Lin model: Y ~ Z + Xc + Z:Xc (fully interacted)
  rhs <- paste0(
    "Z + ",
    paste0("X.", seq_len(ncol(Xc)), collapse = " + "),
    " + ",
    paste0("Z:X.", seq_len(ncol(Xc)), collapse = " + ")
  )
  fm <- stats::as.formula(paste0("Y ~ ", rhs))

  fit <- stats::lm(fm, data = df)

  tau_hat <- summary(fit)$coefficients["Z", "Estimate"]

  # Robust (heteroskedasticity-consistent) SE
  se_ehw <- sqrt(car::hccm(fit, type = "hc2")[2, 2])


  list(tau_hat = tau_hat, se_ehw = se_ehw, fit = fit)
}


#' Calculate Population-Level Statistics (Theoretical Variance)
#'
#' Computes theoretical variance, R-squared, and true ATE from potential outcomes.
#' Used for benchmarking and comparison with conservative sample estimates.
#'
#' @param Y_full Numeric matrix; n x 2 potential outcomes matrix with columns Y(0), Y(1).
#' @param X Numeric matrix; n x K covariate matrix. Optional. Default is NULL.
#' @param r1 Numeric scalar; treatment assignment ratio.
#' @param r0 Numeric scalar; control assignment ratio.
#' @param p_accept Numeric scalar; rerandomization acceptance probability.
#'
#' @return A list containing:
#' \describe{
#' \item{V_tt_true}{Numeric scalar; true sampling variance.}
#' \item{R2_true}{Numeric scalar; true R-squared (proportion of variance explained by covariates).}
#' \item{tau_true}{Numeric scalar; true ATE from potential outcomes.}
#' \item{se_true}{Numeric scalar; true standard error under rerandomization or CRE.}
#' }
#'
#' @keywords internal
#' @noRd
calc_population_stats <- function(Y_full, X = NULL, r1, r0, p_accept = 1) {

  Y0 <- Y_full[, 1]
  Y1 <- Y_full[, 2]
  tau <- Y1 - Y0
  tau_true <- mean(tau)
  n <- length(Y0)

  S_Y0 <- stats::var(Y0)
  S_Y1 <- stats::var(Y1)
  S_tau <- stats::var(tau)

  V_tt_true <- S_Y1 / r1 + S_Y0 / r0 - S_tau

  if (is.null(X)) {
    R2_true <- NULL
    se_true <- sqrt(V_tt_true) / sqrt(n)
  } else {
    S_Y0_given_X <- stats::var(fitted(stats::lm(Y0 ~ X)))
    S_Y1_given_X <- stats::var(fitted(stats::lm(Y1 ~ X)))
    S_tau_given_X <- stats::var(fitted(stats::lm(tau ~ X)))

    R2_true <- (S_Y0_given_X / r0 + S_Y1_given_X / r1 - S_tau_given_X) / V_tt_true

    if (method == "dim") {
      K <- ncol(X)
      a <- stats::qchisq(p = p_accept, df = K)
      v_K_a <- stats::pchisq(q = a, df = K + 2) / p_accept
      se_true <- sqrt((1 - (1 - v_K_a) * R2_true) * V_tt_true) / sqrt(n)
    } else {
      se_true <- sqrt((1 - R2_true) * V_tt_true) / sqrt(n)
    }
  }

  list(
    V_tt_true = as.numeric(V_tt_true),
    R2_true = if (is.null(R2_true)) NULL else as.numeric(R2_true),
    tau_true = as.numeric(tau_true),
    se_true = as.numeric(se_true)
  )
}


#' Generate Quantiles for Estimator under Rerandomization
#'
#' This function simulates the distribution of the estimator under rerandomization to compute quantiles for confidence intervals. It accounts for the variance reduction due to rerandomization and the correlation structure of covariates.
#'
#' @param R2 Numeric scalar; the R-squared value representing the proportion of variance explained by covariates.
#' @param K Integer scalar; the number of covariates used in rerandomization.
#' @param p_accept Numeric scalar; the acceptance probability used in rerandomization.
#' @param alpha Numeric scalar; quantile level of mixed Gaussian distribution (default 0.975 for two-sided 95% CI).
#' @param n_sim Integer scalar; number of simulation points to generate for estimating the quantiles (default 1e5).
#' @param engine Character scalar; computation engine to use, either "R" or "cpp". Default is "cpp".
#'
#' @return Numeric scalar; the estimated quantile of mixed Gaussian distribution.
#' @keywords internal
#' @export
get_quantile <- function(R2, K, p_accept, alpha = 0.975, n_sim = 1e5, engine = "cpp") {

  if (engine == "cpp") {
    quantile <- get_quantile_cpp(R2 = R2, K = K, p_a = p_accept, alpha = alpha, n_sim = n_sim)
  } else {
    # Simulate from the mixed Gaussian distribution
    v_K_a <- stats::pchisq(q = stats::qchisq(p = p_accept, df = K), df = K + 2) / p_accept
    sigma2 <- (1 - (1 - v_K_a) * R2)

    sim_vals <- stats::rnorm(n_sim) * sqrt(sigma2)
    quantile <- stats::quantile(sim_vals, probs = alpha)
  }

  return(quantile)
}
