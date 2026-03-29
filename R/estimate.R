#' Rerandomization Estimator
#'
#' Estimates the ATE using either Difference-in-Means or Lin's (2013) Regression Adjustment.
#' Supports both conservative estimation (using observed data) and theoretical variance calculation (using potential outcomes).
#'
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment (0 = control, 1 = treatment).
#' @param X Numeric matrix; n x K covariate matrix. Required for Lin method and R2 estimation.
#' @param method Character scalar; "dim" for Difference-in-Means, "lin" for Lin's Regression Adjustment.
#' @param p_accept Numeric scalar; Probability of acceptance (default 1, representing CRE). Used to calculate variance reduction factor.
#' @param theoretical Logical scalar; If TRUE, calculates theoretical variance using Y_full.
#' @param Y_full Numeric matrix; potential outcomes matrix with n rows and 2 columns (Y(0), Y(1)). Required if theoretical = TRUE.
#'
#' @return A list of class "rerand_estimate" containing estimates, standard errors, and diagnostic statistics.
#' @importFrom stats lm var cov predict model.matrix
#' @importFrom car hccm
#' @importFrom checkmate assert_numeric assert_matrix assert_subset assert_choice
#' @export
rerand_estimate <- function(Y_obs, Z, X = NULL,
                            method = c("dim", "lin"),
                            p_accept = 1,
                            theoretical = FALSE,
                            Y_full = NULL) {


  # --- 1. Input Validation ---
  method <- match.arg(method)
  n <- length(Z)
  p_a <- p_accept

  checkmate::assert_numeric(Z, len = n, any.missing = FALSE)
  checkmate::assert_subset(unique(Z), choices = c(0, 1))
  checkmate::assert_number(p_a, lower = 0, upper = 1)

  # Check X
  if (!is.null(X)) {
    checkmate::assert_matrix(X, mode = "numeric", nrows = n, any.missing = FALSE)
  } else if (method == "lin" || p_a < 1) {
    stop("Covariate matrix X is required for method 'lin' or when p_a < 1 (to estimate R2).")
  }

  checkmate::assert_numeric(Y_obs, len = n, any.missing = FALSE)
  if (theoretical) {
    checkmate::assert_matrix(Y_full, mode = "numeric", nrows = n, ncols = 2, any.missing = FALSE)
  }

  # --- 2. Sample Statistics Calculation ---
  sample_stats <- calc_sample_stats(Y_obs = Y_obs, Z = Z, X = X, p_accept = p_a)
  r1 <- sample_stats$r1
  r0 <- sample_stats$r0

  if (method == "dim") {
    dim_res <- est_dim(Y_obs = Y_obs, Z = Z, X = X, p_accept = p_a, stats = sample_stats)
    tau_hat <- dim_res$tau_hat
    se_neyman <- dim_res$se_neyman
    se_ding <- dim_res$se_ding

  } else {
    lin_res <- est_lin(Y_obs = Y_obs, Z = Z, X = X, stats = sample_stats)
    tau_hat <- lin_res$tau_hat
    se_ehw <- lin_res$se_ehw
    fit <- lin_res$fit
  }

  # --- 4. Theoretical Variance (Optional) ---

  # --- 4. Theoretical Variance (Optional) ---

  if (theoretical) {
    pop_stats <- calc_population_stats(
      Y_full = Y_full, X = X, r1 = r1, r0 = r0, 
      p_accept = p_a, method = method
    )
    V_tt_true <- pop_stats$V_tt_true
    R2_true <- pop_stats$R2_true
    tau_true <- pop_stats$tau_true
    se_true <- pop_stats$se_true
  }

  # --- 5. Output Construction ---
  res <- list(
    tau_hat = as.numeric(tau_hat),
    tau_true = if (theoretical) as.numeric(tau_true) else NULL,
    se_neyman = if (method == "dim") as.numeric(se_neyman) else NULL,
    se_ding = if (method == "dim" && !is.null(se_ding)) as.numeric(se_ding) else NULL,
    se_ehw = if (method == "lin") as.numeric(se_ehw) else NULL,
    se_true = if (theoretical) as.numeric(se_true) else NULL,

    stats = list(
      n = n,
      method = method,
      p_accept = p_a,
      R2_hat = if (!is.null(sample_stats$R2_hat)) as.numeric(sample_stats$R2_hat) else NULL,
      R2_true = if (theoretical && !is.null(R2_true)) as.numeric(R2_true) else NULL,
      V_tt_hat_1 = as.numeric(sample_stats$V_tt_hat_1),
      V_tt_hat_2 = if (!is.null(sample_stats$V_tt_hat_2)) as.numeric(sample_stats$V_tt_hat_2) else NULL,
      V_tt_true = if (theoretical) as.numeric(V_tt_true) else NULL,
      fit = if (method == "lin") fit else NULL
    )
  )

  return(res)
}
