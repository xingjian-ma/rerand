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
#' @return A list containing:
#' \describe{
#' \item{tau_hat}{Numeric scalar; point estimate of ATE.}
#' \item{se_neyman}{Numeric scalar; Neyman SE for method = "dim"; otherwise NULL.}
#' \item{se_ding}{Numeric scalar; Ding SE for method = "dim" when available; otherwise NULL.}
#' \item{se_ehw}{Numeric scalar; HC2 SE for method = "lin"; otherwise NULL.}
#' \item{fit}{lm object; fitted Lin model when method = "lin"; otherwise NULL.}
#' \item{sample_stats}{List; sample-level diagnostics from `calc_sample_stats`.}
#' \item{pop_stats}{List; population-level diagnostics from `calc_population_stats` when `theoretical = TRUE`; otherwise NULL.}
#' }
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
  checkmate::assert_true(p_a > 0)

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

  if (method == "dim") {
    dim_res <- est_dim(Y_obs = Y_obs, Z = Z, X = X, p_accept = p_a, sample_stats = sample_stats)
    tau_hat <- dim_res$tau_hat
    se_neyman <- dim_res$se_neyman
    se_ding <- dim_res$se_ding

  } else {
    lin_res <- est_lin(Y_obs = Y_obs, Z = Z, X = X, sample_stats = sample_stats)
    tau_hat <- lin_res$tau_hat
    se_ehw <- lin_res$se_ehw
    fit <- lin_res$fit
  }

  # --- 4. Theoretical Variance (Optional) ---

  if (theoretical) {
    pop_stats <- calc_population_stats(
      Y_full = Y_full, X = X, n1 = sample_stats$n1, p_accept = p_a
    )
  }

  # --- 5. Output Construction ---
  res <- list(
    tau_hat = as.numeric(tau_hat),
    se_neyman = if (method == "dim") as.numeric(se_neyman) else NULL,
    se_ding = if (method == "dim" && !is.null(se_ding)) as.numeric(se_ding) else NULL,
    se_ehw = if (method == "lin") as.numeric(se_ehw) else NULL,
    fit = if (method == "lin") fit else NULL,
    sample_stats = sample_stats,
    pop_stats = if (theoretical) pop_stats else NULL
  )

  return(res)
}
