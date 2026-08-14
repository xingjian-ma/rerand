#' Estimate an average treatment effect from a formula and data frame
#'
#' @param formula A two-sided formula. The first right-hand-side term is the
#'   untransformed treatment column; subsequent terms are covariates.
#' @param data Data frame containing the outcome, treatment, and covariates.
#' @param design Optional `rerand_design_result`. For difference-in-means it
#'   supplies the assignment, balance covariates, and rerandomization criterion.
#' @param method Either `"dim"` or `"lin"`.
#' @param treated Treated level for a two-level factor or character treatment.
#' @param accept_prob Acceptance probability when estimating difference-in-means
#'   without `design`. Mutually exclusive with `threshold`.
#' @param threshold Direct Mahalanobis threshold when estimating
#'   difference-in-means without `design`.
#' @param se_type Standard error for difference-in-means: `"neyman"` or
#'   `"ding"`. The default is `"neyman"`.
#' @return An object of class `rerand_estimate_result`.
#' @export
rerand_estimate <- function(formula, data, design = NULL,
                            method = c("dim", "lin"), treated = NULL,
                            accept_prob = NULL, threshold = NULL,
                            se_type = NULL) {
  method <- match.arg(method)
  parsed <- .rerand_prepare_estimate_formula(formula, data, treated = treated)
  Z <- if (is.null(design)) parsed$Z else .rerand_assignment_from_design(design, data)
  if (!identical(as.numeric(parsed$Z), as.numeric(Z))) {
    stop("The treatment column in data does not match the design assignment.",
         call. = FALSE)
  }

  if (method == "lin") {
    if (!parsed$has_covariates) {
      stop("method = 'lin' requires at least one covariate in formula.",
           call. = FALSE)
    }
    if (!is.null(accept_prob) || !is.null(threshold)) {
      stop("Lin estimation does not use accept_prob or threshold.", call. = FALSE)
    }
    X <- .rerand_prepare_estimation_covariates(data, parsed$covariate_formula)
    return(.rerand_estimate_matrix(
      Y_obs = parsed$Y_obs, Z = Z, X = X, method = "lin", se_type = se_type,
      formula = formula, treatment_name = parsed$treatment_name,
      design = design
    ))
  }

  if (!is.null(design)) {
    if (!is.null(accept_prob) || !is.null(threshold)) {
      stop("When design is supplied, its rerandomization criterion is used.",
           call. = FALSE)
    }
    if (parsed$has_covariates) {
      warning(
        "Formula covariates are ignored for method = 'dim' when design is supplied.",
        call. = FALSE
      )
    }
    return(.rerand_estimate_matrix(
      Y_obs = parsed$Y_obs, Z = Z, X = design$spec$X, method = "dim",
      criterion = design$spec$criterion, se_type = se_type, formula = formula,
      treatment_name = parsed$treatment_name, design = design
    ))
  }

  X <- .rerand_prepare_estimation_covariates(data, parsed$covariate_formula)
  .rerand_estimate_matrix(
    Y_obs = parsed$Y_obs, Z = Z, X = X, method = "dim",
    accept_prob = accept_prob, threshold = threshold, se_type = se_type,
    formula = formula, treatment_name = parsed$treatment_name
  )
}

#' Estimate an average treatment effect from vectors and a covariate matrix
#'
#' @param Y_obs Numeric observed-outcome vector.
#' @param Z Binary treatment-assignment vector.
#' @param X Optional numeric covariate matrix.
#' @param method Either `"dim"` or `"lin"`.
#' @param accept_prob Acceptance probability for difference-in-means.
#' @param threshold Direct Mahalanobis threshold for difference-in-means.
#' @param se_type Standard error for difference-in-means: `"neyman"` or
#'   `"ding"`.
#' @return An object of class `rerand_estimate_result`.
#' @export
rerand_estimate_matrix <- function(Y_obs, Z, X = NULL,
                                   method = c("dim", "lin"),
                                   accept_prob = NULL, threshold = NULL,
                                   se_type = NULL) {
  .rerand_estimate_matrix(
    Y_obs = Y_obs, Z = Z, X = X, method = method,
    accept_prob = accept_prob, threshold = threshold, se_type = se_type
  )
}

.rerand_estimate_matrix <- function(Y_obs, Z, X = NULL,
                                    method = c("dim", "lin"),
                                    accept_prob = NULL, threshold = NULL,
                                    se_type = NULL, criterion = NULL,
                                    formula = NULL, treatment_name = NULL,
                                    design = NULL) {
  method <- match.arg(method)
  .rerand_assert_finite_numeric(Y_obs, "Y_obs")
  assignment <- .rerand_validate_assignment(Z, length(Y_obs), min_group_size = 2L)
  if (!is.null(X)) {
    X <- .rerand_validate_matrix(X, n = length(Y_obs))
  }
  if (method == "lin") {
    if (is.null(X)) {
      stop("X is required for method = 'lin'.", call. = FALSE)
    }
    if (!is.null(accept_prob) || !is.null(threshold) || !is.null(criterion)) {
      stop("Lin estimation does not use a rerandomization criterion.",
           call. = FALSE)
    }
    estimate <- .estimate_lin(as.numeric(Y_obs), assignment$Z, X)
    return(.rerand_estimate_result(
      estimate = estimate, method = method, sample_stats = NULL,
      criterion = NULL, se_type = se_type, formula = formula,
      treatment_name = treatment_name, design = design
    ))
  }

  if (is.null(criterion)) {
    if (is.null(X)) {
      if (!is.null(threshold) || is.null(accept_prob) || accept_prob != 1) {
        stop("Without covariates, accept_prob must be explicitly set to 1.",
             call. = FALSE)
      }
      criterion <- .rerand_complete_randomization_criterion()
    } else {
      rank <- .rerand_whiten_covariates(X)$effective_rank
      criterion <- .rerand_resolve_criterion(
        accept_prob = accept_prob, threshold = threshold, K = rank,
        require_criterion = TRUE
      )
    }
  }
  sample_stats <- .calc_sample_stats(
    Y_obs = as.numeric(Y_obs), Z = assignment$Z, X = X, criterion = criterion
  )
  .rerand_estimate_result(
    estimate = .estimate_dim(sample_stats), method = method,
    sample_stats = sample_stats, criterion = criterion, se_type = se_type,
    formula = formula, treatment_name = treatment_name, design = design
  )
}

#' Calculate population-level rerandomization benchmarks
#'
#' @param Y_full Numeric matrix with columns `Y(0)` and `Y(1)`.
#' @param design Optional `rerand_design_result` supplying the covariates,
#'   treatment-group size, and criterion.
#' @param X Optional numeric covariate matrix when `design` is omitted.
#' @param n_treat Number of treated units when `design` is omitted.
#' @param accept_prob Acceptance probability when `design` is omitted.
#' @param threshold Direct Mahalanobis threshold when `design` is omitted.
#' @return A list of population-level treatment-effect and variance quantities.
#' @export
rerand_population_stats <- function(Y_full, design = NULL, X = NULL,
                                    n_treat = NULL, accept_prob = NULL,
                                    threshold = NULL) {
  if (!is.matrix(Y_full) || !is.numeric(Y_full) || any(!is.finite(Y_full)) ||
      ncol(Y_full) != 2L || nrow(Y_full) < 2L) {
    stop("Y_full must be a finite numeric matrix with two columns.",
         call. = FALSE)
  }
  if (!is.null(design)) {
    if (!inherits(design, "rerand_design_result")) {
      stop("design must be a rerand_design_result object.", call. = FALSE)
    }
    if (!is.null(X) || !is.null(n_treat) || !is.null(accept_prob) ||
        !is.null(threshold)) {
      stop("design already supplies X, n_treat, and the criterion.", call. = FALSE)
    }
    X <- design$spec$X
    n_treat <- design$n_treat
    criterion <- design$spec$criterion
  } else if (is.null(X)) {
    n_treat <- .rerand_validate_n_treat(n_treat, nrow(Y_full))
    if (!is.null(threshold) || is.null(accept_prob) || accept_prob != 1) {
      stop("Without covariates, accept_prob must be explicitly set to 1.",
           call. = FALSE)
    }
    criterion <- .rerand_complete_randomization_criterion()
  } else {
    X <- .rerand_validate_matrix(X, n = nrow(Y_full))
    n_treat <- .rerand_validate_n_treat(n_treat, nrow(Y_full))
    criterion <- .rerand_resolve_criterion(
      accept_prob = accept_prob, threshold = threshold,
      K = .rerand_whiten_covariates(X)$effective_rank,
      require_criterion = TRUE
    )
  }
  .calc_population_stats(Y_full, X, n_treat, criterion)
}

.rerand_complete_randomization_criterion <- function() {
  list(
    type = "probability", accept_prob = 1, threshold = Inf, K = 0L,
    acceptance_mass = 1, v_K_a = 1
  )
}

.rerand_estimate_result <- function(estimate, method, sample_stats, criterion,
                                    se_type, formula, treatment_name, design) {
  if (is.null(se_type)) {
    se_type <- if (method == "dim") "neyman" else "hc2"
  }
  valid_se <- if (method == "dim") c("neyman", "ding") else "hc2"
  if (length(se_type) != 1L || !is.character(se_type) || !se_type %in% valid_se) {
    stop(sprintf("se_type must be one of %s.", paste(valid_se, collapse = ", ")),
         call. = FALSE)
  }
  if (method == "dim" && se_type == "ding" && is.null(estimate$se_ding)) {
    stop("se_type = 'ding' requires covariates.", call. = FALSE)
  }
  standard_error <- if (method == "dim") {
    if (se_type == "neyman") estimate$se_neyman else estimate$se_ding
  } else {
    estimate$se_ehw
  }
  unadjusted_standard_error <- if (method == "dim") {
    variance <- if (se_type == "neyman") {
      sample_stats$V_tt_hat_1
    } else {
      sample_stats$V_tt_hat_2
    }
    sqrt(variance) / sqrt(sample_stats$n)
  } else {
    estimate$se_ehw
  }
  result <- list(
    tau_hat = as.numeric(estimate$tau_hat), method = method,
    standard_error = as.numeric(standard_error),
    unadjusted_standard_error = as.numeric(unadjusted_standard_error),
    se_type = se_type,
    se_neyman = if (method == "dim") estimate$se_neyman else NULL,
    se_ding = if (method == "dim") estimate$se_ding else NULL,
    se_ehw = if (method == "lin") estimate$se_ehw else NULL,
    fit = if (method == "lin") estimate$fit else NULL,
    sample_stats = sample_stats,
    criterion_type = if (is.null(criterion)) NULL else criterion$type,
    accept_prob = if (is.null(criterion)) NULL else criterion$accept_prob,
    threshold = if (is.null(criterion)) NULL else criterion$threshold,
    formula = formula, treatment_name = treatment_name, design = design
  )
  class(result) <- c("rerand_estimate_result", "list")
  result
}

#' @export
print.rerand_estimate_result <- function(x, ...) {
  cat("Rerandomization estimate\n")
  cat("  method:        ", x$method, "\n", sep = "")
  cat("  estimate:      ", format(x$tau_hat), "\n", sep = "")
  cat("  ", toupper(x$se_type), " SE: ", format(x$standard_error), "\n",
      sep = "")
  if (!is.null(x$criterion_type)) {
    cat("  criterion:     ", x$criterion_type, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.rerand_estimate_result <- function(object, ...) {
  summary <- object[c(
    "tau_hat", "standard_error", "se_type", "method", "se_neyman",
    "se_ding", "se_ehw", "criterion_type", "accept_prob", "threshold"
  )]
  class(summary) <- c("summary.rerand_estimate_result", "list")
  summary
}

#' @export
print.summary.rerand_estimate_result <- function(x, ...) {
  print.rerand_estimate_result(x, ...)
  invisible(x)
}

#' @export
coef.rerand_estimate_result <- function(object, ...) {
  stats::setNames(object$tau_hat, "treatment")
}

#' @export
vcov.rerand_estimate_result <- function(object, ...) {
  matrix(object$standard_error^2, nrow = 1L, ncol = 1L,
         dimnames = list("treatment", "treatment"))
}

#' @export
confint.rerand_estimate_result <- function(object, parm, level = 0.95, ...) {
  if (missing(parm)) {
    parm <- "treatment"
  }
  if (!identical(parm, "treatment")) {
    stop("Only the treatment coefficient is available.", call. = FALSE)
  }
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be strictly between 0 and 1.", call. = FALSE)
  }
  alpha <- (1 + level) / 2
  uses_rerandomization <- object$method == "dim" &&
    !is.null(object$sample_stats) &&
    !is.null(object$sample_stats$R2_hat) &&
    object$sample_stats$acceptance_mass < 1
  if (uses_rerandomization) {
    quantile_arguments <- list(
      R2 = object$sample_stats$R2_hat,
      K = object$sample_stats$K,
      alpha = alpha
    )
    if (object$sample_stats$criterion_type == "probability") {
      quantile_arguments$accept_prob <- object$sample_stats$accept_prob
    } else {
      quantile_arguments$threshold <- object$sample_stats$threshold
    }
    critical <- do.call(get_quantile, quantile_arguments)
  } else {
    critical <- stats::qnorm(alpha)
  }
  matrix(object$tau_hat + c(-1, 1) * critical * object$unadjusted_standard_error,
         nrow = 1L, dimnames = list("treatment", c("lower", "upper")))
}
