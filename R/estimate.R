#' Estimate an average treatment effect from complete data and an assignment
#'
#' @description
#' Formula mode derives the estimator from the right-hand side: `Y ~ Z` is
#' difference-in-means, `Y ~ Z + x` is additive ANCOVA, and
#' `Y ~ Z * (x1 + x2)` is Lin's fully interacted estimator. Selector mode uses
#' string column names and an explicit `estimator`.
#'
#' @param data Complete data frame containing the observed outcome, treatment,
#'   and covariates.
#' @param assignment A `rerand_assignment` object created by [rerand_assign()].
#' @param formula Optional canonical two-sided formula.
#' @param outcome Outcome column name for selector mode.
#' @param treatment Treatment column name for selector mode.
#' @param covariates Optional adjustment covariate column names for selector
#'   mode.
#' @param estimator Selector-mode estimator: `"dim"`, `"ancova"`, or `"lin"`.
#' @param treated Treated level for a nonnumeric treatment column.
#' @param se_type DIM standard error: `"neyman"` or `"ding"`; regression
#'   estimators use `"hc2"`.
#' @return An object of class `rerand_estimate`.
#' @export
rerand_estimate <- function(data, assignment, formula = NULL, outcome = NULL,
                            treatment = NULL, covariates = NULL,
                            estimator = NULL, treated = NULL, se_type = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  selector_supplied <- !is.null(outcome) || !is.null(treatment) ||
    !is.null(covariates) || !is.null(estimator)
  if (!is.null(formula)) {
    if (selector_supplied) {
      warning("formula takes priority; selector arguments and estimator were ignored.",
              call. = FALSE)
    }
    parsed <- .rerand_parse_formula(formula, data, treated = treated)
  } else {
    parsed <- .rerand_parse_selectors(
      data = data, outcome = outcome, treatment = treatment,
      covariates = covariates, estimator = estimator, treated = treated
    )
  }
  assignment_z <- .rerand_assignment_vector(assignment, data)
  if (!isTRUE(all.equal(as.numeric(parsed$Z), as.numeric(assignment_z)))) {
    stop("The treatment column in data does not match the assignment.",
         call. = FALSE)
  }
  design <- assignment$design
  criterion <- design$criterion
  design_X <- .rerand_design_covariates(assignment, data)
  analysis_X <- .rerand_prepare_estimation_covariates(
    data, parsed$covariate_formula
  )
  if (parsed$estimator == "dim") {
    if (parsed$has_covariates && is.null(formula)) {
      warning("covariates are ignored for estimator = 'dim'; design covariates are used.",
              call. = FALSE)
    }
    sample_stats <- .calc_sample_stats(
      Y_obs = parsed$Y_obs, Z = assignment_z, X = design_X,
      criterion = criterion
    )
    estimate <- .estimate_dim(sample_stats)
  } else if (parsed$estimator == "ancova") {
    estimate <- .estimate_ancova(parsed$Y_obs, assignment_z, analysis_X)
    sample_stats <- NULL
  } else {
    estimate <- .estimate_lin(parsed$Y_obs, assignment_z, analysis_X)
    sample_stats <- NULL
  }
  .rerand_estimate_result(
    estimate = estimate, estimator = parsed$estimator,
    sample_stats = sample_stats, criterion = criterion, se_type = se_type,
    formula = if (is.null(formula)) NULL else formula,
    outcome_name = parsed$outcome_name, treatment_name = parsed$treatment_name,
    assignment = assignment, data = data,
    analysis_covariates = if (is.null(analysis_X)) character() else colnames(analysis_X)
  )
}

.rerand_estimate_result <- function(estimate, estimator, sample_stats, criterion,
                                    se_type, formula, outcome_name,
                                    treatment_name, assignment, data,
                                    analysis_covariates) {
  if (is.null(se_type)) {
    se_type <- if (estimator == "dim") "neyman" else "hc2"
  }
  valid_se <- if (estimator == "dim") c("neyman", "ding") else "hc2"
  if (length(se_type) != 1L || !is.character(se_type) ||
      !se_type %in% valid_se) {
    stop(sprintf("se_type must be one of %s.", paste(valid_se, collapse = ", ")),
         call. = FALSE)
  }
  if (estimator == "dim" && se_type == "ding" &&
      is.null(estimate$se_ding)) {
    stop("se_type = 'ding' requires covariates.", call. = FALSE)
  }
  standard_error <- if (estimator == "dim") {
    if (se_type == "neyman") estimate$se_neyman else estimate$se_ding
  } else {
    estimate$se_ehw
  }
  unadjusted_standard_error <- if (estimator == "dim") {
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
    estimate = as.numeric(estimate$tau_hat),
    estimator = estimator,
    standard_error = as.numeric(standard_error),
    unadjusted_standard_error = as.numeric(unadjusted_standard_error),
    se_type = se_type,
    se_neyman = if (estimator == "dim") estimate$se_neyman else NULL,
    se_ding = if (estimator == "dim") estimate$se_ding else NULL,
    se_hc2 = if (estimator != "dim") estimate$se_ehw else NULL,
    fit = if (estimator != "dim") estimate$fit else NULL,
    sample_stats = sample_stats,
    criterion_type = criterion$type,
    accept_prob = criterion$accept_prob,
    threshold = criterion$threshold,
    design_method = assignment$design_method,
    formula = formula,
    outcome_name = outcome_name,
    treatment_name = treatment_name,
    analysis_covariates = analysis_covariates,
    assignment = assignment,
    data = data,
    n = nrow(data)
  )
  class(result) <- c("rerand_estimate", "list")
  result
}

#' @export
print.rerand_estimate <- function(x, ...) {
  cat("Rerandomization estimate\n")
  cat("  estimator:     ", x$estimator, "\n", sep = "")
  cat("  estimate:      ", format(x$estimate), "\n", sep = "")
  cat("  ", toupper(x$se_type), " SE: ", format(x$standard_error), "\n",
      sep = "")
  cat("  design method: ", x$design_method, "\n", sep = "")
  invisible(x)
}

#' @export
summary.rerand_estimate <- function(object, ...) {
  result <- object[c(
    "estimate", "standard_error", "se_type", "estimator",
    "design_method", "criterion_type", "accept_prob", "threshold"
  )]
  class(result) <- c("summary.rerand_estimate", "list")
  result
}

#' @export
print.summary.rerand_estimate <- function(x, ...) {
  print.rerand_estimate(x, ...)
  invisible(x)
}

#' @export
coef.rerand_estimate <- function(object, ...) {
  stats::setNames(object$estimate, "treatment")
}

#' @export
vcov.rerand_estimate <- function(object, ...) {
  matrix(object$standard_error^2, nrow = 1L, ncol = 1L,
         dimnames = list("treatment", "treatment"))
}
