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

# Internal data and formula preparation helpers.

.rerand_flatten_addition <- function(expression) {
  if (is.call(expression) && identical(expression[[1L]], as.name("+"))) {
    c(
      .rerand_flatten_addition(expression[[2L]]),
      .rerand_flatten_addition(expression[[3L]])
    )
  } else {
    list(expression)
  }
}

.rerand_unwrap_parentheses <- function(expression) {
  while (is.call(expression) && identical(expression[[1L]], as.name("("))) {
    expression <- expression[[2L]]
  }
  expression
}

.rerand_make_covariate_formula <- function(terms, environment) {
  if (length(terms) == 0L) {
    return(NULL)
  }
  rhs <- Reduce(function(left, right) call("+", left, right), terms)
  stats::as.formula(call("~", rhs), env = environment)
}

.rerand_validate_column_selector <- function(selector, name, data) {
  if (length(selector) != 1L || !is.character(selector) || is.na(selector) ||
      !nzchar(selector) || !selector %in% names(data)) {
    stop(sprintf("%s must name one column in data.", name), call. = FALSE)
  }
  selector
}

.rerand_validate_covariate_selectors <- function(covariates, data) {
  if (is.null(covariates)) {
    return(character())
  }
  if (!is.character(covariates) || anyNA(covariates) ||
      any(!nzchar(covariates)) || anyDuplicated(covariates) ||
      any(!covariates %in% names(data))) {
    stop("covariates must contain unique, existing column names.", call. = FALSE)
  }
  covariates
}

.rerand_assignment_vector <- function(assignment, data) {
  if (!inherits(assignment, "rerand_assignment")) {
    stop("assignment must be a rerand_assignment object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  design <- assignment$design
  if (!inherits(design, "rerand_design")) {
    stop("assignment does not contain a valid rerand_design object.", call. = FALSE)
  }
  if (!is.null(design$id_name)) {
    id_name <- design$id_name
    if (!id_name %in% names(data)) {
      stop(sprintf("data must contain the design ID column '%s'.", id_name),
           call. = FALSE)
    }
    ids <- data[[id_name]]
    if (anyNA(ids) || anyDuplicated(ids)) {
      stop("The design ID column in data must be unique and non-missing.",
           call. = FALSE)
    }
    index <- match(ids, design$unit_id)
    if (anyNA(index) || length(index) != length(design$unit_id) ||
        anyDuplicated(index)) {
      stop("data IDs must match the design IDs exactly.", call. = FALSE)
    }
    return(assignment$Z[index])
  }
  if (nrow(data) != length(assignment$Z)) {
    stop("data must have the same number of rows as assignment.", call. = FALSE)
  }
  assignment$Z
}

.rerand_design_covariates <- function(assignment, data) {
  design <- assignment$design
  if (!is.null(design$id_name)) {
    ids <- data[[design$id_name]]
    index <- match(ids, design$unit_id)
    return(design$X[index, , drop = FALSE])
  }
  design$X
}

.rerand_prepare_estimation_covariates <- function(data, covariate_formula) {
  if (is.null(covariate_formula)) {
    return(NULL)
  }
  .rerand_prepare_model_matrix(
    data_frame = data,
    formula = covariate_formula,
    formula_missing = FALSE
  )$matrix
}

.rerand_parse_treatment <- function(treatment, data, treated = NULL) {
  if (is.numeric(treatment)) {
    return(.rerand_validate_assignment(treatment, nrow(data), min_group_size = 2L))
  }
  if (!(is.factor(treatment) || is.character(treatment)) || anyNA(treatment)) {
    stop("The treatment column must be binary numeric, factor, or character.",
         call. = FALSE)
  }
  levels <- unique(as.character(treatment))
  if (length(levels) != 2L) {
    stop("The treatment column must have exactly two observed levels.", call. = FALSE)
  }
  if (is.null(treated) || length(treated) != 1L || is.na(treated) ||
      !as.character(treated) %in% levels) {
    stop("treated must identify one observed treatment level for nonnumeric treatment.",
         call. = FALSE)
  }
  .rerand_validate_assignment(
    as.numeric(as.character(treatment) == as.character(treated)),
    nrow(data), min_group_size = 2L
  )
}

.rerand_parse_formula <- function(formula, data, treated = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must be a two-sided formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  response <- formula[[2L]]
  if (!is.name(response) || !as.character(response) %in% names(data)) {
    stop("The formula response must be one column in data.", call. = FALSE)
  }
  rhs <- .rerand_unwrap_parentheses(formula[[3L]])
  estimator <- NULL
  treatment_name <- NULL
  covariate_terms <- list()
  if (is.name(rhs)) {
    treatment_name <- as.character(rhs)
    estimator <- "dim"
  } else if (is.call(rhs) && identical(rhs[[1L]], as.name("+"))) {
    terms <- .rerand_flatten_addition(rhs)
    if (length(terms) < 1L || !is.name(terms[[1L]])) {
      stop("The first right-hand-side term must be an untransformed treatment column.",
           call. = FALSE)
    }
    treatment_name <- as.character(terms[[1L]])
    covariate_terms <- terms[-1L]
    estimator <- "ancova"
  } else if (is.call(rhs) && identical(rhs[[1L]], as.name("*"))) {
    left <- .rerand_unwrap_parentheses(rhs[[2L]])
    if (!is.name(left) || !is.call(rhs[[3L]]) ||
        !identical(rhs[[3L]][[1L]], as.name("("))) {
      stop("Lin formulas must have the form Y ~ Z * (x1 + x2).", call. = FALSE)
    }
    right <- .rerand_unwrap_parentheses(rhs[[3L]])
    treatment_name <- as.character(left)
    covariate_terms <- .rerand_flatten_addition(right)
    if (length(covariate_terms) == 0L) {
      stop("Lin formulas must include at least one covariate.", call. = FALSE)
    }
    estimator <- "lin"
  } else {
    stop("formula must be Y ~ Z, Y ~ Z + x, or Y ~ Z * (x1 + x2).",
         call. = FALSE)
  }
  if (!treatment_name %in% names(data)) {
    stop("The treatment column is not found in data.", call. = FALSE)
  }
  if (length(covariate_terms) > 0L && any(vapply(covariate_terms, function(term) {
    treatment_name %in% all.vars(term)
  }, logical(1)))) {
    stop("Treatment terms may not appear among the adjustment covariates.",
         call. = FALSE)
  }
  parsed_treatment <- .rerand_parse_treatment(data[[treatment_name]], data, treated)
  Y_obs <- data[[as.character(response)]]
  .rerand_assert_finite_numeric(Y_obs, "The outcome column")
  list(
    Y_obs = as.numeric(Y_obs),
    Z = parsed_treatment$Z,
    outcome_name = as.character(response),
    treatment_name = treatment_name,
    covariate_formula = .rerand_make_covariate_formula(
      covariate_terms, environment(formula)
    ),
    has_covariates = length(covariate_terms) > 0L,
    estimator = estimator
  )
}

.rerand_parse_selectors <- function(data, outcome, treatment, covariates,
                                    estimator, treated = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  if (length(estimator) != 1L || !is.character(estimator) ||
      is.na(estimator) || !estimator %in% c("dim", "ancova", "lin")) {
    stop("estimator must be one of 'dim', 'ancova', or 'lin'.", call. = FALSE)
  }
  outcome <- .rerand_validate_column_selector(outcome, "outcome", data)
  treatment <- .rerand_validate_column_selector(treatment, "treatment", data)
  covariates <- .rerand_validate_covariate_selectors(covariates, data)
  if (outcome == treatment || treatment %in% covariates || outcome %in% covariates) {
    stop("outcome, treatment, and covariates must be distinct columns.",
         call. = FALSE)
  }
  if (estimator %in% c("ancova", "lin") && length(covariates) == 0L) {
    stop(sprintf("estimator = '%s' requires at least one covariate.", estimator),
         call. = FALSE)
  }
  parsed_treatment <- .rerand_parse_treatment(data[[treatment]], data, treated)
  Y_obs <- data[[outcome]]
  .rerand_assert_finite_numeric(Y_obs, "The outcome column")
  covariate_formula <- if (length(covariates) == 0L) {
    NULL
  } else {
    stats::reformulate(covariates)
  }
  list(
    Y_obs = as.numeric(Y_obs),
    Z = parsed_treatment$Z,
    outcome_name = outcome,
    treatment_name = treatment,
    covariate_formula = covariate_formula,
    has_covariates = length(covariates) > 0L,
    estimator = estimator
  )
}

# Internal estimator helpers.

.calc_sample_stats <- function(Y_obs, Z, X, criterion) {
  assignment <- .rerand_validate_assignment(Z, length(Y_obs), min_group_size = 2L)
  Y1 <- Y_obs[assignment$Z == 1]
  Y0 <- Y_obs[assignment$Z == 0]
  n <- length(Y_obs)
  n1 <- assignment$n_treat
  n0 <- assignment$n_control
  r1 <- n1 / n
  r0 <- n0 / n
  tau_dim <- mean(Y1) - mean(Y0)
  V_tt_hat_1 <- stats::var(Y1) / r1 + stats::var(Y0) / r0

  base <- list(
    n = n,
    n1 = n1,
    n0 = n0,
    r1 = r1,
    r0 = r0,
    tau_dim = as.numeric(tau_dim),
    V_tt_hat_1 = .rerand_clamp_variance(V_tt_hat_1),
    V_tt_hat_2 = NULL,
    R2_hat = NULL,
    v_K_a = NULL,
    correction_factor = 1,
    K = criterion$K,
    criterion_type = criterion$type,
    accept_prob = criterion$accept_prob,
    acceptance_mass = criterion$acceptance_mass,
    threshold = criterion$threshold
  )

  if (is.null(X)) {
    return(base)
  }

  covariance <- .rerand_covariance(X)
  X1 <- X[assignment$Z == 1, , drop = FALSE]
  X0 <- X[assignment$Z == 0, , drop = FALSE]
  S_inv <- covariance$inverse

  S_Y1X <- matrix(stats::cov(Y1, X1), nrow = 1L)
  S_Y0X <- matrix(stats::cov(Y0, X0), nrow = 1L)
  S_tauX <- as.numeric((S_Y1X - S_Y0X) %*% S_inv %*% t(S_Y1X - S_Y0X))
  V_tt_hat_2 <- .rerand_clamp_variance(base$V_tt_hat_1 - S_tauX)

  S_Y1_given_X <- as.numeric(S_Y1X %*% S_inv %*% t(S_Y1X))
  S_Y0_given_X <- as.numeric(S_Y0X %*% S_inv %*% t(S_Y0X))
  numerator <- S_Y0_given_X / r0 + S_Y1_given_X / r1 - S_tauX
  R2_hat <- if (V_tt_hat_2 == 0) {
    if (abs(numerator) < 1e-10) 0 else stop("R2_hat is undefined.", call. = FALSE)
  } else {
    numerator / V_tt_hat_2
  }
  if (!is.finite(R2_hat) || R2_hat < -1e-8 || R2_hat > 1 + 1e-8) {
    stop("The sample R2 estimate is outside its valid range.", call. = FALSE)
  }
  R2_hat <- min(1, max(0, R2_hat))
  correction_factor <- .rerand_correction_factor(R2_hat, criterion)

  base$V_tt_hat_2 <- V_tt_hat_2
  base$R2_hat <- as.numeric(R2_hat)
  base$v_K_a <- criterion$v_K_a
  base$correction_factor <- as.numeric(correction_factor)
  base
}

.estimate_dim <- function(sample_stats) {
  se_neyman <- sqrt(sample_stats$V_tt_hat_1 * sample_stats$correction_factor) /
    sqrt(sample_stats$n)
  se_ding <- if (is.null(sample_stats$V_tt_hat_2)) {
    NULL
  } else {
    sqrt(sample_stats$V_tt_hat_2 * sample_stats$correction_factor) /
      sqrt(sample_stats$n)
  }
  list(
    tau_hat = sample_stats$tau_dim,
    se_neyman = as.numeric(se_neyman),
    se_ding = if (is.null(se_ding)) NULL else as.numeric(se_ding)
  )
}

.estimate_ancova <- function(Y_obs, Z, X) {
  if (is.null(X)) {
    stop("X is required for method 'ancova'.", call. = FALSE)
  }
  X_names <- paste0("X", seq_len(ncol(X)))
  model_data <- data.frame(Y = as.numeric(Y_obs), Z = as.numeric(Z), X)
  names(model_data) <- c("Y", "Z", X_names)
  formula <- stats::as.formula(paste(
    "Y ~ Z +", paste(X_names, collapse = " + ")
  ))
  fit <- stats::lm(formula, data = model_data)
  coefficient <- stats::coef(fit)["Z"]
  if (is.na(coefficient)) {
    stop("The ANCOVA model could not identify the treatment coefficient.",
         call. = FALSE)
  }
  covariance <- sandwich::vcovHC(fit, type = "HC2")
  se_ehw <- sqrt(covariance["Z", "Z"])
  list(tau_hat = coefficient, se_ehw = as.numeric(se_ehw), fit = fit)
}

.estimate_lin <- function(Y_obs, Z, X) {
  if (is.null(X)) {
    stop("X is required for method 'lin'.", call. = FALSE)
  }
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_names <- paste0("X", seq_len(ncol(X_centered)))
  model_data <- data.frame(Y = as.numeric(Y_obs), Z = as.numeric(Z), X_centered)
  names(model_data) <- c("Y", "Z", X_names)
  formula <- stats::as.formula(paste(
    "Y ~ Z +",
    paste(X_names, collapse = " + "),
    "+",
    paste(paste0("Z:", X_names), collapse = " + ")
  ))
  fit <- stats::lm(formula, data = model_data)
  coefficient <- stats::coef(fit)["Z"]
  if (is.na(coefficient)) {
    stop("The Lin model could not identify the treatment coefficient.", call. = FALSE)
  }
  covariance <- sandwich::vcovHC(fit, type = "HC2")
  se_ehw <- sqrt(covariance["Z", "Z"])
  list(tau_hat = coefficient, se_ehw = as.numeric(se_ehw), fit = fit)
}
