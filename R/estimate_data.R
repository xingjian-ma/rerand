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
    right <- .rerand_unwrap_parentheses(rhs[[3L]])
    if (!is.name(left)) {
      stop("Lin formulas must have the form Y ~ Z * (x1 + x2).", call. = FALSE)
    }
    treatment_name <- as.character(left)
    covariate_terms <- .rerand_flatten_addition(right)
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
