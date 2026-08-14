# Internal helpers for formula-based estimation.

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

.rerand_make_covariate_formula <- function(terms, environment) {
  if (length(terms) == 0L) {
    return(NULL)
  }
  rhs <- Reduce(function(left, right) call("+", left, right), terms)
  stats::as.formula(call("~", rhs), env = environment)
}

.rerand_prepare_estimate_formula <- function(formula, data, treated = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("formula must be a two-sided formula such as Y ~ Z + x.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }
  response <- formula[[2L]]
  if (!is.name(response) || !as.character(response) %in% names(data)) {
    stop("The formula response must be one column in data.", call. = FALSE)
  }
  rhs_terms <- .rerand_flatten_addition(formula[[3L]])
  if (length(rhs_terms) == 0L || !is.name(rhs_terms[[1L]])) {
    stop("The first right-hand-side term must be an untransformed treatment column.",
         call. = FALSE)
  }
  treatment_name <- as.character(rhs_terms[[1L]])
  if (!treatment_name %in% names(data)) {
    stop("The treatment column is not found in data.", call. = FALSE)
  }
  covariate_terms <- rhs_terms[-1L]
  if (any(vapply(covariate_terms, function(term) {
    treatment_name %in% all.vars(term)
  }, logical(1)))) {
    stop("Treatment interactions and repeated treatment terms are not supported.",
         call. = FALSE)
  }
  Y_obs <- data[[as.character(response)]]
  .rerand_assert_finite_numeric(Y_obs, "The outcome column")
  treatment <- data[[treatment_name]]
  if (is.numeric(treatment)) {
    assignment <- .rerand_validate_assignment(treatment, nrow(data),
                                               min_group_size = 2L)
  } else {
    if (!(is.factor(treatment) || is.character(treatment)) || anyNA(treatment)) {
      stop("The treatment column must be binary numeric, factor, or character.",
           call. = FALSE)
    }
    levels <- unique(as.character(treatment))
    if (length(levels) != 2L) {
      stop("The treatment column must have exactly two observed levels.",
           call. = FALSE)
    }
    if (is.null(treated) || length(treated) != 1L || is.na(treated) ||
        !as.character(treated) %in% levels) {
      stop("treated must identify one observed treatment level for nonnumeric treatment.",
           call. = FALSE)
    }
    assignment <- .rerand_validate_assignment(
      as.numeric(as.character(treatment) == as.character(treated)), nrow(data),
      min_group_size = 2L
    )
  }
  list(
    Y_obs = as.numeric(Y_obs),
    Z = assignment$Z,
    treatment_name = treatment_name,
    covariate_formula = .rerand_make_covariate_formula(
      covariate_terms, environment(formula)
    ),
    has_covariates = length(covariate_terms) > 0L
  )
}

.rerand_assignment_from_design <- function(design, data) {
  if (!inherits(design, "rerand_design_result")) {
    stop("design must be a rerand_design_result object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data frame when design is supplied.", call. = FALSE)
  }
  spec <- design$spec
  if (!is.null(spec$id_name)) {
    id_name <- spec$id_name
    if (!id_name %in% names(data)) {
      stop(sprintf("data must contain the design ID column '%s'.", id_name),
           call. = FALSE)
    }
    ids <- data[[id_name]]
    if (anyNA(ids) || anyDuplicated(ids)) {
      stop("The design ID column in data must be unique and non-missing.",
           call. = FALSE)
    }
    index <- match(ids, spec$unit_id)
    if (anyNA(index) || length(index) != length(spec$unit_id) ||
        anyDuplicated(index)) {
      stop("data IDs must match the design IDs exactly.", call. = FALSE)
    }
    return(design$Z[index])
  }
  if (nrow(data) != length(design$Z)) {
    stop("data must have the same number of rows as design.", call. = FALSE)
  }
  design$Z
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
