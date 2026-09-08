# Public-input validation helpers.

.validate_finite_numeric <- function(x, name, length = NULL) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  if (!is.null(length) && length(x) != length) {
    stop(sprintf("%s must have length %d.", name, length), call. = FALSE)
  }
  invisible(x)
}

.validate_matrix <- function(X, name = "X", n = NULL) {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop(sprintf("%s must be a finite numeric matrix.", name), call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L) {
    stop(sprintf("%s must have at least two rows and one column.", name),
         call. = FALSE)
  }
  if (!is.null(n) && nrow(X) != n) {
    stop(sprintf("%s must have %d rows.", name, n), call. = FALSE)
  }
  X
}

.validate_n_treat <- function(n_treat, n) {
  if (length(n_treat) != 1L || is.na(n_treat) ||
      !is.numeric(n_treat) || !is.finite(n_treat) ||
      n_treat != as.integer(n_treat)) {
    stop("n_treat must be one integer.", call. = FALSE)
  }
  n_treat <- as.integer(n_treat)
  if (n_treat < 1L || n_treat >= n) {
    stop("n_treat must be between 1 and the number of units minus 1.",
         call. = FALSE)
  }
  n_treat
}

.validate_tol <- function(tol) {
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) ||
      tol <= 0 || tol >= 1) {
    stop("tol must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(tol)
}

.validate_engine <- function(engine) {
  if (length(engine) != 1L || is.na(engine) ||
      !engine %in% c("cpp", "R")) {
    stop("engine must be one of 'cpp' or 'R'.", call. = FALSE)
  }
  engine
}

.validate_assignment <- function(Z, n, min_group_size = 1L) {
  .validate_finite_numeric(Z, "Z", n)
  if (any(!Z %in% c(0, 1))) {
    stop("Z must contain only 0 and 1.", call. = FALSE)
  }
  n_treat <- sum(Z == 1)
  n_control <- sum(Z == 0)
  if (n_treat < min_group_size || n_control < min_group_size) {
    stop("Z must contain enough observations in both treatment groups.",
         call. = FALSE)
  }
  list(Z = as.numeric(Z), n_treat = n_treat, n_control = n_control)
}

.validate_simulation_inputs <- function(R2, K, alpha, n_sim) {
  if (length(R2) != 1L || !is.numeric(R2) || !is.finite(R2) ||
      R2 < 0 || R2 > 1) {
    stop("R2 must be a finite number between 0 and 1.", call. = FALSE)
  }
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) || K < 1 ||
      K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha <= 0.5 || alpha >= 1) {
    stop("alpha must be strictly between 0.5 and 1.", call. = FALSE)
  }
  if (length(n_sim) != 1L || !is.numeric(n_sim) || !is.finite(n_sim) ||
      n_sim < 1 || n_sim != as.integer(n_sim)) {
    stop("n_sim must be a positive integer.", call. = FALSE)
  }
  list(R2 = as.numeric(R2), K = as.integer(K), alpha = as.numeric(alpha),
       n_sim = as.integer(n_sim))
}

.validate_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite integer or NULL.", call. = FALSE)
  }
  as.integer(seed)
}

.validate_n_draws <- function(n_draws) {
  if (length(n_draws) != 1L || !is.numeric(n_draws) ||
      !is.finite(n_draws) || n_draws < 1 ||
      n_draws != as.integer(n_draws)) {
    stop("n_draws must be one positive integer.", call. = FALSE)
  }
  as.integer(n_draws)
}

.validate_id <- function(data, id) {
  n <- nrow(data)
  if (is.null(id)) return(list(name = NULL, values = seq_len(n)))
  if (!is.data.frame(data)) {
    stop("id can only name a column when data is a data frame.", call. = FALSE)
  }
  if (length(id) != 1L || !is.character(id) || is.na(id) || !nzchar(id) ||
      !id %in% names(data)) {
    stop("id must name one column in data.", call. = FALSE)
  }
  values <- data[[id]]
  if (is.list(values) || is.matrix(values) || anyNA(values)) {
    stop("The id column must be an atomic vector without missing values.",
         call. = FALSE)
  }
  if (anyDuplicated(values)) {
    stop("The id column must contain unique values.", call. = FALSE)
  }
  list(name = id, values = values)
}

.validate_one_sided_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("formula must be a one-sided formula.", call. = FALSE)
  }
  formula
}

.validate_column_selector <- function(selector, name, data) {
  if (length(selector) != 1L || !is.character(selector) || is.na(selector) ||
      !nzchar(selector) || !selector %in% names(data)) {
    stop(sprintf("%s must name one column in data.", name), call. = FALSE)
  }
  selector
}

.validate_covariate_selectors <- function(covariates, data) {
  if (is.null(covariates)) return(character())
  if (!is.character(covariates) || anyNA(covariates) ||
      any(!nzchar(covariates)) || anyDuplicated(covariates) ||
      any(!covariates %in% names(data))) {
    stop("covariates must contain unique, existing column names.",
         call. = FALSE)
  }
  covariates
}

.validate_data_frame <- function(data, name = "data") {
  if (!is.data.frame(data)) {
    stop(sprintf("%s must be a data frame.", name), call. = FALSE)
  }
  invisible(data)
}

.validate_inference_inputs <- function(estimate, level, integration_tol) {
  if (!inherits(estimate, "rerand_estimate")) {
    stop("estimate must be a rerand_estimate object.", call. = FALSE)
  }
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(integration_tol) != 1L || !is.numeric(integration_tol) ||
      !is.finite(integration_tol) || integration_tol <= 0 ||
      integration_tol >= 1) {
    stop("integration_tol must be a finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  invisible(TRUE)
}
