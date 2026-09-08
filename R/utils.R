# Shared non-validation helpers.

.or_else <- function(x, y) {
  if (is.null(x)) y else x
}

.with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.supported_column <- function(x) {
  is.numeric(x) || is.logical(x) || is.factor(x) || is.character(x)
}

.prepare_model_matrix <- function(data_frame, formula, id_name = NULL,
                                  formula_missing = FALSE) {
  if (formula_missing) {
    selected <- if (is.null(id_name)) names(data_frame) else {
      setdiff(names(data_frame), id_name)
    }
    if (length(selected) == 0L) {
      stop("data must contain at least one covariate column.", call. = FALSE)
    }
    bad <- !vapply(data_frame[selected], .supported_column, logical(1))
    if (any(bad)) {
      stop(sprintf("Unsupported covariate column type: %s.",
                   paste(selected[bad], collapse = ", ")), call. = FALSE)
    }
    formula <- stats::reformulate(selected)
    environment(formula) <- parent.frame()
  } else {
    formula <- .validate_one_sided_formula(formula)
    variables <- all.vars(formula)
    missing_variables <- setdiff(variables, names(data_frame))
    if (length(missing_variables) > 0L) {
      stop(sprintf("Variables not found in data: %s.",
                   paste(missing_variables, collapse = ", ")), call. = FALSE)
    }
    if (!is.null(id_name) && id_name %in% variables) {
      stop("The id column cannot be used as a covariate.", call. = FALSE)
    }
    bad <- !vapply(data_frame[variables], .supported_column, logical(1))
    if (any(bad)) {
      stop(sprintf("Unsupported covariate column type: %s.",
                   paste(variables[bad], collapse = ", ")), call. = FALSE)
    }
  }
  model_frame <- tryCatch(
    stats::model.frame(formula, data = data_frame, na.action = stats::na.fail),
    error = function(error) stop(conditionMessage(error), call. = FALSE)
  )
  matrix <- stats::model.matrix(formula, data = model_frame)
  assign <- attr(matrix, "assign")
  if (!is.null(assign)) {
    matrix <- matrix[, assign != 0L, drop = FALSE]
  } else if ("(Intercept)" %in% colnames(matrix)) {
    matrix <- matrix[, colnames(matrix) != "(Intercept)", drop = FALSE]
  }
  storage.mode(matrix) <- "double"
  if (ncol(matrix) < 1L) {
    stop("The covariate formula must produce at least one column.", call. = FALSE)
  }
  if (anyDuplicated(colnames(matrix))) {
    stop("The encoded covariate columns must have unique names.", call. = FALSE)
  }
  if (any(!is.finite(matrix))) {
    stop("The encoded covariate matrix must contain only finite values.",
         call. = FALSE)
  }
  xlevels <- lapply(model_frame[vapply(model_frame, is.factor, logical(1))], levels)
  list(matrix = matrix, model_formula = formula, terms = stats::terms(formula),
       xlevels = xlevels, contrasts = attr(matrix, "contrasts"))
}

.prepare_covariates <- function(data, formula = NULL, id = NULL) {
  if (is.matrix(data)) {
    data <- .validate_matrix(data, name = "data")
    if (!is.null(id)) {
      stop("id is only supported when data is a data frame.", call. = FALSE)
    }
    if (is.null(formula)) {
      X <- unname(data)
      colnames(X) <- .or_else(colnames(data), paste0("X", seq_len(ncol(data))))
      return(list(data = data, X = X, unit_id = seq_len(nrow(data)),
                  id_name = NULL, formula = NULL, model_formula = NULL,
                  terms = NULL, xlevels = list(), contrasts = NULL))
    }
    if (is.null(colnames(data)) || any(!nzchar(colnames(data))) ||
        anyDuplicated(colnames(data))) {
      stop("Matrix data must have unique column names when formula is supplied.",
           call. = FALSE)
    }
    prepared <- .prepare_model_matrix(as.data.frame(data, optional = TRUE),
                                      formula, formula_missing = FALSE)
    return(c(list(data = data, X = prepared$matrix, unit_id = seq_len(nrow(data)),
                  id_name = NULL, formula = formula),
             prepared[c("model_formula", "terms", "xlevels", "contrasts")]))
  }
  if (!is.data.frame(data)) {
    stop("data must be a numeric matrix or a data frame.", call. = FALSE)
  }
  if (nrow(data) < 2L || ncol(data) < 1L) {
    stop("data must have at least two rows and one column.", call. = FALSE)
  }
  if (any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("data must have unique, non-empty column names.", call. = FALSE)
  }
  id_info <- .validate_id(data, id)
  prepared <- .prepare_model_matrix(data, formula, id_info$name,
                                    formula_missing = is.null(formula))
  c(list(data = data, X = prepared$matrix, unit_id = id_info$values,
         id_name = id_info$name, formula = formula),
    prepared[c("model_formula", "terms", "xlevels", "contrasts")])
}

.whiten_covariates <- function(X, tol = 1e-10) {
  X <- .validate_matrix(X)
  tol <- .validate_tol(tol)
  centered <- scale(X, center = TRUE, scale = FALSE)
  decomposition <- svd(centered, nu = 0L, nv = min(dim(centered)))
  if (length(decomposition$d) == 0L || decomposition$d[1L] == 0) {
    stop("X must contain at least one non-constant covariate direction.",
         call. = FALSE)
  }
  keep <- decomposition$d > tol * decomposition$d[1L]
  if (!any(keep)) {
    stop("X must contain at least one non-constant covariate direction.",
         call. = FALSE)
  }
  singular_values <- decomposition$d[keep]
  loadings <- decomposition$v[, keep, drop = FALSE]
  whitening <- sweep(loadings, 2L, singular_values, "/") * sqrt(nrow(X) - 1)
  whitened <- centered %*% whitening
  colnames(whitened) <- paste0("W", seq_len(ncol(whitened)))
  list(centered = centered, whitened = whitened, whitening = whitening,
       singular_values = singular_values, effective_rank = as.integer(sum(keep)),
       tol = tol)
}

.resolve_criterion <- function(accept_prob = NULL, threshold = NULL,
                               K, require_criterion = TRUE) {
  .validate_criterion_inputs(accept_prob, threshold, K, require_criterion)
  K <- as.integer(K)
  if (is.null(accept_prob) && is.null(threshold)) {
    if (require_criterion) {
      stop("Exactly one of accept_prob or threshold must be supplied.",
           call. = FALSE)
    }
    accept_prob <- 1
  }
  if (!is.null(accept_prob) && !is.null(threshold)) {
    stop("accept_prob and threshold are mutually exclusive.", call. = FALSE)
  }
  if (!is.null(accept_prob)) {
    threshold <- stats::qchisq(accept_prob, df = K)
    acceptance_mass <- accept_prob
    type <- "probability"
  } else {
    acceptance_mass <- stats::pchisq(threshold, df = K)
    if (!is.finite(acceptance_mass) || acceptance_mass <= 0) {
      stop("threshold does not define a positive acceptance region.",
           call. = FALSE)
    }
    type <- "threshold"
  }
  v_K_a <- if (acceptance_mass == 1) 1 else {
    stats::pchisq(threshold, df = K + 2) / acceptance_mass
  }
  list(type = type,
       accept_prob = if (type == "probability") as.numeric(accept_prob) else NULL,
       threshold = as.numeric(threshold), K = K,
       acceptance_mass = as.numeric(acceptance_mass), v_K_a = as.numeric(v_K_a))
}
