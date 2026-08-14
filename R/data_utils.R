# Internal data and formula preparation helpers.

.rerand_supported_column <- function(x) {
  is.numeric(x) || is.logical(x) || is.factor(x) || is.character(x)
}

.rerand_validate_id <- function(data, id) {
  n <- nrow(data)
  if (is.null(id)) {
    return(list(name = NULL, values = seq_len(n)))
  }
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

.rerand_validate_one_sided_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("formula must be a one-sided formula.", call. = FALSE)
  }
  formula
}

.rerand_prepare_model_matrix <- function(data_frame, formula, id_name = NULL,
                                         formula_missing = FALSE) {
  if (formula_missing) {
    selected <- if (is.null(id_name)) names(data_frame) else {
      setdiff(names(data_frame), id_name)
    }
    if (length(selected) == 0L) {
      stop("data must contain at least one covariate column.", call. = FALSE)
    }
    bad <- !vapply(data_frame[selected], .rerand_supported_column, logical(1))
    if (any(bad)) {
      stop(
        sprintf(
          "Unsupported covariate column type: %s.",
          paste(selected[bad], collapse = ", ")
        ),
        call. = FALSE
      )
    }
    formula <- stats::reformulate(selected)
    environment(formula) <- parent.frame()
  } else {
    formula <- .rerand_validate_one_sided_formula(formula)
    variables <- all.vars(formula)
    missing_variables <- setdiff(variables, names(data_frame))
    if (length(missing_variables) > 0L) {
      stop(
        sprintf(
          "Variables not found in data: %s.",
          paste(missing_variables, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(id_name) && id_name %in% variables) {
      stop("The id column cannot be used as a covariate.", call. = FALSE)
    }
    bad <- !vapply(data_frame[variables], .rerand_supported_column, logical(1))
    if (any(bad)) {
      stop(
        sprintf(
          "Unsupported covariate column type: %s.",
          paste(variables[bad], collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  model_frame <- tryCatch(
    stats::model.frame(formula, data = data_frame, na.action = stats::na.fail),
    error = function(error) {
      stop(conditionMessage(error), call. = FALSE)
    }
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
    stop("The covariate formula must produce at least one column.",
         call. = FALSE)
  }
  if (anyDuplicated(colnames(matrix))) {
    stop("The encoded covariate columns must have unique names.", call. = FALSE)
  }
  if (any(!is.finite(matrix))) {
    stop("The encoded covariate matrix must contain only finite values.",
         call. = FALSE)
  }

  xlevels <- lapply(model_frame[vapply(model_frame, is.factor, logical(1))], levels)
  list(
    matrix = matrix,
    model_formula = formula,
    terms = stats::terms(formula),
    xlevels = xlevels,
    contrasts = attr(matrix, "contrasts")
  )
}

.rerand_prepare_covariates <- function(data, formula = NULL, id = NULL) {
  if (is.matrix(data)) {
    data <- .rerand_validate_matrix(data, name = "data")
    if (!is.null(id)) {
      stop("id is only supported when data is a data frame.", call. = FALSE)
    }
    if (is.null(formula)) {
      X <- unname(data)
      colnames(X) <- colnames(data) %||% paste0("X", seq_len(ncol(data)))
      return(list(
        data = data,
        X = X,
        unit_id = seq_len(nrow(data)),
        id_name = NULL,
        formula = NULL,
        model_formula = NULL,
        terms = NULL,
        xlevels = list(),
        contrasts = NULL
      ))
    }
    if (is.null(colnames(data)) || any(!nzchar(colnames(data))) ||
        anyDuplicated(colnames(data))) {
      stop("Matrix data must have unique column names when formula is supplied.",
           call. = FALSE)
    }
    data_frame <- as.data.frame(data, optional = TRUE)
    prepared <- .rerand_prepare_model_matrix(
      data_frame = data_frame,
      formula = formula,
      formula_missing = FALSE
    )
    return(c(
      list(
        data = data,
        X = prepared$matrix,
        unit_id = seq_len(nrow(data)),
        id_name = NULL,
        formula = formula
      ),
      prepared[c("model_formula", "terms", "xlevels", "contrasts")]
    ))
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
  id_info <- .rerand_validate_id(data, id)
  prepared <- .rerand_prepare_model_matrix(
    data_frame = data,
    formula = formula,
    id_name = id_info$name,
    formula_missing = is.null(formula)
  )
  c(
    list(
      data = data,
      X = prepared$matrix,
      unit_id = id_info$values,
      id_name = id_info$name,
      formula = formula
    ),
    prepared[c("model_formula", "terms", "xlevels", "contrasts")]
  )
}

.rerand_whiten_covariates <- function(X, tol = 1e-10) {
  X <- .rerand_validate_matrix(X)
  tol <- .rerand_validate_tol(tol)
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
  list(
    centered = centered,
    whitened = whitened,
    whitening = whitening,
    singular_values = singular_values,
    effective_rank = as.integer(sum(keep)),
    tol = tol
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
