#' Prepare a reusable rerandomization specification
#'
#' @description
#' Encodes covariates, determines their effective rank, constructs a whitened
#' covariate matrix, and resolves the rerandomization acceptance criterion.
#' The returned object can be reused by `rerand_draw()`.
#'
#' @param data Numeric matrix or data frame of pretreatment covariates.
#' @param n_treat Number of treated units.
#' @param formula Optional one-sided formula selecting and constructing
#'   covariates. With `NULL`, all columns except `id` are used.
#' @param accept_prob Target acceptance probability. Mutually exclusive with
#'   `threshold`.
#' @param threshold Direct Mahalanobis threshold. Mutually exclusive with
#'   `accept_prob`.
#' @param id Optional name of a unique ID column in data-frame input.
#' @param tol Relative tolerance used to determine the effective covariate rank.
#' @return An object of class `rerand_spec`.
#' @export
rerand_spec <- function(data, n_treat, formula = NULL,
                        accept_prob = NULL, threshold = NULL,
                        id = NULL, tol = 1e-10) {
  prepared <- .rerand_prepare_covariates(data, formula = formula, id = id)
  n_treat <- .rerand_validate_n_treat(n_treat, nrow(prepared$X))
  transformed <- .rerand_whiten_covariates(prepared$X, tol = tol)
  criterion <- .rerand_resolve_criterion(
    accept_prob = accept_prob,
    threshold = threshold,
    K = transformed$effective_rank,
    require_criterion = TRUE
  )
  result <- c(
    prepared,
    transformed,
    list(
      n = nrow(prepared$X),
      n_treat = n_treat,
      n_control = nrow(prepared$X) - n_treat,
      criterion = criterion
    )
  )
  class(result) <- c("rerand_spec", "list")
  result
}

#' @export
print.rerand_spec <- function(x, ...) {
  cat("Rerandomization specification\n")
  cat("  units:          ", x$n, "\n", sep = "")
  cat("  treated:        ", x$n_treat, "\n", sep = "")
  cat("  control:        ", x$n_control, "\n", sep = "")
  cat("  encoded columns:", ncol(x$X), "\n")
  cat("  effective rank: ", x$effective_rank, "\n", sep = "")
  cat("  criterion:      ", x$criterion$type, "\n", sep = "")
  cat("  threshold:      ", format(x$criterion$threshold), "\n", sep = "")
  invisible(x)
}

#' @export
summary.rerand_spec <- function(object, ...) {
  result <- list(
    n = object$n,
    n_treat = object$n_treat,
    n_control = object$n_control,
    n_covariates = ncol(object$X),
    effective_rank = object$effective_rank,
    criterion_type = object$criterion$type,
    accept_prob = object$criterion$accept_prob,
    implied_accept_prob = object$criterion$acceptance_mass,
    threshold = object$criterion$threshold,
    covariate_names = colnames(object$X)
  )
  class(result) <- c("summary.rerand_spec", "list")
  result
}

#' @export
print.summary.rerand_spec <- function(x, ...) {
  cat("Rerandomization specification summary\n")
  cat("  units:          ", x$n, "\n", sep = "")
  cat("  treated:        ", x$n_treat, "\n", sep = "")
  cat("  control:        ", x$n_control, "\n", sep = "")
  cat("  encoded columns:", x$n_covariates, "\n")
  cat("  effective rank: ", x$effective_rank, "\n", sep = "")
  cat("  criterion:      ", x$criterion_type, "\n", sep = "")
  invisible(x)
}
