# Internal inference collection and comparison helpers.

.inference_provenance <- function(inference) {
  estimate <- inference$estimate
  assignment <- estimate$assignment
  design <- assignment$design
  list(
    n = estimate$n,
    Z = as.numeric(assignment$Z),
    unit_id = design$unit_id,
    design_method = estimate$design_method,
    criterion_type = estimate$criterion_type,
    accept_prob = estimate$accept_prob,
    threshold = estimate$threshold,
    n_treat = design$n_treat
  )
}

#' Combine inference objects for comparison
#'
#' @param ... Named `rerand_inference` objects.
#' @param recursive Ignored; included for compatibility with the base generic.
#' @return A named inference collection.
#' @export
c.rerand_inference <- function(..., recursive = FALSE) {
  dots <- list(...)
  if (length(dots) == 1L && inherits(dots[[1L]], "rerand_inference_collection")) {
    return(dots[[1L]])
  }
  .validate_inference_objects(dots)
  object_names <- names(dots)
  if (is.null(object_names)) {
    object_names <- rep("", length(dots))
  }
  result <- list(inferences = dots)
  names(result$inferences) <- object_names
  class(result) <- c("rerand_inference_collection", "list")
  result
}

#' Compare inference results from different estimators
#'
#' @param inferences A named collection created with [c.rerand_inference()].
#' @return An object of class `rerand_compare`.
#' @export
rerand_compare <- function(inferences) {
  inferences <- .validate_inference_collection(inferences)
  reference <- .inference_provenance(inferences[[1L]])
  levels <- vapply(inferences, function(x) x$level, numeric(1))
  .validate_confidence_levels(levels)
  for (inference in inferences[-1L]) {
    .validate_provenance_compatible(
      reference, .inference_provenance(inference)
    )
  }
  object_names <- names(inferences)
  table <- data.frame(
    method = object_names,
    estimator = vapply(inferences, `[[`, character(1), "estimator"),
    estimate = vapply(inferences, `[[`, numeric(1), "estimate_value"),
    standard_error = vapply(inferences, `[[`, numeric(1), "standard_error"),
    lower = vapply(inferences, function(x) x$interval[1L, 1L], numeric(1)),
    upper = vapply(inferences, function(x) x$interval[1L, 2L], numeric(1)),
    design_method = vapply(inferences, `[[`, character(1), "design_method"),
    reference_distribution = vapply(
      inferences, `[[`, character(1), "reference_distribution"
    ),
    level = levels,
    row.names = object_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  result <- list(
    results = table,
    inferences = inferences,
    level = levels[[1L]],
    design_method = reference$design_method,
    criterion_type = reference$criterion_type,
    accept_prob = reference$accept_prob,
    threshold = reference$threshold,
    n = reference$n
  )
  class(result) <- c("rerand_compare", "list")
  result
}

#' @export
print.rerand_compare <- function(x, ...) {
  cat("Rerandomization inference comparison\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
summary.rerand_compare <- function(object, ...) {
  result <- object[c(
    "results", "level", "design_method", "criterion_type", "accept_prob",
    "threshold", "n"
  )]
  class(result) <- c("summary.rerand_compare", "list")
  result
}

#' @export
print.summary.rerand_compare <- function(x, ...) {
  print.rerand_compare(x, ...)
  invisible(x)
}
