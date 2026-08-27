# Internal inference collection and comparison helpers.

.rerand_validate_inference_collection <- function(collection) {
  if (!inherits(collection, "rerand_inference_collection")) {
    stop("inferences must be created with c.rerand_inference().", call. = FALSE)
  }
  inferences <- collection$inferences
  if (!is.list(inferences) || length(inferences) < 2L ||
      any(!vapply(inferences, inherits, logical(1), "rerand_inference"))) {
    stop("at least two rerand_inference objects are required.", call. = FALSE)
  }
  object_names <- names(inferences)
  if (is.null(object_names) || anyNA(object_names) || any(!nzchar(object_names)) ||
      anyDuplicated(object_names)) {
    stop("inferences must be named with unique, non-empty method names.",
         call. = FALSE)
  }
  inferences
}

.rerand_inference_provenance <- function(inference) {
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

.rerand_assert_provenance_compatible <- function(reference, candidate) {
  fields <- c(
    "n", "Z", "unit_id", "design_method", "criterion_type",
    "accept_prob", "threshold", "n_treat"
  )
  for (field in fields) {
    if (!isTRUE(all.equal(reference[[field]], candidate[[field]]))) {
      stop("inference objects must share assignment and design provenance.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
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
  if (length(dots) == 0L || any(!vapply(dots, inherits, logical(1),
                                        "rerand_inference"))) {
    stop("all objects must be rerand_inference objects.", call. = FALSE)
  }
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
  inferences <- .rerand_validate_inference_collection(inferences)
  reference <- .rerand_inference_provenance(inferences[[1L]])
  levels <- vapply(inferences, function(x) x$level, numeric(1))
  if (any(abs(levels - levels[[1L]]) > 1e-12)) {
    stop("inference objects must use the same confidence level.", call. = FALSE)
  }
  for (inference in inferences[-1L]) {
    .rerand_assert_provenance_compatible(
      reference, .rerand_inference_provenance(inference)
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
