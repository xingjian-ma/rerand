#' Evaluate inference against finite-population potential outcomes
#'
#' @description
#' Computes the finite-population average treatment effect from a separate
#' two-column data frame and reports estimation error and interval coverage.
#'
#' @param result A `rerand_inference` or `rerand_compare` object.
#' @param newdata Two-column data frame containing the potential outcomes.
#' @param potential_outcomes A named character vector mapping `control` and
#'   `treated` to the two columns in `newdata`.
#' @return An object of class `rerand_evaluation`.
#' @export
rerand_evaluate <- function(result, newdata, potential_outcomes) {
  inferences <- if (inherits(result, "rerand_inference")) {
    named <- list(result)
    names(named) <- result$estimator
    named
  } else if (inherits(result, "rerand_compare")) {
    result$inferences
  } else {
    stop("result must be a rerand_inference or rerand_compare object.",
         call. = FALSE)
  }
  if (!is.data.frame(newdata) || ncol(newdata) != 2L) {
    stop("newdata must be a data frame with exactly two potential-outcome columns.",
         call. = FALSE)
  }
  if (!is.character(potential_outcomes) || length(potential_outcomes) != 2L ||
      is.null(names(potential_outcomes)) ||
      !setequal(names(potential_outcomes), c("control", "treated")) ||
      anyNA(potential_outcomes) || any(!nzchar(potential_outcomes)) ||
      anyDuplicated(potential_outcomes) ||
      any(!potential_outcomes %in% names(newdata))) {
    stop("potential_outcomes must map control and treated to distinct newdata columns.",
         call. = FALSE)
  }
  if (any(!vapply(newdata, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(newdata)))) {
    stop("newdata potential outcomes must be finite numeric columns.",
         call. = FALSE)
  }
  if (length(inferences) < 1L || any(!vapply(
    inferences, inherits, logical(1), "rerand_inference"
  ))) {
    stop("result does not contain valid rerand_inference objects.", call. = FALSE)
  }
  if (any(nrow(newdata) != vapply(
    inferences, function(x) x$estimate$n, integer(1)
  ))) {
    stop("newdata must have one row for each design unit.", call. = FALSE)
  }

  control <- newdata[[potential_outcomes[["control"]]]]
  treated <- newdata[[potential_outcomes[["treated"]]]]
  true_effect <- mean(treated - control)
  object_names <- names(inferences)
  if (is.null(object_names) || any(!nzchar(object_names))) {
    object_names <- vapply(inferences, `[[`, character(1), "estimator")
  }
  table <- data.frame(
    method = object_names,
    estimator = vapply(inferences, `[[`, character(1), "estimator"),
    estimate = vapply(inferences, `[[`, numeric(1), "estimate_value"),
    signed_error = vapply(inferences, function(x) x$estimate_value - true_effect,
                          numeric(1)),
    absolute_error = vapply(
      inferences, function(x) abs(x$estimate_value - true_effect), numeric(1)
    ),
    squared_error = vapply(
      inferences, function(x) (x$estimate_value - true_effect)^2, numeric(1)
    ),
    lower = vapply(inferences, function(x) x$interval[1L, 1L], numeric(1)),
    upper = vapply(inferences, function(x) x$interval[1L, 2L], numeric(1)),
    covered = vapply(
      inferences,
      function(x) isTRUE(true_effect >= x$interval[1L, 1L] &&
                          true_effect <= x$interval[1L, 2L]),
      logical(1)
    ),
    row.names = object_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  evaluation <- list(
    results = table,
    true_effect = as.numeric(true_effect),
    potential_outcome_names = potential_outcomes,
    n = nrow(newdata),
    source = result
  )
  class(evaluation) <- c("rerand_evaluation", "list")
  evaluation
}

#' @export
print.rerand_evaluation <- function(x, ...) {
  cat("Rerandomization evaluation\n")
  cat("  true effect: ", format(x$true_effect), "\n", sep = "")
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
summary.rerand_evaluation <- function(object, ...) {
  result <- object[c(
    "results", "true_effect", "potential_outcome_names", "n"
  )]
  class(result) <- c("summary.rerand_evaluation", "list")
  result
}

#' @export
print.summary.rerand_evaluation <- function(x, ...) {
  print.rerand_evaluation(x, ...)
  invisible(x)
}
