#' Generate a realized treatment assignment
#'
#' @description
#' Generates an assignment with exactly `n_treat` treated units and accepts it
#' when its Mahalanobis distance is below a prespecified criterion.
#'
#' @param design A `rerand_design` object created by [rerand_design()].
#' @param n_draws Number of assignments to generate.
#' @param max_tries Maximum number of assignments to draw. By default it is
#'   determined from `failure_prob` and the implied acceptance probability.
#' @param failure_prob Upper bound on the probability of exhausting the default
#'   number of draws without an accepted assignment.
#' @param seed Optional integer seed.
#' @param engine Either `"cpp"` or `"R"`.
#' @param on_failure Whether to error or warn after exhausting `max_tries`.
#' @return An object of class `rerand_assignment`.
#' @export
rerand_assign <- function(design, n_draws = 1L, max_tries = NULL,
                          failure_prob = 1e-6, seed = NULL,
                          engine = c("cpp", "R"),
                          on_failure = c("error", "warn")) {
  if (!inherits(design, "rerand_design")) {
    stop("design must be a rerand_design object.", call. = FALSE)
  }
  n_draws <- .rerand_validate_n_draws(n_draws)
  max_tries <- .rerand_resolve_max_tries(
    max_tries, design$criterion$acceptance_mass, failure_prob
  )
  seed <- .rerand_validate_seed(seed)
  engine <- .rerand_validate_engine(match.arg(engine))
  on_failure <- match.arg(on_failure)

  draw_one <- function() {
    if (engine == "R") {
      .design_r(design$whitened, design$n_treat, design$criterion$threshold,
                max_tries)
    } else {
      design_cpp_core(design$whitened, design$n_treat,
                      design$criterion$threshold, max_tries)
    }
  }
  draws <- .rerand_with_seed(seed, lapply(seq_len(n_draws), function(...) draw_one()))
  accepted <- vapply(draws, function(draw) isTRUE(draw$accepted), logical(1))
  if (any(!accepted)) {
    message <- sprintf(
      "Maximum tries (%d) exceeded for %d of %d assignment(s).",
      max_tries, sum(!accepted), n_draws
    )
    if (on_failure == "error") {
      stop(message, call. = FALSE)
    }
    warning(message, call. = FALSE)
  }

  assignments <- do.call(cbind, lapply(draws, `[[`, "Z"))
  colnames(assignments) <- paste0("draw", seq_len(n_draws))
  diagnostics <- data.frame(
    draw = seq_len(n_draws),
    accepted = accepted,
    tries = vapply(draws, `[[`, integer(1), "tries"),
    mahalanobis = vapply(draws, function(draw) as.numeric(draw$M), numeric(1))
  )
  result <- list(
    Z = as.numeric(assignments[, 1L]),
    pool = list(assignments = assignments, diagnostics = diagnostics),
    design = design,
    tries = diagnostics$tries[1L],
    mahalanobis = diagnostics$mahalanobis[1L],
    threshold = design$criterion$threshold,
    accept_prob = design$criterion$accept_prob,
    criterion_type = design$criterion$type,
    design_method = design$design_method,
    accepted = diagnostics$accepted[1L],
    engine = engine,
    n_treat = design$n_treat,
    n_control = design$n_control
  )
  class(result) <- c("rerand_assignment", "list")
  result
}

#' @export
print.rerand_assignment <- function(x, ...) {
  cat("Rerandomization assignment\n")
  cat("  engine:        ", x$engine, "\n", sep = "")
  cat("  criterion:     ", x$criterion_type, "\n", sep = "")
  cat("  design method: ", x$design_method, "\n", sep = "")
  cat("  accepted:      ", x$accepted, "\n", sep = "")
  cat("  tries:         ", x$tries, "\n", sep = "")
  cat("  Mahalanobis M: ", format(x$mahalanobis), "\n", sep = "")
  cat("  threshold:     ", format(x$threshold), "\n", sep = "")
  cat("  treated:       ", x$n_treat, "\n", sep = "")
  cat("  control:       ", x$n_control, "\n", sep = "")
  cat("  draws:         ", ncol(x$pool$assignments), "\n", sep = "")
  invisible(x)
}

#' @export
summary.rerand_assignment <- function(object, ...) {
  summary <- object[c(
    "accepted", "tries", "mahalanobis", "threshold", "accept_prob",
    "criterion_type", "engine", "n_treat", "n_control"
  )]
  summary$design_method <- object$design_method
  class(summary) <- c("summary.rerand_assignment", "list")
  summary
}

#' @export
print.summary.rerand_assignment <- function(x, ...) {
  print.rerand_assignment(x, ...)
  invisible(x)
}

#' Return design data with a treatment-assignment column
#'
#' @param x A `rerand_assignment` object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param treatment Name of the treatment-assignment column.
#' @param overwrite Whether to replace an existing column named `treatment`.
#' @param ... Ignored.
#' @export
as.data.frame.rerand_assignment <- function(x, row.names = NULL,
                                            optional = FALSE,
                                            treatment = "Z",
                                            overwrite = FALSE, ...) {
  if (length(treatment) != 1L || !is.character(treatment) || is.na(treatment) ||
      !nzchar(treatment)) {
    stop("treatment must be one non-empty column name.", call. = FALSE)
  }
  data <- as.data.frame(x$design$data)
  if (treatment %in% names(data) && !isTRUE(overwrite)) {
    stop("treatment already exists in the design data; use overwrite = TRUE.",
         call. = FALSE)
  }
  data[[treatment]] <- x$Z
  data
}
