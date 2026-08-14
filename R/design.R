#' Generate a rerandomized treatment assignment
#'
#' @description
#' Generates an assignment with exactly `n_treat` treated units and accepts it
#' when its Mahalanobis distance is below a prespecified criterion.
#'
#' @param data Numeric matrix or data frame of pretreatment covariates.
#' @param n_treat Number of treated units.
#' @param accept_prob Target acceptance probability used to derive the
#'   Mahalanobis threshold. Mutually exclusive with `threshold`.
#' @param threshold Direct Mahalanobis threshold. Mutually exclusive with
#'   `accept_prob`.
#' @param formula Optional one-sided formula selecting balance covariates.
#' @param id Optional name of a unique ID column in data-frame input.
#' @param max_tries Maximum number of assignments to draw. By default it is
#'   determined from `failure_prob` and the implied acceptance probability.
#' @param failure_prob Upper bound on the probability of exhausting the default
#'   number of draws without an accepted assignment.
#' @param seed Optional integer seed.
#' @param engine Either `"cpp"` or `"R"`.
#' @param on_failure Whether to error or warn after exhausting `max_tries`.
#' @param tol Relative tolerance used to determine the effective covariate rank.
#' @return An object of class `rerand_design_result`.
#' @export
rerand_design <- function(data, n_treat, formula = NULL, accept_prob = NULL,
                          threshold = NULL, id = NULL, max_tries = NULL,
                          failure_prob = 1e-6, seed = NULL,
                          engine = c("cpp", "R"),
                          on_failure = c("error", "warn"), tol = 1e-10) {
  spec <- rerand_spec(
    data = data, n_treat = n_treat, formula = formula,
    accept_prob = accept_prob, threshold = threshold, id = id, tol = tol
  )
  rerand_draw(
    spec = spec, max_tries = max_tries, failure_prob = failure_prob,
    seed = seed, engine = engine, on_failure = on_failure
  )
}

#' Draw assignments from a rerandomization specification
#'
#' @param spec A `rerand_spec` object created by [rerand_spec()].
#' @param n_draws Number of assignments to generate.
#' @param max_tries Maximum number of attempts for each assignment.
#' @param failure_prob Failure bound used to choose `max_tries` when omitted.
#' @param seed Optional integer seed. Supplying a seed does not change the
#'   caller's random-number-generator state.
#' @param engine Either `"cpp"` or `"R"`.
#' @param on_failure Whether to error or warn after an unsuccessful draw.
#' @return An object of class `rerand_design_result`.
#' @export
rerand_draw <- function(spec, n_draws = 1L, max_tries = NULL,
                        failure_prob = 1e-6, seed = NULL,
                        engine = c("cpp", "R"),
                        on_failure = c("error", "warn")) {
  if (!inherits(spec, "rerand_spec")) {
    stop("spec must be a rerand_spec object.", call. = FALSE)
  }
  n_draws <- .rerand_validate_n_draws(n_draws)
  max_tries <- .rerand_resolve_max_tries(
    max_tries, spec$criterion$acceptance_mass, failure_prob
  )
  seed <- .rerand_validate_seed(seed)
  engine <- .rerand_validate_engine(match.arg(engine))
  on_failure <- match.arg(on_failure)

  draw_one <- function() {
    if (engine == "R") {
      .design_r(spec$whitened, spec$n_treat, spec$criterion$threshold,
                max_tries)
    } else {
      design_cpp_core(spec$whitened, spec$n_treat,
                      spec$criterion$threshold, max_tries)
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
    spec = spec,
    tries = diagnostics$tries[1L],
    mahalanobis = diagnostics$mahalanobis[1L],
    threshold = spec$criterion$threshold,
    accept_prob = spec$criterion$accept_prob,
    criterion_type = spec$criterion$type,
    accepted = diagnostics$accepted[1L],
    engine = engine,
    n_treat = spec$n_treat,
    n_control = spec$n_control
  )
  class(result) <- c("rerand_design_result", "list")
  result
}

#' @export
print.rerand_design_result <- function(x, ...) {
  cat("Rerandomization design\n")
  cat("  engine:        ", x$engine, "\n", sep = "")
  cat("  criterion:     ", x$criterion_type, "\n", sep = "")
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
summary.rerand_design_result <- function(object, ...) {
  summary <- object[c(
    "accepted", "tries", "mahalanobis", "threshold", "accept_prob",
    "criterion_type", "engine", "n_treat", "n_control"
  )]
  class(summary) <- c("summary.rerand_design_result", "list")
  summary
}

#' @export
print.summary.rerand_design_result <- function(x, ...) {
  print.rerand_design_result(x, ...)
  invisible(x)
}

#' Return design data with a treatment-assignment column
#'
#' @param x A `rerand_design_result` object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param treatment Name of the treatment-assignment column.
#' @param overwrite Whether to replace an existing column named `treatment`.
#' @param ... Ignored.
#' @export
as.data.frame.rerand_design_result <- function(x, row.names = NULL,
                                                optional = FALSE,
                                                treatment = "Z",
                                                overwrite = FALSE, ...) {
  if (length(treatment) != 1L || !is.character(treatment) || is.na(treatment) ||
      !nzchar(treatment)) {
    stop("treatment must be one non-empty column name.", call. = FALSE)
  }
  data <- as.data.frame(x$spec$data)
  if (treatment %in% names(data) && !isTRUE(overwrite)) {
    stop("treatment already exists in the design data; use overwrite = TRUE.",
         call. = FALSE)
  }
  data[[treatment]] <- x$Z
  data
}

#' Summarize covariate balance in a rerandomized design
#'
#' @param x A `rerand_design_result` object.
#' @param ... Ignored.
#' @return A data frame with treatment and control means, mean differences, and
#'   standardized mean differences for each encoded covariate.
#' @export
balance <- function(x, ...) {
  UseMethod("balance")
}

#' @export
balance.rerand_design_result <- function(x, ...) {
  X <- x$spec$X
  Z <- x$Z
  treated <- colMeans(X[Z == 1, , drop = FALSE])
  control <- colMeans(X[Z == 0, , drop = FALSE])
  standard_deviation <- apply(X, 2L, stats::sd)
  data.frame(
    covariate = colnames(X),
    treated_mean = unname(treated),
    control_mean = unname(control),
    difference = unname(treated - control),
    standardized_difference = unname((treated - control) / standard_deviation),
    row.names = NULL
  )
}
