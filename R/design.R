#' Generate a rerandomized treatment assignment
#'
#' @description
#' Generates an assignment with exactly `n_treat` treated units and accepts it
#' when its Mahalanobis distance is below a prespecified criterion.
#'
#' @param X Numeric covariate matrix.
#' @param n_treat Number of treated units.
#' @param accept_prob Target acceptance probability used to derive the
#'   Mahalanobis threshold. Mutually exclusive with `threshold`.
#' @param threshold Direct Mahalanobis threshold. Mutually exclusive with
#'   `accept_prob`.
#' @param max_tries Maximum number of assignments to draw.
#' @param seed Optional integer seed.
#' @param engine Either `"cpp"` or `"R"`.
#' @param p_accept Deprecated alias for `accept_prob`.
#' @return An object of class `rerand_design_result`.
#' @export
rerand_design <- function(X, n_treat, accept_prob = NULL, threshold = NULL,
                          max_tries = 10000L, seed = NULL,
                          engine = c("cpp", "R"), p_accept = NULL) {
  if (!is.null(p_accept)) {
    if (!is.null(accept_prob)) {
      stop("Use only one of accept_prob and p_accept.", call. = FALSE)
    }
    warning("p_accept is deprecated; use accept_prob instead.", call. = FALSE)
    accept_prob <- p_accept
  }
  inputs <- .rerand_validate_design_inputs(X, n_treat, max_tries)
  engine <- .rerand_validate_engine(match.arg(engine))
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
        seed != as.integer(seed)) {
      stop("seed must be one finite integer or NULL.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  covariance <- .rerand_covariance(inputs$X)
  criterion <- .rerand_resolve_criterion(
    accept_prob = accept_prob,
    threshold = threshold,
    K = covariance$rank,
    require_criterion = TRUE
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }

  core <- if (engine == "R") {
    .design_r(
      X = inputs$X,
      n_treat = inputs$n_treat,
      threshold = criterion$threshold,
      max_tries = inputs$max_tries,
      S_inv = covariance$inverse
    )
  } else {
    design_cpp_core(
      X = inputs$X,
      n1 = inputs$n_treat,
      a = criterion$threshold,
      max_tries = inputs$max_tries,
      S_inv = covariance$inverse
    )
  }

  result <- list(
    Z = as.numeric(core$Z),
    tries = as.integer(core$tries),
    mahalanobis = as.numeric(core$M),
    threshold = criterion$threshold,
    accept_prob = criterion$accept_prob,
    criterion_type = criterion$type,
    accepted = isTRUE(core$accepted),
    engine = engine,
    n_treat = inputs$n_treat,
    n_control = nrow(inputs$X) - inputs$n_treat
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
