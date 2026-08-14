#' Estimate an average treatment effect under rerandomization
#'
#' @param Y_obs Observed outcomes.
#' @param Z Binary treatment assignment.
#' @param X Numeric covariate matrix. Required for `method = "lin"` and for
#'   rerandomization adjustment.
#' @param method Either `"dim"` for difference-in-means or `"lin"` for Lin's
#'   fully interacted regression adjustment.
#' @param accept_prob Acceptance probability. Mutually exclusive with
#'   `threshold`.
#' @param threshold Direct Mahalanobis threshold. Mutually exclusive with
#'   `accept_prob`.
#' @param theoretical Whether to calculate population-level quantities.
#' @param Y_full Potential outcomes matrix with columns Y(0) and Y(1).
#' @param p_accept Deprecated alias for `accept_prob`.
#' @return An object of class `rerand_estimate_result`.
#' @export
rerand_estimate <- function(Y_obs, Z, X = NULL,
                            method = c("dim", "lin"),
                            accept_prob = NULL, threshold = NULL,
                            theoretical = FALSE, Y_full = NULL,
                            p_accept = NULL) {
  if (!is.null(p_accept)) {
    if (!is.null(accept_prob)) {
      stop("Use only one of accept_prob and p_accept.", call. = FALSE)
    }
    warning("p_accept is deprecated; use accept_prob instead.", call. = FALSE)
    accept_prob <- p_accept
  }
  method <- match.arg(method)
  inputs <- .rerand_validate_estimate_inputs(
    Y_obs = Y_obs,
    Z = Z,
    X = X,
    method = method,
    theoretical = theoretical,
    Y_full = Y_full
  )

  if (is.null(X)) {
    if (!is.null(threshold) || is.null(accept_prob) || accept_prob != 1) {
      stop("Without X, accept_prob must be explicitly set to 1.", call. = FALSE)
    }
    criterion <- list(
      type = "probability",
      accept_prob = 1,
      threshold = Inf,
      K = 0L,
      acceptance_mass = 1,
      v_K_a = 1
    )
  } else {
    covariance <- .rerand_covariance(X)
    criterion <- .rerand_resolve_criterion(
      accept_prob = accept_prob,
      threshold = threshold,
      K = covariance$rank,
      require_criterion = TRUE
    )
  }

  sample_stats <- .calc_sample_stats(
    Y_obs = inputs$Y_obs,
    Z = inputs$assignment$Z,
    X = inputs$X,
    criterion = criterion
  )

  estimate <- if (method == "dim") {
    .estimate_dim(sample_stats)
  } else {
    .estimate_lin(inputs$Y_obs, inputs$assignment$Z, inputs$X)
  }

  population_stats <- if (theoretical) {
    .calc_population_stats(
      Y_full = inputs$Y_full,
      X = inputs$X,
      n_treat = inputs$assignment$n_treat,
      criterion = criterion
    )
  } else {
    NULL
  }

  result <- c(
    list(
      tau_hat = as.numeric(estimate$tau_hat),
      method = method,
      se_neyman = if (method == "dim") estimate$se_neyman else NULL,
      se_ding = if (method == "dim") estimate$se_ding else NULL,
      se_ehw = if (method == "lin") estimate$se_ehw else NULL,
      fit = if (method == "lin") estimate$fit else NULL
    ),
    list(
      sample_stats = sample_stats,
      pop_stats = population_stats,
      criterion_type = criterion$type,
      accept_prob = criterion$accept_prob,
      threshold = criterion$threshold
    )
  )
  class(result) <- c("rerand_estimate_result", "list")
  result
}

#' @export
print.rerand_estimate_result <- function(x, ...) {
  cat("Rerandomization estimate\n")
  cat("  method:        ", x$method, "\n", sep = "")
  cat("  estimate:      ", format(x$tau_hat), "\n", sep = "")
  if (x$method == "dim") {
    cat("  Neyman SE:     ", format(x$se_neyman), "\n", sep = "")
    if (!is.null(x$se_ding)) {
      cat("  Ding SE:       ", format(x$se_ding), "\n", sep = "")
    }
  } else {
    cat("  HC2 SE:        ", format(x$se_ehw), "\n", sep = "")
  }
  cat("  criterion:     ", x$criterion_type, "\n", sep = "")
  invisible(x)
}

#' @export
summary.rerand_estimate_result <- function(object, ...) {
  summary <- object[c(
    "tau_hat", "method", "se_neyman", "se_ding", "se_ehw",
    "criterion_type", "accept_prob", "threshold"
  )]
  class(summary) <- c("summary.rerand_estimate_result", "list")
  summary
}

#' @export
print.summary.rerand_estimate_result <- function(x, ...) {
  print.rerand_estimate_result(x, ...)
  invisible(x)
}
