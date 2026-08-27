#' Construct design-aware inference for an estimate
#'
#' @description
#' The assignment design determines the reference distribution. Complete
#' randomization (`cre`) uses a normal critical value. Rerandomization (`rem`)
#' uses the estimator's rerandomization mixture quantile.
#'
#' @param estimate A `rerand_estimate` object.
#' @param level Confidence level, strictly between 0 and 1.
#' @param integration_tol Relative tolerance for the ReM quantile calculation.
#' @return An object of class `rerand_inference`.
#' @export
rerand_inference <- function(estimate, level = 0.95, integration_tol = 1e-8) {
  if (!inherits(estimate, "rerand_estimate")) {
    stop("estimate must be a rerand_estimate object.", call. = FALSE)
  }
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("level must be strictly between 0 and 1.", call. = FALSE)
  }
  if (length(integration_tol) != 1L || !is.numeric(integration_tol) ||
      !is.finite(integration_tol) || integration_tol <= 0 ||
      integration_tol >= 1) {
    stop("integration_tol must be a finite number strictly between 0 and 1.",
         call. = FALSE)
  }

  level <- as.numeric(level)
  alpha <- (1 + level) / 2
  design_method <- estimate$design_method
  if (!design_method %in% c("cre", "rem")) {
    stop("estimate has an unsupported design method.", call. = FALSE)
  }

  r2_stats <- if (design_method == "rem" && estimate$estimator != "dim") {
    assignment <- estimate$assignment
    design_X <- .rerand_design_covariates(assignment, estimate$data)
    .calc_sample_stats(
      Y_obs = estimate$data[[estimate$outcome_name]],
      Z = .rerand_assignment_vector(assignment, estimate$data),
      X = design_X,
      criterion = assignment$design$criterion
    )
  } else {
    estimate$sample_stats
  }

  if (design_method == "cre") {
    reference_distribution <- "normal"
    critical_value <- stats::qnorm(alpha)
    standard_error <- estimate$standard_error
  } else {
    if (is.null(r2_stats) || is.null(r2_stats$R2_hat)) {
      stop("A ReM inference requires a finite sample R2 estimate.",
           call. = FALSE)
    }
    reference_distribution <- "rem"
    criterion <- estimate$assignment$design$criterion
    critical_value <- .rerand_quantile_integration(
      R2 = r2_stats$R2_hat,
      K = criterion$K,
      criterion = criterion,
      alpha = alpha,
      integration_tol = integration_tol
    )
    standard_error <- estimate$unadjusted_standard_error
  }

  interval <- matrix(
    estimate$estimate + c(-1, 1) * critical_value * standard_error,
    nrow = 1L,
    dimnames = list("treatment", c("lower", "upper"))
  )
  result <- list(
    estimate = estimate,
    estimate_value = estimate$estimate,
    estimator = estimate$estimator,
    design_method = design_method,
    reference_distribution = reference_distribution,
    level = level,
    critical_value = as.numeric(critical_value),
    standard_error = as.numeric(standard_error),
    interval = interval,
    R2_hat = if (is.null(r2_stats)) NULL else r2_stats$R2_hat,
    criterion_type = estimate$criterion_type,
    accept_prob = estimate$accept_prob,
    threshold = estimate$threshold
  )
  class(result) <- c("rerand_inference", "list")
  result
}

#' @export
confint.rerand_inference <- function(object, parm, level = NULL, ...) {
  if (!missing(parm) && !identical(parm, "treatment")) {
    stop("Only the treatment coefficient is available.", call. = FALSE)
  }
  if (!is.null(level)) {
    if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
        level <= 0 || level >= 1) {
      stop("level must be strictly between 0 and 1.", call. = FALSE)
    }
    if (!isTRUE(all.equal(as.numeric(level), object$level))) {
      stop("level is fixed when the inference object is created.", call. = FALSE)
    }
  }
  object$interval
}

#' @export
print.rerand_inference <- function(x, ...) {
  cat("Rerandomization inference\n")
  cat("  estimator:     ", x$estimator, "\n", sep = "")
  cat("  design method: ", x$design_method, "\n", sep = "")
  cat("  reference:     ", x$reference_distribution, "\n", sep = "")
  cat("  level:         ", format(x$level), "\n", sep = "")
  cat("  critical value:", format(x$critical_value), "\n", sep = "")
  cat("  interval:      [", format(x$interval[1, 1]), ", ",
      format(x$interval[1, 2]), "]\n", sep = "")
  invisible(x)
}

#' @export
summary.rerand_inference <- function(object, ...) {
  result <- object[c(
    "estimate_value", "estimator", "design_method",
    "reference_distribution", "level", "critical_value", "standard_error",
    "interval", "R2_hat"
  )]
  class(result) <- c("summary.rerand_inference", "list")
  result
}

#' @export
print.summary.rerand_inference <- function(x, ...) {
  print.rerand_inference(x, ...)
  invisible(x)
}

# Internal rerandomization reference-quantile helpers.

.rerand_truncated_d1_density <- function(x, K, threshold, acceptance_mass) {
  inside <- abs(x) <= sqrt(threshold)
  result <- numeric(length(x))
  if (!any(inside)) {
    return(result)
  }
  remaining_mass <- if (K == 1L) {
    rep(1, sum(inside))
  } else {
    stats::pchisq(threshold - x[inside]^2, df = K - 1L)
  }
  result[inside] <- stats::dnorm(x[inside]) * remaining_mass / acceptance_mass
  result
}

.rerand_integrate_d1 <- function(function_to_integrate, K, criterion,
                                 integration_tol) {
  limit <- sqrt(criterion$threshold)
  stats::integrate(
    function(x) {
      function_to_integrate(x) * .rerand_truncated_d1_density(
        x, K = K, threshold = criterion$threshold,
        acceptance_mass = criterion$acceptance_mass
      )
    },
    lower = -limit,
    upper = limit,
    rel.tol = integration_tol,
    subdivisions = 2000L
  )$value
}

.rerand_quantile_integration <- function(R2, K, criterion, alpha,
                                         integration_tol) {
  if (criterion$acceptance_mass == 1) {
    return(stats::qnorm(alpha))
  }
  if (R2 == 0) {
    return(stats::qnorm(alpha))
  }
  if (R2 == 1) {
    limit <- sqrt(criterion$threshold)
    cdf <- function(value) {
      if (value <= -limit) return(0)
      if (value >= limit) return(1)
      .rerand_integrate_d1(
        function(x) as.numeric(x <= value), K, criterion, integration_tol
      )
    }
    return(stats::uniroot(
      function(value) cdf(value) - alpha,
      interval = c(-limit, limit), tol = integration_tol
    )$root)
  }
  standard_deviation <- sqrt(1 - R2)
  cdf <- function(value) {
    .rerand_integrate_d1(
      function(x) stats::pnorm((value - sqrt(R2) * x) / standard_deviation),
      K, criterion, integration_tol
    )
  }
  upper <- max(8, sqrt(R2 * criterion$threshold) + 8 * standard_deviation)
  stats::uniroot(
    function(value) cdf(value) - alpha,
    interval = c(-upper, upper), tol = integration_tol
  )$root
}

.rerand_quantile_simulation <- function(R2, K, criterion, alpha, n_sim,
                                        engine) {
  if (engine == "cpp") {
    return(as.numeric(get_quantile_cpp(
      R2 = R2, K = K, a = criterion$threshold, alpha = alpha, n_sim = n_sim
    )))
  }
  values <- numeric(n_sim)
  accepted <- 0L
  while (accepted < n_sim) {
    batch_size <- max(1000L, n_sim - accepted)
    d1 <- stats::rnorm(batch_size)
    remaining <- if (K > 1L) stats::rchisq(batch_size, K - 1L) else 0
    keep <- d1^2 + remaining <= criterion$threshold
    if (any(keep)) {
      values_to_add <- sqrt(R2) * d1[keep] +
        sqrt(1 - R2) * stats::rnorm(sum(keep))
      take <- min(length(values_to_add), n_sim - accepted)
      values[(accepted + 1L):(accepted + take)] <-
        abs(values_to_add[seq_len(take)])
      accepted <- accepted + take
    }
  }
  as.numeric(stats::quantile(values, probs = 2 * alpha - 1,
                             names = FALSE, type = 7))
}

.rerand_quantile <- function(R2, K, accept_prob = NULL, threshold = NULL,
                             alpha = 0.975,
                             method = c("integration", "simulation"),
                             integration_tol = 1e-8, n_sim = 100000L,
                             seed = NULL, engine = c("cpp", "R")) {
  inputs <- .rerand_validate_simulation_inputs(R2, K, alpha, n_sim)
  if (length(integration_tol) != 1L || !is.numeric(integration_tol) ||
      !is.finite(integration_tol) || integration_tol <= 0 ||
      integration_tol >= 1) {
    stop("integration_tol must be a finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  method <- match.arg(method)
  engine <- .rerand_validate_engine(match.arg(engine))
  criterion <- .rerand_resolve_criterion(
    accept_prob = accept_prob, threshold = threshold, K = inputs$K,
    require_criterion = TRUE
  )
  if (method == "integration") {
    return(.rerand_quantile_integration(
      inputs$R2, inputs$K, criterion, inputs$alpha, integration_tol
    ))
  }
  .rerand_with_seed(
    .rerand_validate_seed(seed),
    .rerand_quantile_simulation(
      inputs$R2, inputs$K, criterion, inputs$alpha, inputs$n_sim, engine
    )
  )
}
