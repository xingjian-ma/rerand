# Internal estimator and inference implementations.

.calc_sample_stats <- function(Y_obs, Z, X, criterion) {
  assignment <- .rerand_validate_assignment(Z, length(Y_obs), min_group_size = 2L)
  Y1 <- Y_obs[assignment$Z == 1]
  Y0 <- Y_obs[assignment$Z == 0]
  n <- length(Y_obs)
  n1 <- assignment$n_treat
  n0 <- assignment$n_control
  r1 <- n1 / n
  r0 <- n0 / n
  tau_dim <- mean(Y1) - mean(Y0)
  V_tt_hat_1 <- stats::var(Y1) / r1 + stats::var(Y0) / r0

  base <- list(
    n = n,
    n1 = n1,
    n0 = n0,
    r1 = r1,
    r0 = r0,
    tau_dim = as.numeric(tau_dim),
    V_tt_hat_1 = .rerand_clamp_variance(V_tt_hat_1),
    V_tt_hat_2 = NULL,
    R2_hat = NULL,
    v_K_a = NULL,
    correction_factor = 1,
    K = criterion$K,
    criterion_type = criterion$type,
    accept_prob = criterion$accept_prob,
    acceptance_mass = criterion$acceptance_mass,
    threshold = criterion$threshold
  )

  if (is.null(X)) {
    return(base)
  }

  covariance <- .rerand_covariance(X)
  X1 <- X[assignment$Z == 1, , drop = FALSE]
  X0 <- X[assignment$Z == 0, , drop = FALSE]
  S_inv <- covariance$inverse

  S_Y1X <- matrix(stats::cov(Y1, X1), nrow = 1L)
  S_Y0X <- matrix(stats::cov(Y0, X0), nrow = 1L)
  S_tauX <- as.numeric((S_Y1X - S_Y0X) %*% S_inv %*% t(S_Y1X - S_Y0X))
  V_tt_hat_2 <- .rerand_clamp_variance(base$V_tt_hat_1 - S_tauX)

  S_Y1_given_X <- as.numeric(S_Y1X %*% S_inv %*% t(S_Y1X))
  S_Y0_given_X <- as.numeric(S_Y0X %*% S_inv %*% t(S_Y0X))
  numerator <- S_Y0_given_X / r0 + S_Y1_given_X / r1 - S_tauX
  R2_hat <- if (V_tt_hat_2 == 0) {
    if (abs(numerator) < 1e-10) 0 else stop("R2_hat is undefined.", call. = FALSE)
  } else {
    numerator / V_tt_hat_2
  }
  if (!is.finite(R2_hat) || R2_hat < -1e-8 || R2_hat > 1 + 1e-8) {
    stop("The sample R2 estimate is outside its valid range.", call. = FALSE)
  }
  R2_hat <- min(1, max(0, R2_hat))
  correction_factor <- .rerand_correction_factor(R2_hat, criterion)

  base$V_tt_hat_2 <- V_tt_hat_2
  base$R2_hat <- as.numeric(R2_hat)
  base$v_K_a <- criterion$v_K_a
  base$correction_factor <- as.numeric(correction_factor)
  base
}

.estimate_dim <- function(sample_stats) {
  se_neyman <- sqrt(sample_stats$V_tt_hat_1 * sample_stats$correction_factor) /
    sqrt(sample_stats$n)
  se_ding <- if (is.null(sample_stats$V_tt_hat_2)) {
    NULL
  } else {
    sqrt(sample_stats$V_tt_hat_2 * sample_stats$correction_factor) /
      sqrt(sample_stats$n)
  }
  list(
    tau_hat = sample_stats$tau_dim,
    se_neyman = as.numeric(se_neyman),
    se_ding = if (is.null(se_ding)) NULL else as.numeric(se_ding)
  )
}

.estimate_ancova <- function(Y_obs, Z, X) {
  if (is.null(X)) {
    stop("X is required for method 'ancova'.", call. = FALSE)
  }
  X_names <- paste0("X", seq_len(ncol(X)))
  model_data <- data.frame(Y = as.numeric(Y_obs), Z = as.numeric(Z), X)
  names(model_data) <- c("Y", "Z", X_names)
  formula <- stats::as.formula(paste(
    "Y ~ Z +", paste(X_names, collapse = " + ")
  ))
  fit <- stats::lm(formula, data = model_data)
  coefficient <- stats::coef(fit)["Z"]
  if (is.na(coefficient)) {
    stop("The ANCOVA model could not identify the treatment coefficient.",
         call. = FALSE)
  }
  covariance <- sandwich::vcovHC(fit, type = "HC2")
  se_ehw <- sqrt(covariance["Z", "Z"])
  list(tau_hat = coefficient, se_ehw = as.numeric(se_ehw), fit = fit)
}

.estimate_lin <- function(Y_obs, Z, X) {
  if (is.null(X)) {
    stop("X is required for method 'lin'.", call. = FALSE)
  }
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  X_names <- paste0("X", seq_len(ncol(X_centered)))
  data <- data.frame(Y = as.numeric(Y_obs), Z = as.numeric(Z), X_centered)
  names(data) <- c("Y", "Z", X_names)
  formula <- stats::as.formula(paste(
    "Y ~ Z +",
    paste(X_names, collapse = " + "),
    "+",
    paste(paste0("Z:", X_names), collapse = " + ")
  ))
  fit <- stats::lm(formula, data = data)
  coefficient <- stats::coef(fit)["Z"]
  if (is.na(coefficient)) {
    stop("The Lin model could not identify the treatment coefficient.", call. = FALSE)
  }
  covariance <- sandwich::vcovHC(fit, type = "HC2")
  se_ehw <- sqrt(covariance["Z", "Z"])
  list(tau_hat = coefficient, se_ehw = as.numeric(se_ehw), fit = fit)
}

.project_onto_covariates <- function(y, X) {
  design <- cbind(1, X)
  beta <- MASS::ginv(crossprod(design)) %*% crossprod(design, y)
  as.numeric(design %*% beta)
}

.calc_population_stats <- function(Y_full, X, n_treat, criterion) {
  Y0 <- Y_full[, 1]
  Y1 <- Y_full[, 2]
  tau <- Y1 - Y0
  n <- length(Y0)
  n0 <- n - n_treat
  r1 <- n_treat / n
  r0 <- n0 / n
  V_tt_true <- .rerand_clamp_variance(
    stats::var(Y1) / r1 + stats::var(Y0) / r0 - stats::var(tau)
  )
  result <- list(
    tau_true = as.numeric(mean(tau)),
    R2_true = NULL,
    V_tt_true = V_tt_true,
    se_cre_dim = sqrt(V_tt_true) / sqrt(n),
    se_rem_dim = NULL,
    se_cre_lin = NULL,
    se_rem_lin = NULL,
    criterion_type = criterion$type,
    accept_prob = criterion$accept_prob,
    threshold = criterion$threshold
  )

  if (is.null(X)) {
    return(result)
  }

  fitted_Y0 <- .project_onto_covariates(Y0, X)
  fitted_Y1 <- .project_onto_covariates(Y1, X)
  fitted_tau <- fitted_Y1 - fitted_Y0
  R2_true <- (
    stats::var(fitted_Y0) / r0 + stats::var(fitted_Y1) / r1 - stats::var(fitted_tau)
  ) / V_tt_true
  if (!is.finite(R2_true) || R2_true < -1e-8 || R2_true > 1 + 1e-8) {
    stop("The population R2 is outside its valid range.", call. = FALSE)
  }
  R2_true <- min(1, max(0, R2_true))
  correction_factor <- .rerand_correction_factor(R2_true, criterion)
  result$R2_true <- as.numeric(R2_true)
  result$se_rem_dim <- sqrt(correction_factor * V_tt_true) / sqrt(n)
  result$se_cre_lin <- sqrt(.rerand_clamp_variance((1 - R2_true) * V_tt_true)) /
    sqrt(n)
  result$se_rem_lin <- result$se_cre_lin
  result
}

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
