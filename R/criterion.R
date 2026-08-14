# Internal criterion and covariance helpers.

.rerand_covariance <- function(X, tol = 1e-10) {
  X <- .rerand_validate_matrix(X)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  rank <- qr(X_centered, tol = tol)$rank
  if (rank < 1L) {
    stop("X must contain at least one non-constant covariate direction.", call. = FALSE)
  }
  covariance <- stats::cov(X)
  inverse <- MASS::ginv(covariance, tol = tol)
  list(
    X = X,
    centered = X_centered,
    covariance = covariance,
    inverse = inverse,
    rank = rank
  )
}

.rerand_resolve_criterion <- function(accept_prob = NULL, threshold = NULL,
                                     K, require_criterion = TRUE) {
  if (length(K) != 1L || !is.numeric(K) || K < 1 || K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  K <- as.integer(K)
  if (is.null(accept_prob) && is.null(threshold)) {
    if (require_criterion) {
      stop("Exactly one of accept_prob or threshold must be supplied.", call. = FALSE)
    }
    accept_prob <- 1
  }
  if (!is.null(accept_prob) && !is.null(threshold)) {
    stop("accept_prob and threshold are mutually exclusive.", call. = FALSE)
  }
  if (!is.null(accept_prob)) {
    if (length(accept_prob) != 1L || !is.numeric(accept_prob) ||
        !is.finite(accept_prob) || accept_prob <= 0 || accept_prob > 1) {
      stop("accept_prob must be a finite number in (0, 1].", call. = FALSE)
    }
    threshold <- stats::qchisq(accept_prob, df = K)
    acceptance_mass <- accept_prob
    type <- "probability"
  } else {
    if (length(threshold) != 1L || !is.numeric(threshold) ||
        !is.finite(threshold) || threshold <= 0) {
      stop("threshold must be a finite positive number.", call. = FALSE)
    }
    acceptance_mass <- stats::pchisq(threshold, df = K)
    if (!is.finite(acceptance_mass) || acceptance_mass <= 0) {
      stop("threshold does not define a positive acceptance region.", call. = FALSE)
    }
    type <- "threshold"
  }
  v_K_a <- if (acceptance_mass == 1) {
    1
  } else {
    stats::pchisq(threshold, df = K + 2) / acceptance_mass
  }
  list(
    type = type,
    accept_prob = if (type == "probability") as.numeric(accept_prob) else NULL,
    threshold = as.numeric(threshold),
    K = K,
    acceptance_mass = as.numeric(acceptance_mass),
    v_K_a = as.numeric(v_K_a)
  )
}

.rerand_correction_factor <- function(R2, criterion) {
  factor <- 1 - (1 - criterion$v_K_a) * R2
  if (!is.finite(factor) || factor < -1e-10) {
    stop("The rerandomization correction factor is invalid.", call. = FALSE)
  }
  max(0, factor)
}

.rerand_clamp_variance <- function(x, tolerance = 1e-10) {
  if (!is.finite(x) || x < -tolerance) {
    stop("A variance calculation produced an invalid value.", call. = FALSE)
  }
  max(0, as.numeric(x))
}
