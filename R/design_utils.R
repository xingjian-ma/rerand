# Internal rerandomization assignment helpers.

.rerand_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite integer or NULL.", call. = FALSE)
  }
  as.integer(seed)
}

.rerand_with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.rerand_resolve_max_tries <- function(max_tries, acceptance_mass,
                                      failure_prob) {
  if (is.null(max_tries)) {
    if (length(failure_prob) != 1L || !is.numeric(failure_prob) ||
        !is.finite(failure_prob) || failure_prob <= 0 || failure_prob >= 1) {
      stop("failure_prob must be a finite number strictly between 0 and 1.",
           call. = FALSE)
    }
    if (acceptance_mass >= 1) {
      return(1L)
    }
    max_tries <- ceiling(log(failure_prob) / log1p(-acceptance_mass))
  }
  if (length(max_tries) != 1L || !is.numeric(max_tries) ||
      !is.finite(max_tries) || max_tries < 1 ||
      max_tries != as.integer(max_tries)) {
    stop("max_tries must be one positive integer or NULL.", call. = FALSE)
  }
  as.integer(max_tries)
}

.rerand_validate_n_draws <- function(n_draws) {
  if (length(n_draws) != 1L || !is.numeric(n_draws) || !is.finite(n_draws) ||
      n_draws < 1 || n_draws != as.integer(n_draws)) {
    stop("n_draws must be one positive integer.", call. = FALSE)
  }
  as.integer(n_draws)
}

.design_r <- function(X, n_treat, threshold, max_tries) {
  n <- nrow(X)
  n_control <- n - n_treat
  scale_factor <- n_treat * n_control / n
  total <- colSums(X)
  Z <- numeric(n)
  M <- Inf

  for (tries in seq_len(max_tries)) {
    Z[] <- 0
    Z[sample.int(n, n_treat)] <- 1
    treated_total <- colSums(X[Z == 1, , drop = FALSE])
    difference <- treated_total / n_treat - (total - treated_total) / n_control
    M <- as.numeric(crossprod(difference)) * scale_factor
    M <- max(0, M)

    if (M <= threshold) {
      return(list(Z = Z, M = M, tries = tries, accepted = TRUE))
    }
  }

  list(Z = Z, M = M, tries = max_tries, accepted = FALSE)
}
