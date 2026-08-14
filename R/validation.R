# Internal validation helpers shared by the public APIs.

.rerand_assert_finite_numeric <- function(x, name, length = NULL) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(sprintf("%s must be a finite numeric vector.", name), call. = FALSE)
  }
  if (!is.null(length) && length(x) != length) {
    stop(sprintf("%s must have length %d.", name, length), call. = FALSE)
  }
  invisible(x)
}

.rerand_validate_matrix <- function(X, name = "X", n = NULL) {
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop(sprintf("%s must be a finite numeric matrix.", name), call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L) {
    stop(sprintf("%s must have at least two rows and one column.", name), call. = FALSE)
  }
  if (!is.null(n) && nrow(X) != n) {
    stop(sprintf("%s must have %d rows.", name, n), call. = FALSE)
  }
  X
}

.rerand_validate_n_treat <- function(n_treat, n) {
  if (length(n_treat) != 1L || is.na(n_treat) ||
      !is.numeric(n_treat) || !is.finite(n_treat) ||
      n_treat != as.integer(n_treat)) {
    stop("n_treat must be one integer.", call. = FALSE)
  }
  n_treat <- as.integer(n_treat)
  if (n_treat < 1L || n_treat >= n) {
    stop("n_treat must be between 1 and the number of units minus 1.",
         call. = FALSE)
  }
  n_treat
}

.rerand_validate_tol <- function(tol) {
  if (length(tol) != 1L || !is.numeric(tol) || !is.finite(tol) ||
      tol <= 0 || tol >= 1) {
    stop("tol must be one finite number strictly between 0 and 1.",
         call. = FALSE)
  }
  as.numeric(tol)
}

.rerand_validate_engine <- function(engine) {
  if (length(engine) != 1L || is.na(engine) || !engine %in% c("cpp", "R")) {
    stop("engine must be one of 'cpp' or 'R'.", call. = FALSE)
  }
  engine
}

.rerand_validate_assignment <- function(Z, n, min_group_size = 1L) {
  .rerand_assert_finite_numeric(Z, "Z", n)
  if (any(!Z %in% c(0, 1))) {
    stop("Z must contain only 0 and 1.", call. = FALSE)
  }
  n_treat <- sum(Z == 1)
  n_control <- sum(Z == 0)
  if (n_treat < min_group_size || n_control < min_group_size) {
    stop("Z must contain enough observations in both treatment groups.", call. = FALSE)
  }
  list(Z = as.numeric(Z), n_treat = n_treat, n_control = n_control)
}

.rerand_validate_simulation_inputs <- function(R2, K, alpha, n_sim) {
  if (length(R2) != 1L || !is.numeric(R2) || !is.finite(R2) || R2 < 0 || R2 > 1) {
    stop("R2 must be a finite number between 0 and 1.", call. = FALSE)
  }
  if (length(K) != 1L || !is.numeric(K) || !is.finite(K) || K < 1 ||
      K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
      alpha <= 0.5 || alpha >= 1) {
    stop("alpha must be strictly between 0.5 and 1.", call. = FALSE)
  }
  if (length(n_sim) != 1L || !is.numeric(n_sim) || !is.finite(n_sim) ||
      n_sim < 1 || n_sim != as.integer(n_sim)) {
    stop("n_sim must be a positive integer.", call. = FALSE)
  }
  list(R2 = as.numeric(R2), K = as.integer(K), alpha = as.numeric(alpha),
       n_sim = as.integer(n_sim))
}
