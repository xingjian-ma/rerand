#' Deprecated compatibility wrappers
#'
#' These functions are retained for callers of previous development versions.
#' New code should use `rerand_design()`, `rerand_estimate()`, and
#' `get_quantile()`.
#'
#' @param X Numeric covariate matrix.
#' @param n1 Number of treated units.
#' @param p_accept Deprecated acceptance-probability argument.
#' @param threshold Direct Mahalanobis threshold.
#' @param max_tries Maximum number of assignments to draw.
#' @param seed Optional integer seed.
#' @param engine Computation engine.
#' @param a Direct Mahalanobis threshold for low-level wrappers.
#' @param Y_obs Observed outcomes.
#' @param Z Binary treatment assignment.
#' @param method Estimation method.
#' @param X Numeric covariate matrix.
#' @param ... Additional arguments accepted for compatibility.
#' @param R2 Covariate-explained variance proportion.
#' @param K Number of covariate directions.
#' @param alpha Quantile level.
#' @param n_sim Number of simulation draws.
#' @param theoretical Whether to calculate theoretical statistics.
#' @param Y_full Potential outcomes matrix.
#' @name rerand-compatibility
#' @keywords internal
NULL

.rerand_legacy_design_result <- function(result, p_accept = NULL) {
  list(
    Z = result$Z,
    tries = result$tries,
    M = result$mahalanobis,
    threshold = result$threshold,
    p_accept = if (is.null(p_accept)) result$accept_prob else p_accept,
    accepted = result$accepted,
    engine = result$engine
  )
}

#' @rdname rerand-compatibility
#' @export
rerand.design <- function(X, n1, p_accept = 0.1, threshold = NULL,
                          max_tries = 10000, seed = NULL, engine = "cpp") {
  .Deprecated("rerand_design")
  if (is.null(threshold)) {
    result <- rerand_design(
      X = X, n_treat = n1, accept_prob = p_accept,
      max_tries = max_tries, seed = seed, engine = engine
    )
  } else {
    result <- rerand_design(
      X = X, n_treat = n1, threshold = threshold,
      max_tries = max_tries, seed = seed, engine = engine
    )
  }
  .rerand_legacy_design_result(result, p_accept = p_accept)
}

#' @rdname rerand-compatibility
#' @export
design <- function(X, n1, p_accept = 0.1, threshold = NULL,
                   max_tries = 10000, seed = NULL, engine = "cpp") {
  .Deprecated("rerand_design")
  rerand.design(
    X = X, n1 = n1, p_accept = p_accept, threshold = threshold,
    max_tries = max_tries, seed = seed, engine = engine
  )
}

#' @rdname rerand-compatibility
#' @export
design_R <- function(X, n1, a, max_tries) {
  .Deprecated("rerand_design")
  inputs <- .rerand_validate_design_inputs(X, n1, max_tries)
  covariance <- .rerand_covariance(inputs$X)
  .design_r(
    X = inputs$X,
    n_treat = inputs$n_treat,
    threshold = a,
    max_tries = inputs$max_tries,
    S_inv = covariance$inverse
  )
}

#' @rdname rerand-compatibility
#' @export
design.R <- function(X, n1, a, max_tries) {
  .Deprecated("rerand_design")
  design_R(X = X, n1 = n1, a = a, max_tries = max_tries)
}

#' @rdname rerand-compatibility
#' @export
design_cpp <- function(X, n1, a, max_tries) {
  .Deprecated("rerand_design")
  inputs <- .rerand_validate_design_inputs(X, n1, max_tries)
  covariance <- .rerand_covariance(inputs$X)
  design_cpp_core(
    X = inputs$X,
    n1 = inputs$n_treat,
    a = a,
    max_tries = inputs$max_tries,
    S_inv = covariance$inverse
  )
}

#' @rdname rerand-compatibility
#' @export
estimate.dim <- function(Y_obs, Z, X = NULL, p_accept = 1, ...) {
  .Deprecated("rerand_estimate")
  est_dim(Y_obs = Y_obs, Z = Z, X = X, p_accept = p_accept, ...)
}

#' @rdname rerand-compatibility
#' @export
estimate.lin <- function(Y_obs, Z, X, ...) {
  .Deprecated("rerand_estimate")
  est_lin(Y_obs = Y_obs, Z = Z, X = X, ...)
}

#' @rdname rerand-compatibility
#' @export
est_dim <- function(Y_obs, Z, X = NULL, p_accept = 1, ...) {
  .Deprecated("rerand_estimate")
  result <- rerand_estimate(
    Y_obs = Y_obs, Z = Z, X = X, method = "dim", accept_prob = p_accept
  )
  result[c("tau_hat", "se_neyman", "se_ding", "sample_stats")]
}

#' @rdname rerand-compatibility
#' @export
est_lin <- function(Y_obs, Z, X, ...) {
  .Deprecated("rerand_estimate")
  result <- rerand_estimate(
    Y_obs = Y_obs, Z = Z, X = X, method = "lin", accept_prob = 1
  )
  result[c("tau_hat", "se_ehw", "fit")]
}

#' @rdname rerand-compatibility
#' @export
get.quantile <- function(R2, K, p_accept, alpha = 0.975, n_sim = 1e5,
                         engine = "cpp") {
  .Deprecated("get_quantile")
  get_quantile(
    R2 = R2, K = K, accept_prob = p_accept, alpha = alpha,
    n_sim = n_sim, engine = engine
  )
}

#' @rdname rerand-compatibility
#' @export
rerand.estimate <- function(Y_obs, Z, X = NULL, method = c("dim", "lin"),
                            p_accept = 1, theoretical = FALSE, Y_full = NULL) {
  .Deprecated("rerand_estimate")
  rerand_estimate(
    Y_obs = Y_obs, Z = Z, X = X, method = method,
    accept_prob = p_accept, theoretical = theoretical, Y_full = Y_full
  )
}
