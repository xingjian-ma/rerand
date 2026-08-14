# Internal R reference implementation for the rerandomization core.

.design_r <- function(X, n_treat, threshold, max_tries, S_inv) {
  n <- nrow(X)
  n_control <- n - n_treat
  scale_factor <- n_treat * n_control / n
  Z <- numeric(n)
  M <- Inf

  for (tries in seq_len(max_tries)) {
    Z[] <- 0
    Z[sample.int(n, n_treat)] <- 1
    Xbar_treat <- colMeans(X[Z == 1, , drop = FALSE])
    Xbar_control <- colMeans(X[Z == 0, , drop = FALSE])
    difference <- Xbar_treat - Xbar_control
    M <- as.numeric(crossprod(difference, S_inv %*% difference)) * scale_factor
    M <- max(0, M)

    if (M <= threshold) {
      return(list(Z = Z, M = M, tries = tries, accepted = TRUE))
    }
  }

  warning(
    "Maximum tries exceeded without reaching threshold; returning the last assignment.",
    call. = FALSE
  )
  list(Z = Z, M = M, tries = max_tries, accepted = FALSE)
}
