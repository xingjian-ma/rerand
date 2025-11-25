# internal R reference implementation
rem_core_R <- function(X, Y, n_1, a, max_tries) {

  n <- nrow(X)
  K <- ncol(X)
  n_0 <- n - n_1
  
  S_inv <- solve(stats::cov(X))
  
  Z <- rep(0, n)
  M <- Inf
  accepted <- FALSE

  for (t in seq_len(max_tries)) {
    Z <- sample(c(rep(1, n_1), rep(0, n_0)))
    
    Xbar_1 <- colMeans(X[Z == 1, , drop = FALSE])
    Xbar_0 <- colMeans(X[Z == 0, , drop = FALSE])
    diff <- Xbar_1 - Xbar_0
    
    M <- as.numeric(t(diff) %*% S_inv %*% diff) * (n_1 * n_0 / n)
    
    if (M <= a) {
      accepted <- TRUE
      
      Y_obs <- Y[, 1] * (1 - Z) + Y[, 2] * Z

      return(list(Z = Z,
                  Y_obs = Y_obs,
                  M = M,
                  tries = t,
                  accepted = accepted))
    }
  }
  
  warning("Maximum tries exceeded without reaching threshold. Returning last assignment anyway.")
  Y_obs <- Y[, 1] * (1 - Z) + Y[, 2] * Z
  return(list(Z = Z,
              Y_obs = Y_obs,
              M = M,
              tries = max_tries,
              accepted = accepted))
}
