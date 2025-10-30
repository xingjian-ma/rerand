#' Rerandomization Design Function
#'
#' This function implements a rerandomization design based on Mahalanobis distance.
#' @param X An n × K covariate matrix (n rows, K cols).
#' @param n_1 Number of units to assign to treatment.










ReM <- function(X,
                n_1,
                p_a = 0.1,
                max_tries = 10000,
                seed = NULL) {
  # X        : an n × K covariate matrix (n rows, K cols)
  # n_1      : number of units to assign to treatment
  # p_a      : acceptance probability (between 0 and 1)
  # max_tries: maximum number of random draws before giving up
  # seed     : optional random seed for reproducibility

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(X)
  K <- ncol(X)
  n_0 <- n - n_1

  # compute covariance estimate and its inverse
  S   <- cov(X)
  S_inv <- solve(S)

  # approximate threshold a from p_a and chi-square distribution with df = K
  # Since under complete randomization (and large n) M² ~ χ²_K approximately.
  a <- qchisq(p_a, df = K)

  for (t in seq_len(max_tries)) {
    # draw an assignment vector Z
    Z <- sample(c(rep(1, n_1), rep(0, n_0)))

    # compute group means
    Xbar_T <- colMeans(X[Z == 1, , drop = FALSE])
    Xbar_C <- colMeans(X[Z == 0, , drop = FALSE])

    # compute Mahalanobis distance
    diff   <- Xbar_T - Xbar_C
    M      <- as.numeric(t(diff) %*% S_inv %*% diff)

    # check acceptance
    if (M <= a) {
      # accepted
      return(list(Z = Z,
                  tries = t,
                  M     = M,
                  a     = a,
                  p_a   = p_a))
    }
  }

  warning("Maximum tries exceeded without reaching threshold. Returning last assignment anyway.")
  return(list(Z     = Z,
              tries = max_tries,
              M     = M,
              a     = a,
              p_a   = p_a))
}


# Example usage:
# set.seed(123)
X <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)  # 100 units, 3 covariates
result <- ReM(X, n_1 = 50, p_a = 0.1, max_tries = 10000, seed = 42)

