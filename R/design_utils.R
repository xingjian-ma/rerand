#' Internal R Reference Implementation for Rerandomization Core
#'
#' @param X Numeric matrix of covariates.
#' @param n1 Integer number of treated units.
#' @param a Numeric Mahalanobis threshold.
#' @param max_tries Integer maximum number of draws.
#'
#' @return A list with accepted assignment and diagnostics.
#' @keywords internal
#' @export
design_R <- function(X, n1, a, max_tries) {

  n <- nrow(X)
  n0 <- n - n1

  S_inv <- solve(stats::cov(X))

  Z <- rep(0, n)
  M <- Inf
  accepted <- FALSE

  for (t in seq_len(max_tries)) {
    Z <- sample(c(rep(1, n1), rep(0, n0)))

    Xbar_1 <- colMeans(X[Z == 1, , drop = FALSE])
    Xbar_0 <- colMeans(X[Z == 0, , drop = FALSE])
    diff <- Xbar_1 - Xbar_0

    M <- as.numeric(t(diff) %*% S_inv %*% diff) * (n1 * n0 / n)

    if (M <= a) {
      accepted <- TRUE

      return(list(
        Z = Z,
        M = M,
        tries = t,
        accepted = accepted
      ))
    }
  }

  warning("Maximum tries exceeded without reaching threshold. Returning last assignment anyway.")
  return(list(
    Z = Z,
    M = M,
    tries = max_tries,
    accepted = accepted
  ))
}
