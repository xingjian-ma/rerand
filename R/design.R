#' Rerandomization Design Function
#'
#'@description
#' Provides high-performance tools for rerandomization in randomized
#' experiments, including Mahalanobis-distance-based acceptance rules,
#' efficient generation of treatment assignments, and covariate-balance
#' diagnostics. Core algorithms are implemented in Rcpp and
#' RcppArmadillo for computational efficiency.
#'
#' This function implements a rerandomization design based on Mahalanobis distance.
#' @param X Numerical matrix; covariate matrix with n rows (units) and K columns (covariates). n >= 2 and K >= 1.
#' @param Y Numerical matrix; potential outcomes matrix with n rows and 2 columns (Y(0), Y(1)).
#' @param n_1 Integer scalar; number of units to assign to treatment. n_1 > 0 and n_1 < nrow(X).
#' @param p_accept Numeric scalar; acceptance probability in (0, 1], default is 0.1.
#' @param threshold Numeric scalar; threshold for Mahalanobis distance. If NULL, it is computed from p_accept. Default is NULL. If provided, p_accept will be ignored.
#' @param max_tries Integer scalar; maximum number of random draws before giving up, default is 10000.
#' @param seed Integer scalar; optional random seed for reproducibility, default is NULL.
#' @param engine Character scalar; computation engine to use, either "R" or "cpp". Default is "cpp".
#'
#' @return A list containing:
#' \describe{
#'   \item{Z}{Numeric vector; accepted treatment assignment with length n and values 0 (control) or 1 (treatment).}
#'   \item{Y_obs}{Numeric vector; observed outcomes corresponding to the accepted assignment.}
#'   \item{tries}{Integer scalar; number of random draws made until acceptance.}
#'   \item{M}{Numeric scalar; Mahalanobis dsistance of the accepted assignment.}
#'   \item{threshold}{Numeric scalar; threshold Mahalanobis distance used for acceptance.}
#'   \item{p_accept}{Numeric scalar; acceptance probability used.}
#'   \item{accepted}{Logical scalar; indicates whether an acceptable assignment was found within max_tries.}
#'   \item{engine}{Character scalar; computation engine used ("R" or "cpp").}
#' }
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)  # 100 units, 3 covariates
#' Y <- matrix(rnorm(100 * 2), nrow = 100, ncol = 2)  # potential outcomes
#' result <- ReM(X, Y, n_1 = 50, p_a = 0.1, max_tries = 10000)
#'
#' @export
ReM <- function(X,
                Y,
                n_1,
                p_accept = 0.1,
                threshold = NULL,
                max_tries = 10000,
                seed = NULL,
                engine = "cpp") {

  # Check inputs
  checkmate::assert_matrix(X, mode = "numeric", min.rows = 2, min.cols = 1, any.missing = FALSE)
  checkmate::assert_matrix(Y, mode = "numeric", ncols = 2, nrows = nrow(X), any.missing = FALSE)
  checkmate::assert_count(n_1)


  p_a <- p_accept
  a <- threshold

  checkmate::assert_numeric(p_a, lower = 0, upper = 1, len = 1, any.missing = FALSE)
  checkmate::assert_true(p_a > 0)


  if (!is.null(a)){
    checkmate::assert_numeric(a, lower = 0, len = 1, any.missing = FALSE)
    checkmate::assert_true(a > 0)
  }else{
    # compute threshold a from p_a
    K <- ncol(X)
    a <- stats::qchisq(p = p_a, df = K)
  }

  checkmate::assert_count(max_tries)
  checkmate::assert_count(seed, null.ok = TRUE)

  n_1 <- as.integer(n_1)
  max_tries <- as.integer(max_tries)

  checkmate::assert_integer(n_1, lower = 1, upper = nrow(X) - 1, len = 1, any.missing = FALSE)
  checkmate::assert_integer(max_tries, lower = 1, len = 1, any.missing = FALSE)

  # Set parameters
  if (!is.null(seed)) set.seed(seed)

  # get dimensions
  n <- nrow(X)
  K <- ncol(X)
  n_0 <- n - n_1

  # run core computation
  if (engine == "R") {
    res <- rem_core_R(X = X,
                      Y = Y,
                      n_1 = n_1,
                      a = a,
                      max_tries = max_tries)
  } else if (engine == "cpp") {
    res <- rem_core_cpp(X = X,
                        Y = Y,
                        n1 = n_1,
                        a = a,
                        max_tries = max_tries)
  }

  # return results
  return(list(Z = res$Z,
              Y_obs = res$Y_obs,
              tries = res$tries,
              M = res$M,
              threshold = a,
              p_accept = p_a,
              accepted = res$accepted,
              engine = engine))

}

