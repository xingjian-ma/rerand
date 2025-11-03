#' Unadjusted Diffference-in-Means Estimator Function
#'
#' Calculates the unadjusted difference-in-means estimator for the average treatment effect (ATE) along with its standard error and 95% confidence interval.
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment with length n (0 = control, 1 = treatment).
#' 
#' @return A list containing:
#' \describe{
#' \item{tau_hat}{Numeric scalar; estimated ATE.}
#' \item{se}{Numeric scalar; standard error of the estimate.}
#' \item{ci}{Numeric vector; 95% confidence interval for the ATE.}
#' \item{method}{Character scalar; method used ("difference-in-means unadjusted").}
#' }
#'
#' @examples
#' res <- est_diff_unadj(Y_obs = rnorm(100), Z = rbinom(100, 1, 0.5))
#' 
#' @export
est_diff_unadj <- function(Y_obs, Z) {
  
  # Check inputs
  checkmate::assert_numeric(Y_obs, len = length(Z), any.missing = FALSE)
  checkmate::assert_numeric(Z, len = length(Y_obs), any.missing = FALSE)
  checkmate::assert_subset(unique(Z), choices = c(0, 1))
  
  # Calculate difference-in-means estimator
  Y1 <- Y_obs[Z == 1]
  Y0 <- Y_obs[Z == 0]
  n1 <- length(Y1)
  n0 <- length(Y0)

  tau_hat <- mean(Y1) - mean(Y0)
  se_hat  <- sqrt(var(Y1) / n1 + var(Y0) / n0)
  ci      <- c(tau_hat - 1.96 * se_hat, tau_hat + 1.96 * se_hat)

  return(list(
    tau_hat = tau_hat,
    se      = se_hat,
    ci      = ci,
    method  = "difference-in-means unadjusted"
  ))
}


#' Lin (2013) Regression-Adjusted Estimator Function
#' 
#' Calculates the Lin (2013) regression-adjusted estimator for the average treatment effect (ATE) using a fully interacted linear model with covariates, along with heteroskedasticity-consistent standard errors.
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment with length0 = control, 1 = treatment).
#' @param X Numeric matrix; n × K covariate matrix.
#' 
#' @return A list containing:
#' \describe{
#' \item{tau_hat}{Numeric scalar; estimated ATE.}
#' \item{se}{Numeric scalar; standard error of the estimate.}
#' \item{ci{Numeric vector; 95% confidence interval for the ATE.}
#' \item{method}{Character scalar; method used ("lin adjusted").}
#' \item{fit}{lm object; fitted linear model.}
#' }
#' 
#' @examples
#' X <- matrix(rnorm(200), nrow = 100, ncol = 2)
#' res <- est_adjusted(Y_obs = rnorm(100), Z = rbinom(100, 1, 0.5), X = X)
#' 
#' @export
est_lin_adjusted <- function(Y_obs, Z, X, hc_type = "HC3") {
  
  # Check inputs
  checkmate::assert_numeric(Y_obs, len = length(Z), any.missing = FALSE)
  checkmate::assert_numeric(Z, len = length(Y_obs), any.missing = FALSE)
  checkmate::assert_subset(unique(Z), choices = c(0, 1))
  checkmate::assert_matrix(X, mode = "numeric", nrows = length(Z), any.missing = FALSE)

  n <- length(Z)
  X_means <- colMeans(X)
  Xc <- sweep(X, 2, X_means, FUN = "-")  # Center covariates
  
  # Build data.frame with explicit column names
  df <- data.frame(Y = as.numeric(Y_obs), Z = as.numeric(Z), Xc)
  colnames(df) <- c("Y", "Z", paste0("X.", seq_len(ncol(Xc))))

  # Lin model: Y ~ Z + Xc + Z:Xc  (fully interacted)
  rhs <- paste0("Z + ",
                paste0("X.", seq_len(ncol(Xc)), collapse = " + "),
                " + ",
                paste0("Z:X.", seq_len(ncol(Xc)), collapse = " + "))
  fm  <- stats::as.formula(paste0("Y ~ ", rhs))

  fit <- stats::lm(fm, data = df)

  # Robust (heteroskedasticity-consistent) SE, default HC3
  Vhc <- sandwich::vcovHC(fit, type = hc_type)
  ct  <- lmtest::coeftest(fit, vcov. = Vhc)


  tau_hat <- as.numeric(ct["Z", "Estimate"])
  se_hat  <- as.numeric(ct["Z", "Std. Error"])
  ci    <- c(tau_hat - 1.96 * se_hat, tau_hat + 1.96 * se_hat)

  list(
    tau_hat  = tau_hat,
    se       = se_hat,
    ci       = ci,
    method   = "lin adjusted",
    fit      = fit
  )
}


