#' Unadjusted Diffference-in-Means Estimator Function
#'
#' Calculates the unadjusted difference-in-means estimator for the average treatment effect (ATE) along with its standard error and 95% confidence interval.
#' @param Y_obs Numeric vector; observed outcomes with length n.
#' @param Z Numeric vector; treatment assignment with length n (0 = control, 1 = treatment).
#'





est_diff_unadj <- function(Y_obs, Z) {
  stopifnot(length(Y) == length(Z))

  Y1 <- Y[Z == 1]; Y0 <- Y[Z == 0]
  n1 <- length(Y1); n0 <- length(Y0)

  tau_hat <- mean(Y1) - mean(Y0)
  se_hat  <- sqrt(var(Y1) / n1 + var(Y0) / n0)
  ci      <- c(tau_hat - 1.96 * se_hat, tau_hat + 1.96 * se_hat)

  list(tau_hat = tau_hat, se = se_hat, ci95 = ci,
       n1 = n1, n0 = n0, method = "unadjusted")
}

# ==========================================================
# Lin (2013) regression-adjusted estimator with HC SE
#   - Centers covariates X (column-wise)
#   - Fits: Y ~ W + Xc + W:Xc  (fully interacted OLS)
#   - Returns the coefficient on W (ATE) with robust SE
# ==========================================================
est_adjusted <- function(Y, Z, X, center = TRUE, hc_type = "HC3") {
  stopifnot(length(Y) == length(Z), nrow(X) == length(Z))

  n <- length(Z)
  X <- as.matrix(X)

  # Center covariates (recommended in Lin, 2013)
  if (center) {
    x_means <- colMeans(X)
    Xc <- sweep(X, 2, x_means, FUN = "-")
  } else {
    Xc <- X
    x_means <- rep(0, ncol(X))
  }

  # Build data.frame with explicit column names
  df <- data.frame(Y = as.numeric(Y), Z = as.numeric(Z), Xc)
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

  # The ATE is the coefficient on Z (main effect)
  if (!"Z" %in% rownames(ct)) {
    stop("Coefficient 'Z' not found (check design matrix / collinearity).")
  }

  tau_hat <- as.numeric(ct["Z", "Estimate"])
  se_hat  <- as.numeric(ct["Z", "Std. Error"])
  ci95    <- c(tau_hat - 1.96 * se_hat, tau_hat + 1.96 * se_hat)

  list(
    tau_hat  = tau_hat,
    se       = se_hat,
    ci95     = ci95,
    method   = "lin",
    centered = center,
    x_means  = x_means,
    hc_type  = hc_type,
    coef_table = ct,
    fit      = fit
  )
}


