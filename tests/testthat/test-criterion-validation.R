test_that("design criterion is explicit and mutually exclusive", {
  X <- matrix(seq_len(40), nrow = 20, ncol = 2)
  X[, 2] <- X[, 2] + rev(X[, 1])

  expect_error(
    rerand_design(X, n_treat = 10),
    "Exactly one of accept_prob or threshold"
  )
  expect_error(
    rerand_design(X, n_treat = 10, accept_prob = 0.5, threshold = 1),
    "mutually exclusive"
  )
  expect_error(
    rerand_design(X, n_treat = 10, accept_prob = 0),
    "accept_prob"
  )
  expect_error(
    rerand_design(X, n_treat = 10, threshold = 0),
    "threshold"
  )
})

test_that("rank-deficient covariates are handled consistently", {
  set.seed(10)
  x <- rnorm(40)
  X <- cbind(x, 2 * x)
  result <- rerand_design(X, n_treat = 20, accept_prob = 1)

  expect_s3_class(result, "rerand_design_result")
  expect_equal(result$n_treat, 20)
  expect_equal(sum(result$Z), 20)
  expect_true(is.finite(result$mahalanobis))

  X_constant <- matrix(1, nrow = 40, ncol = 1)
  expect_error(
    rerand_design(X_constant, n_treat = 20, accept_prob = 1),
    "non-constant"
  )
})

test_that("estimation validates assignment and criterion inputs", {
  Y <- seq_len(10)
  data <- data.frame(Y = Y, Z = rep(0:1, each = 5), x = seq_len(10))
  expect_error(
    rerand_estimate(Y ~ Z, data),
    "accept_prob must be explicitly set to 1"
  )
  expect_error(
    rerand_estimate_matrix(Y, c(rep(1, 9), 0), accept_prob = 1),
    "enough observations"
  )
  expect_error(
    rerand_estimate(Y ~ Z, data, method = "lin"),
    "requires at least one covariate"
  )
})
