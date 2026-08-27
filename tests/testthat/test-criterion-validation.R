test_that("design criterion is explicit and mutually exclusive", {
  data <- data.frame(x1 = seq_len(20), x2 = rev(seq_len(20)))
  expect_error(
    rerand_design(data, n_treat = 10),
    "Exactly one of accept_prob or threshold"
  )
  expect_error(
    rerand_design(data, n_treat = 10, accept_prob = 0.5, threshold = 1),
    "mutually exclusive"
  )
  expect_error(
    rerand_design(data, n_treat = 10, accept_prob = 0),
    "accept_prob"
  )
  expect_error(
    rerand_design(data, n_treat = 10, threshold = 0),
    "threshold"
  )
})

test_that("rank-deficient covariates are handled consistently", {
  set.seed(10)
  x <- rnorm(40)
  data <- data.frame(x1 = x, x2 = 2 * x)
  result <- rerand_design(data, n_treat = 20, accept_prob = 1)
  assignment <- rerand_assign(result)

  expect_s3_class(result, "rerand_design")
  expect_equal(result$n_treat, 20)
  expect_equal(sum(assignment$Z), 20)
  expect_true(is.finite(assignment$mahalanobis))

  constant <- data.frame(x = rep(1, 40))
  expect_error(
    rerand_design(constant, n_treat = 20, accept_prob = 1),
    "non-constant"
  )
})

test_that("estimation validates assignment and method inputs", {
  data <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  design <- rerand_design(data, n_treat = 20, accept_prob = 1)
  assignment <- rerand_assign(design, seed = 22)
  data$Z <- assignment$Z
  data$Y <- rnorm(nrow(data)) + data$Z
  expect_error(
    rerand_estimate(data, design, formula = Y ~ Z),
    "rerand_assignment"
  )
  expect_error(
    rerand_estimate(data, assignment, outcome = "Y",
                    treatment = "Z", estimator = "unknown"),
    "estimator"
  )
  expect_error(
    rerand_estimate(data, assignment, outcome = "Y",
                    treatment = "Z", estimator = "ancova"),
    "requires"
  )
})
