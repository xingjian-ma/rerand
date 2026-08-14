make_estimate_data <- function() {
  set.seed(12)
  X <- matrix(rnorm(120), nrow = 40, ncol = 3)
  Z <- rep(c(1, 0), each = 20)
  Y0 <- 0.5 * X[, 1] + rnorm(40)
  Y1 <- Y0 + 1 + 0.5 * X[, 2]
  list(X = X, Z = Z, Y0 = Y0, Y1 = Y1, Y = ifelse(Z == 1, Y1, Y0))
}

test_that("difference-in-means supports CRE and rerandomization criteria", {
  data <- make_estimate_data()
  cre <- rerand_estimate(data$Y, data$Z, accept_prob = 1)
  rem <- rerand_estimate(data$Y, data$Z, data$X, accept_prob = 0.2)

  expect_s3_class(cre, "rerand_estimate_result")
  expect_equal(cre$method, "dim")
  expect_true(is.finite(cre$tau_hat))
  expect_true(is.finite(cre$se_neyman))
  expect_null(cre$se_ding)
  expect_true(is.finite(rem$se_ding))
  expect_true(rem$sample_stats$R2_hat >= 0)
})

test_that("Lin estimator and theoretical statistics are returned", {
  data <- make_estimate_data()
  result <- rerand_estimate(
    data$Y, data$Z, data$X, method = "lin", accept_prob = 0.2,
    theoretical = TRUE, Y_full = cbind(data$Y0, data$Y1)
  )

  expect_s3_class(result, "rerand_estimate_result")
  expect_true(is.finite(result$tau_hat))
  expect_true(is.finite(result$se_ehw))
  expect_s3_class(result$fit, "lm")
  expect_equal(result$pop_stats$tau_true, mean(data$Y1 - data$Y0))
  expect_true(is.finite(result$pop_stats$se_rem_dim))
  expect_output(print(result), "Rerandomization estimate")
  expect_s3_class(summary(result), "summary.rerand_estimate_result")
})

test_that("R and C++ quantile engines target the same distribution", {
  set.seed(13)
  q_r <- get_quantile(0.4, 2, accept_prob = 0.2, n_sim = 2000, engine = "R")
  set.seed(13)
  q_cpp <- get_quantile(0.4, 2, accept_prob = 0.2, n_sim = 2000, engine = "cpp")

  expect_true(is.finite(q_r))
  expect_true(is.finite(q_cpp))
  expect_lt(abs(q_r - q_cpp), 0.2)
})
