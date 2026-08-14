make_estimate_data <- function() {
  set.seed(12)
  X <- matrix(rnorm(120), nrow = 40, ncol = 3,
              dimnames = list(NULL, c("x1", "x2", "x3")))
  Z <- rep(c(1, 0), each = 20)
  Y0 <- 0.5 * X[, 1] + rnorm(40)
  Y1 <- Y0 + 1 + 0.5 * X[, 2]
  data.frame(
    id = paste0("u", seq_len(40)), Y = ifelse(Z == 1, Y1, Y0), Z = Z,
    X, Y0 = Y0, Y1 = Y1
  )
}

test_that("formula difference-in-means supports CRE and rerandomization", {
  data <- make_estimate_data()
  cre <- rerand_estimate(Y ~ Z, data, accept_prob = 1)
  rem <- rerand_estimate(Y ~ Z + x1 + x2, data, accept_prob = 0.2)

  expect_s3_class(cre, "rerand_estimate_result")
  expect_equal(cre$method, "dim")
  expect_true(is.finite(cre$tau_hat))
  expect_true(is.finite(cre$se_neyman))
  expect_null(cre$se_ding)
  expect_true(is.finite(rem$se_ding))
  expect_true(rem$sample_stats$R2_hat >= 0)
  expect_equal(unname(coef(rem)), rem$tau_hat)
  expect_equal(dim(vcov(rem)), c(1, 1))
  expect_equal(dim(confint(rem)), c(1, 2))
})

test_that("design objects supply the ReM assignment and criterion", {
  data <- make_estimate_data()
  design <- rerand_design(data, n_treat = 20, formula = ~ x1 + x2,
                          id = "id", accept_prob = 1, seed = 21)
  data$Z <- design$Z
  result <- rerand_estimate(Y ~ Z, data, design = design)

  expect_equal(result$accept_prob, 1)
  expect_equal(result$sample_stats$criterion_type, "probability")
  expect_warning(
    rerand_estimate(Y ~ Z + x3, data, design = design),
    "ignored"
  )
  data_reordered <- data[rev(seq_len(nrow(data))), ]
  expect_equal(
    rerand_estimate(Y ~ Z, data_reordered, design = design)$tau_hat,
    result$tau_hat
  )
})

test_that("Lin follows the formula and requires adjustment covariates", {
  data <- make_estimate_data()
  result <- rerand_estimate(Y ~ Z + x1 + x2, data, method = "lin")

  expect_s3_class(result, "rerand_estimate_result")
  expect_true(is.finite(result$tau_hat))
  expect_true(is.finite(result$se_ehw))
  expect_s3_class(result$fit, "lm")
  expect_error(
    rerand_estimate(Y ~ Z, data, method = "lin"),
    "requires at least one covariate"
  )
  expect_error(
    rerand_estimate(Y ~ Z + x1, data, method = "lin", accept_prob = 0.2),
    "does not use"
  )
  design <- rerand_design(data, n_treat = 20, formula = ~ x1,
                          id = "id", accept_prob = 1)
  data$Z <- design$Z
  expect_error(
    rerand_estimate(Y ~ Z, data, design = design, method = "lin"),
    "requires at least one covariate"
  )
})

test_that("formula treatment coding and low-level matrix estimation are explicit", {
  data <- make_estimate_data()
  data$group <- ifelse(data$Z == 1, "treated", "control")
  factor_result <- rerand_estimate(Y ~ group + x1, data,
                                   treated = "treated", accept_prob = 0.2)
  matrix_result <- rerand_estimate_matrix(
    data$Y, data$Z, as.matrix(data[, c("x1", "x2")]), accept_prob = 0.2
  )

  expect_true(is.finite(factor_result$tau_hat))
  expect_true(is.finite(matrix_result$tau_hat))
  expect_error(rerand_estimate(Y ~ Z + Z:x1, data), "interactions")
  expect_error(rerand_estimate(Y ~ group + x1, data, accept_prob = 0.2),
               "treated must")
})

test_that("population statistics are a separate calculation", {
  data <- make_estimate_data()
  design <- rerand_design(data, n_treat = 20, formula = ~ x1 + x2,
                          id = "id", accept_prob = 1)
  population <- rerand_population_stats(cbind(data$Y0, data$Y1), design = design)

  expect_equal(population$tau_true, mean(data$Y1 - data$Y0))
  expect_true(is.finite(population$se_rem_dim))
})
