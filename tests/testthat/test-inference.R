make_inference_data <- function(accept_prob = 1) {
  set.seed(12)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  id <- paste0("u", seq_len(n))
  design_data <- data.frame(id = id, x1 = x1, x2 = x2)
  design <- rerand_design(
    design_data, n_treat = n / 2, id = "id", accept_prob = accept_prob
  )
  assignment <- rerand_assign(design, seed = 21)
  Y0 <- 0.5 * x1 + rnorm(n)
  Y1 <- Y0 + 1 + 0.5 * x2
  data <- data.frame(
    id = id, Y = ifelse(assignment$Z == 1, Y1, Y0),
    Z = assignment$Z, x1 = x1, x2 = x2
  )
  list(data = data, assignment = assignment)
}

test_that("CRE inference uses normal reference for every estimator", {
  inputs <- make_inference_data(accept_prob = 1)
  estimates <- list(
    dim = rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z),
    ancova = rerand_estimate(
      inputs$data, inputs$assignment, formula = Y ~ Z + x1 + x2
    ),
    lin = rerand_estimate(
      inputs$data, inputs$assignment, formula = Y ~ Z * (x1 + x2)
    )
  )

  for (estimate in estimates) {
    inference <- rerand_inference(estimate)
    expect_s3_class(inference, "rerand_inference")
    expect_equal(inference$design_method, "cre")
    expect_equal(inference$reference_distribution, "normal")
    expect_equal(inference$critical_value, stats::qnorm(0.975))
    expect_equal(confint(inference), inference$interval)
    expect_true(all(is.finite(inference$interval)))
  }
})

test_that("ReM inference uses estimator-specific ReM intervals", {
  inputs <- make_inference_data(accept_prob = 0.2)
  estimates <- list(
    dim = rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z),
    ancova = rerand_estimate(
      inputs$data, inputs$assignment, formula = Y ~ Z + x1 + x2
    ),
    lin = rerand_estimate(
      inputs$data, inputs$assignment, formula = Y ~ Z * (x1 + x2)
    )
  )

  for (estimate in estimates) {
    inference <- rerand_inference(estimate)
    expect_equal(inference$design_method, "rem")
    expect_equal(inference$reference_distribution, "rem")
    expect_true(is.finite(inference$R2_hat))
    expect_true(is.finite(inference$critical_value))
    expect_true(all(is.finite(inference$interval)))
    expect_equal(inference$standard_error, estimate$unadjusted_standard_error)
  }
})

test_that("inference validates object type and confidence level", {
  inputs <- make_inference_data()
  estimate <- rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z)
  expect_error(rerand_inference(inputs$assignment), "rerand_estimate")
  expect_error(rerand_inference(estimate, level = 1), "level")
  expect_error(confint(rerand_inference(estimate), level = 0.9), "fixed")
})
