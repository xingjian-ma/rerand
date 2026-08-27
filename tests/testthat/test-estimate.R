make_estimate_data <- function(accept_prob = 1) {
  set.seed(12)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
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
    Z = assignment$Z, x1 = x1, x2 = x2, x3 = x3,
    Y0 = Y0, Y1 = Y1
  )
  list(data = data, design = design, assignment = assignment)
}

test_that("formula mode selects DIM, ANCOVA, and Lin", {
  inputs <- make_estimate_data()
  dim <- rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z)
  ancova <- rerand_estimate(
    inputs$data, inputs$assignment, formula = Y ~ Z + x1 + x2
  )
  lin <- rerand_estimate(
    inputs$data, inputs$assignment, formula = Y ~ Z * (x1 + x2)
  )

  expect_s3_class(dim, "rerand_estimate")
  expect_equal(dim$estimator, "dim")
  expect_true(is.finite(dim$estimate))
  expect_true(is.finite(dim$standard_error))
  expect_equal(dim$se_type, "neyman")
  expect_equal(ancova$estimator, "ancova")
  expect_s3_class(ancova$fit, "lm")
  expect_equal(ancova$se_type, "hc2")
  expect_equal(lin$estimator, "lin")
  expect_s3_class(lin$fit, "lm")
  expect_true(is.finite(lin$standard_error))
  expect_equal(unname(coef(lin)), lin$estimate)
  expect_equal(dim(vcov(lin)), c(1, 1))
})

test_that("selector mode requires explicit estimator and supports strings", {
  inputs <- make_estimate_data()
  result <- rerand_estimate(
    inputs$data, inputs$assignment, outcome = "Y", treatment = "Z",
    covariates = c("x1", "x2"), estimator = "lin"
  )
  expect_equal(result$estimator, "lin")
  expect_equal(result$outcome_name, "Y")
  expect_equal(result$treatment_name, "Z")
  expect_equal(result$analysis_covariates, c("x1", "x2"))

  expect_warning(
    dim <- rerand_estimate(
      inputs$data, inputs$assignment, outcome = "Y", treatment = "Z",
      covariates = "x1", estimator = "dim"
    ),
    "covariates are ignored"
  )
  expect_equal(dim$estimator, "dim")
  expect_error(
    rerand_estimate(inputs$data, inputs$assignment, outcome = "Y",
                    treatment = "Z"),
    "estimator"
  )
})

test_that("formula mode takes priority over selector arguments", {
  inputs <- make_estimate_data()
  expect_warning(
    result <- rerand_estimate(
      inputs$data, inputs$assignment, formula = Y ~ Z,
      outcome = "Y", treatment = "Z", estimator = "lin", covariates = "x1"
    ),
    "formula takes priority"
  )
  expect_equal(result$estimator, "dim")
})

test_that("assignment and data IDs are aligned explicitly", {
  inputs <- make_estimate_data()
  data_reordered <- inputs$data[rev(seq_len(nrow(inputs$data))), ]
  result <- rerand_estimate(
    inputs$data, inputs$assignment, formula = Y ~ Z
  )
  reordered <- rerand_estimate(
    data_reordered, inputs$assignment, formula = Y ~ Z
  )
  expect_equal(reordered$estimate, result$estimate)
  expect_error(
    rerand_estimate(inputs$data, inputs$design, formula = Y ~ Z),
    "rerand_assignment"
  )
  mismatched <- inputs$data
  mismatched$Z <- rev(mismatched$Z)
  expect_error(
    rerand_estimate(mismatched, inputs$assignment, formula = Y ~ Z),
    "does not match the assignment"
  )
})

test_that("formula grammar and treatment coding are validated", {
  inputs <- make_estimate_data()
  expect_error(
    rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z:x1),
    "formula must be"
  )
  expect_error(
    rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z * x1),
    "Lin formulas"
  )
  expect_error(
    rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z + Z:x1),
    "Treatment terms"
  )
  factor_data <- inputs$data
  factor_data$group <- ifelse(factor_data$Z == 1, "treated", "control")
  factor_data$Z <- NULL
  expect_true(is.finite(rerand_estimate(
    factor_data, inputs$assignment, formula = Y ~ group, treated = "treated"
  )$estimate))
  expect_error(
    rerand_estimate(factor_data, inputs$assignment, formula = Y ~ group),
    "treated must"
  )
})
