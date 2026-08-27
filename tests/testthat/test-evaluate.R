make_evaluation_inputs <- function() {
  set.seed(44)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  id <- paste0("u", seq_len(n))
  design_data <- data.frame(id = id, x1 = x1, x2 = x2)
  design <- rerand_design(
    design_data, n_treat = n / 2, id = "id", accept_prob = 1
  )
  assignment <- rerand_assign(design, seed = 13)
  Y0 <- 0.4 * x1 + rnorm(n)
  Y1 <- Y0 + 1 + 0.3 * x2
  data <- data.frame(
    id = id, Y = ifelse(assignment$Z == 1, Y1, Y0), Z = assignment$Z,
    x1 = x1, x2 = x2, Y0 = Y0, Y1 = Y1
  )
  estimates <- list(
    dim = rerand_estimate(data, assignment, formula = Y ~ Z),
    lin = rerand_estimate(data, assignment, formula = Y ~ Z * (x1 + x2))
  )
  list(data = data, estimates = estimates)
}

test_that("evaluate reports finite-population error and interval coverage", {
  inputs <- make_evaluation_inputs()
  inferences <- lapply(inputs$estimates, rerand_inference)
  comparison <- rerand_compare(c(dim = inferences$dim, lin = inferences$lin))
  potential <- inputs$data[, c("Y0", "Y1")]
  evaluation <- rerand_evaluate(
    comparison, potential,
    potential_outcomes = c(control = "Y0", treated = "Y1")
  )
  expect_s3_class(evaluation, "rerand_evaluation")
  expect_equal(evaluation$true_effect, mean(inputs$data$Y1 - inputs$data$Y0))
  expect_equal(nrow(evaluation$results), 2)
  expect_true(all(c(
    "signed_error", "absolute_error", "squared_error", "covered"
  ) %in% names(evaluation$results)))
  expect_output(print(evaluation), "Rerandomization evaluation")

  single <- rerand_evaluate(
    inferences$dim, potential,
    potential_outcomes = c(control = "Y0", treated = "Y1")
  )
  expect_equal(nrow(single$results), 1)

  expect_error(
    rerand_evaluate(inputs$estimates$dim, potential,
                    potential_outcomes = c(control = "Y0", treated = "Y1")),
    "rerand_inference or rerand_compare"
  )
  expect_error(
    rerand_evaluate(inferences$dim, potential,
                    potential_outcomes = c(control = "Y0", treated = "Y0")),
    "distinct"
  )
  expect_error(
    rerand_evaluate(inferences$dim, potential[, 1, drop = FALSE],
                    potential_outcomes = c(control = "Y0", treated = "Y1")),
    "exactly two"
  )
  expect_error(
    rerand_evaluate(inferences$dim, potential[-1, ],
                    potential_outcomes = c(control = "Y0", treated = "Y1")),
    "one row"
  )
})
