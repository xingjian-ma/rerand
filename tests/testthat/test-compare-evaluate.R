make_compare_inputs <- function(accept_prob = 1) {
  set.seed(44)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  id <- paste0("u", seq_len(n))
  design_data <- data.frame(id = id, x1 = x1, x2 = x2)
  design <- rerand_design(
    design_data, n_treat = n / 2, id = "id", accept_prob = accept_prob
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
    ancova = rerand_estimate(data, assignment, formula = Y ~ Z + x1 + x2),
    lin = rerand_estimate(data, assignment, formula = Y ~ Z * (x1 + x2))
  )
  list(data = data, assignment = assignment, estimates = estimates)
}

test_that("named inference objects can be compared", {
  inputs <- make_compare_inputs()
  inferences <- lapply(inputs$estimates, rerand_inference)
  collection <- c(
    dim = inferences$dim, ancova = inferences$ancova, lin = inferences$lin
  )
  expect_s3_class(collection, "rerand_inference_collection")
  comparison <- rerand_compare(collection)
  expect_s3_class(comparison, "rerand_compare")
  expect_equal(comparison$results$method, c("dim", "ancova", "lin"))
  expect_equal(comparison$results$design_method, rep("cre", 3))
  expect_equal(nrow(comparison$results), 3)
  expect_output(print(comparison), "Rerandomization inference comparison")
})

test_that("comparison validates names, levels, and provenance", {
  inputs <- make_compare_inputs()
  dim <- rerand_inference(inputs$estimates$dim)
  ancova <- rerand_inference(inputs$estimates$ancova)
  expect_error(rerand_compare(c(dim, ancova)), "named")
  expect_error(rerand_compare(c(dim = dim)), "at least two")

  different_level <- rerand_inference(inputs$estimates$ancova, level = 0.9)
  expect_error(
    rerand_compare(c(dim = dim, ancova = different_level)),
    "same confidence level"
  )
  rem_inputs <- make_compare_inputs(accept_prob = 0.2)
  rem <- rerand_inference(rem_inputs$estimates$dim)
  expect_error(
    rerand_compare(c(cre = dim, rem = rem)),
    "provenance"
  )
})

test_that("evaluation uses a separate mapped two-column potential-outcome frame", {
  inputs <- make_compare_inputs()
  inferences <- lapply(inputs$estimates[c("dim", "lin")], rerand_inference)
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
})

test_that("evaluation rejects raw estimates, bad mappings, and bad row counts", {
  inputs <- make_compare_inputs()
  potential <- inputs$data[, c("Y0", "Y1")]
  inference <- rerand_inference(inputs$estimates$dim)
  expect_error(
    rerand_evaluate(inputs$estimates$dim, potential,
                    potential_outcomes = c(control = "Y0", treated = "Y1")),
    "rerand_inference or rerand_compare"
  )
  expect_error(
    rerand_evaluate(inference, potential,
                    potential_outcomes = c(control = "Y0", treated = "Y0")),
    "distinct"
  )
  expect_error(
    rerand_evaluate(inference, potential[, 1, drop = FALSE],
                    potential_outcomes = c(control = "Y0", treated = "Y1")),
    "exactly two"
  )
  expect_error(
    rerand_evaluate(inference, potential[-1, ],
                    potential_outcomes = c(control = "Y0", treated = "Y1")),
    "one row"
  )
})
