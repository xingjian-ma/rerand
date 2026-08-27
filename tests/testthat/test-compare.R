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

test_that("compare combines compatible named inference objects", {
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

  expect_error(rerand_compare(c(inferences$dim, inferences$ancova)), "named")
  expect_error(rerand_compare(c(dim = inferences$dim)), "at least two")
  different_level <- rerand_inference(inputs$estimates$ancova, level = 0.9)
  expect_error(
    rerand_compare(c(dim = inferences$dim, ancova = different_level)),
    "same confidence level"
  )
  rem_inputs <- make_compare_inputs(accept_prob = 0.2)
  rem <- rerand_inference(rem_inputs$estimates$dim)
  expect_error(
    rerand_compare(c(cre = inferences$dim, rem = rem)),
    "provenance"
  )
})
