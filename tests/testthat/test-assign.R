make_assignment_design <- function(threshold = NULL, accept_prob = 1) {
  set.seed(11)
  data <- data.frame(matrix(rnorm(120), nrow = 40, ncol = 3,
                            dimnames = list(NULL, c("x1", "x2", "x3"))))
  rerand_design(
    data, n_treat = 20, accept_prob = accept_prob, threshold = threshold
  )
}

test_that("assignment generation preserves counts, diagnostics, and RNG", {
  design <- make_assignment_design()
  cpp <- rerand_assign(design, engine = "cpp")
  r_engine <- rerand_assign(design, engine = "R")

  expect_s3_class(cpp, "rerand_assignment")
  expect_s3_class(r_engine, "rerand_assignment")
  expect_equal(length(cpp$Z), design$n)
  expect_equal(sum(cpp$Z), design$n_treat)
  expect_equal(sum(r_engine$Z), design$n_treat)
  expect_true(cpp$accepted)
  expect_true(r_engine$accepted)
  expect_lte(cpp$mahalanobis, cpp$threshold)
  expect_lte(r_engine$mahalanobis, r_engine$threshold)
  expect_equal(cpp$criterion_type, "probability")
  expect_equal(cpp$design_method, "cre")
  expect_equal(cpp$design, design)

  threshold_design <- make_assignment_design(threshold = 100, accept_prob = NULL)
  threshold_assignment <- rerand_assign(threshold_design)
  expect_equal(threshold_assignment$criterion_type, "threshold")
  expect_null(threshold_assignment$accept_prob)
  expect_equal(threshold_assignment$threshold, 100)
  expect_equal(threshold_assignment$design_method, "rem")
  expect_true(threshold_assignment$accepted)

  strict_design <- make_assignment_design(threshold = 1e-15, accept_prob = NULL)
  expect_error(
    rerand_assign(strict_design, max_tries = 3),
    "Maximum tries"
  )
  expect_warning(
    exhausted <- rerand_assign(strict_design, max_tries = 3,
                               on_failure = "warn"),
    "Maximum tries"
  )
  expect_false(exhausted$accepted)
  expect_equal(exhausted$tries, 3)
  expect_equal(sum(exhausted$Z), design$n_treat)

  pooled <- rerand_assign(design, n_draws = 3, seed = 42)
  expect_equal(dim(pooled$pool$assignments), c(design$n, 3))
  expect_true(all(colSums(pooled$pool$assignments) == design$n_treat))
  expect_equal(nrow(pooled$pool$diagnostics), 3)
  expect_equal(nrow(as.data.frame(pooled)), design$n)
  expect_true("Z" %in% names(as.data.frame(pooled)))
  expect_output(print(pooled), "Rerandomization assignment")
  expect_s3_class(summary(pooled), "summary.rerand_assignment")

  set.seed(91)
  before <- .Random.seed
  rerand_assign(design, seed = 11)
  expect_identical(.Random.seed, before)
})
