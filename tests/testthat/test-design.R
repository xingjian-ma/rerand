make_design_data <- function() {
  set.seed(11)
  matrix(rnorm(120), nrow = 40, ncol = 3)
}

test_that("both design engines satisfy the assignment contract", {
  X <- make_design_data()
  cpp <- rerand_design(X, n_treat = 20, accept_prob = 1, engine = "cpp")
  r_engine <- rerand_design(X, n_treat = 20, accept_prob = 1, engine = "R")

  expect_s3_class(cpp, "rerand_design_result")
  expect_s3_class(r_engine, "rerand_design_result")
  expect_equal(length(cpp$Z), nrow(X))
  expect_equal(sum(cpp$Z), 20)
  expect_equal(sum(r_engine$Z), 20)
  expect_true(cpp$accepted)
  expect_true(r_engine$accepted)
  expect_lte(cpp$mahalanobis, cpp$threshold)
  expect_lte(r_engine$mahalanobis, r_engine$threshold)
  expect_equal(cpp$criterion_type, "probability")
})

test_that("explicit thresholds are recorded without an inferred probability", {
  X <- make_design_data()
  result <- rerand_design(X, n_treat = 20, threshold = 100)

  expect_equal(result$criterion_type, "threshold")
  expect_null(result$accept_prob)
  expect_equal(result$threshold, 100)
  expect_true(result$accepted)
})

test_that("design reports exhausted attempts according to on_failure", {
  X <- make_design_data()
  expect_error(
    rerand_design(X, n_treat = 20, threshold = 1e-15, max_tries = 3),
    "Maximum tries"
  )
  expect_warning(
    result <- rerand_design(X, n_treat = 20, threshold = 1e-15, max_tries = 3,
                            on_failure = "warn"),
    "Maximum tries"
  )
  expect_false(result$accepted)
  expect_equal(result$tries, 3)
  expect_equal(sum(result$Z), 20)
})

test_that("design result has print and summary methods", {
  X <- make_design_data()
  result <- rerand_design(X, n_treat = 20, accept_prob = 1)
  expect_output(print(result), "Rerandomization design")
  expect_s3_class(summary(result), "summary.rerand_design_result")
})

test_that("draw pools, balance summaries, and data extraction are available", {
  X <- make_design_data()
  spec <- rerand_spec(X, n_treat = 20, accept_prob = 1)
  result <- rerand_draw(spec, n_draws = 3, seed = 42)

  expect_equal(dim(result$pool$assignments), c(40, 3))
  expect_true(all(colSums(result$pool$assignments) == 20))
  expect_equal(nrow(result$pool$diagnostics), 3)
  expect_equal(nrow(balance(result)), 3)
  expect_equal(nrow(as.data.frame(result)), 40)
  expect_true("Z" %in% names(as.data.frame(result)))
})

test_that("a seeded draw does not change the caller RNG state", {
  X <- make_design_data()
  spec <- rerand_spec(X, n_treat = 20, accept_prob = 1)
  set.seed(91)
  before <- .Random.seed
  rerand_draw(spec, seed = 11)
  expect_identical(.Random.seed, before)
})
