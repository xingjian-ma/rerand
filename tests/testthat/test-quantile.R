make_quantile_data <- function(accept_prob = 0.2) {
  set.seed(31)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  id <- paste0("u", seq_len(n))
  design <- rerand_design(
    data.frame(id = id, x1 = x1, x2 = x2), n_treat = n / 2,
    id = "id", accept_prob = accept_prob
  )
  assignment <- rerand_assign(design, seed = 9)
  data <- data.frame(
    id = id, Z = assignment$Z,
    Y = rnorm(n) + assignment$Z,
    x1 = x1, x2 = x2
  )
  list(data = data, assignment = assignment)
}

test_that("ReM quantiles are deterministic and finite", {
  inputs <- make_quantile_data()
  estimate <- rerand_estimate(inputs$data, inputs$assignment, formula = Y ~ Z)
  first <- rerand_inference(estimate)
  second <- rerand_inference(estimate)
  expect_identical(first$critical_value, second$critical_value)
  expect_true(is.finite(first$critical_value))
  expect_true(all(is.finite(confint(first))))
})

test_that("simulation quantiles preserve the caller RNG state", {
  criterion <- rerand_design(
    data.frame(x1 = rnorm(20), x2 = rnorm(20)),
    n_treat = 10, accept_prob = 0.2
  )$criterion
  set.seed(22)
  before <- .Random.seed
  simulated <- rerand:::.rerand_quantile(
    R2 = 0.4, K = criterion$K, threshold = criterion$threshold,
    alpha = 0.975, method = "simulation", n_sim = 1000,
    seed = 9, engine = "R"
  )
  expect_identical(.Random.seed, before)
  expect_true(is.finite(simulated))
})
