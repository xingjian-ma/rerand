test_that("integration is deterministic and handles limiting cases", {
  q1 <- get_quantile(0.4, 2, accept_prob = 0.2)
  q2 <- get_quantile(0.4, 2, accept_prob = 0.2)

  expect_identical(q1, q2)
  expect_equal(q1, 1.571342, tolerance = 1e-5)
  expect_equal(get_quantile(0, 2, accept_prob = 0.2), stats::qnorm(0.975))
  expect_equal(get_quantile(0.4, 2, accept_prob = 1), stats::qnorm(0.975))
  expect_true(is.finite(get_quantile(1, 1, accept_prob = 0.2)))
})

test_that("simulation is an explicit seeded fallback", {
  set.seed(22)
  before <- .Random.seed
  simulated <- get_quantile(
    0.4, 2, accept_prob = 0.2, method = "simulation", n_sim = 5000,
    seed = 9, engine = "R"
  )

  expect_identical(.Random.seed, before)
  expect_lt(abs(simulated - get_quantile(0.4, 2, accept_prob = 0.2)), 0.08)
})

test_that("difference-in-means confidence intervals use ReM quantiles", {
  set.seed(31)
  data <- data.frame(
    Y = rnorm(40), Z = rep(c(1, 0), each = 20),
    x1 = rnorm(40), x2 = rnorm(40)
  )
  result <- rerand_estimate(Y ~ Z + x1 + x2, data, accept_prob = 0.2)
  interval <- confint(result)

  expect_equal(dim(interval), c(1, 2))
  expect_true(all(is.finite(interval)))
})
