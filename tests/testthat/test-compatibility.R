test_that("historical design API forwards with a deprecation warning", {
  set.seed(14)
  X <- matrix(rnorm(80), nrow = 40, ncol = 2)
  expect_warning(
    result <- rerand.design(X, n1 = 20, p_accept = 1),
    "deprecated"
  )
  expect_true(all(c("Z", "tries", "M", "threshold", "p_accept", "accepted") %in%
                    names(result)))
  expect_equal(sum(result$Z), 20)
})

test_that("historical estimator API forwards with compatible fields", {
  set.seed(15)
  Y <- rnorm(40)
  Z <- rep(c(1, 0), each = 20)
  expect_warning(
    result <- est_dim(Y, Z, p_accept = 1),
    "deprecated"
  )
  expect_true(all(c("tau_hat", "se_neyman") %in% names(result)))
  expect_true(is.finite(result$tau_hat))
})
