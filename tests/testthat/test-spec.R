test_that("design accepts data-frame inputs", {
  set.seed(21)
  X <- matrix(rnorm(120), nrow = 40, ncol = 3,
              dimnames = list(NULL, c("x1", "x2", "x3")))
  frame_design <- rerand_design(as.data.frame(X), n_treat = 20,
                                accept_prob = 0.1)

  expect_s3_class(frame_design, "rerand_design")
  expect_equal(unname(crossprod(frame_design$whitened) / 39), diag(3),
               tolerance = 1e-10)
  expect_output(print(frame_design), "Rerandomization design")
  expect_s3_class(summary(frame_design), "summary.rerand_design")
  expect_error(rerand_design(X, n_treat = 20, accept_prob = 0.1),
               "data must be a data frame")
})

test_that("spec formulas select and encode covariates", {
  data <- data.frame(
    unit = paste0("u", seq_len(20)),
    x = seq_len(20),
    group = rep(c("a", "b"), 10),
    unused = rnorm(20)
  )
  design <- rerand_design(
    data,
    n_treat = 10,
    formula = ~ x + group,
    id = "unit",
    accept_prob = 0.2
  )

  expect_equal(design$unit_id, data$unit)
  expect_equal(colnames(design$X), c("x", "groupb"))
  expect_equal(ncol(design$X), 2)
  expect_equal(design$id_name, "unit")
  expect_identical(design$data, data)
})

test_that("spec without a formula uses all non-ID columns", {
  data <- data.frame(
    id = seq_len(20),
    x = rnorm(20),
    flag = rep(c(TRUE, FALSE), 10),
    category = rep(c("a", "b"), 10)
  )
  design <- rerand_design(data, n_treat = 10, id = "id", accept_prob = 0.2)

  expect_true(all(c("x", "flagTRUE", "categoryb") %in% colnames(design$X)))
  expect_false("id" %in% colnames(design$X))
})

test_that("spec validates IDs, formulas, and unsupported columns", {
  data <- data.frame(id = c(1, 1, 2, 3), x = 1:4)
  expect_error(
    rerand_design(data, n_treat = 2, id = "id", accept_prob = 0.2),
    "unique"
  )
  expect_error(
    rerand_design(data.frame(id = 1:4, x = 1:4), n_treat = 2,
                formula = ~ id + x, id = "id", accept_prob = 0.2),
    "id column"
  )
  expect_error(
    rerand_design(data.frame(x = 1:4, when = as.Date("2020-01-01") + 0:3),
                n_treat = 2, accept_prob = 0.2),
    "Unsupported"
  )
  expect_error(
    rerand_design(matrix(1:8, nrow = 4), n_treat = 2,
                formula = ~ x1, accept_prob = 0.2),
    "data must be a data frame"
  )
})

test_that("spec uses one rank definition for collinear covariates", {
  x <- seq_len(20)
  design <- rerand_design(data.frame(x = x, twice_x = 2 * x),
                          n_treat = 10, accept_prob = 0.2)

  expect_equal(design$effective_rank, 1L)
  expect_equal(ncol(design$whitened), 1L)
  expect_error(
    rerand_design(data.frame(x = rep(1, 20), y = rep(1, 20)), n_treat = 10,
                accept_prob = 0.2),
    "non-constant"
  )
})
