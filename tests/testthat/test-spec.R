test_that("spec accepts matrix and data-frame inputs", {
  set.seed(21)
  X <- matrix(rnorm(120), nrow = 40, ncol = 3,
              dimnames = list(NULL, c("x1", "x2", "x3")))
  matrix_spec <- rerand_spec(X, n_treat = 20, accept_prob = 0.1)
  frame_spec <- rerand_spec(as.data.frame(X), n_treat = 20,
                            accept_prob = 0.1)

  expect_s3_class(matrix_spec, "rerand_spec")
  expect_equal(matrix_spec$X, frame_spec$X, ignore_attr = TRUE)
  expect_equal(unname(matrix_spec$whitened), unname(frame_spec$whitened),
               tolerance = 1e-10)
  expect_equal(unname(crossprod(matrix_spec$whitened) / 39), diag(3),
               tolerance = 1e-10)
  expect_output(print(matrix_spec), "Rerandomization specification")
  expect_s3_class(summary(matrix_spec), "summary.rerand_spec")
})

test_that("spec formulas select and encode covariates", {
  data <- data.frame(
    unit = paste0("u", seq_len(20)),
    x = seq_len(20),
    group = rep(c("a", "b"), 10),
    unused = rnorm(20)
  )
  spec <- rerand_spec(
    data,
    n_treat = 10,
    formula = ~ x + group,
    id = "unit",
    accept_prob = 0.2
  )

  expect_equal(spec$unit_id, data$unit)
  expect_equal(colnames(spec$X), c("x", "groupb"))
  expect_equal(ncol(spec$X), 2)
  expect_equal(spec$id_name, "unit")
  expect_identical(spec$data, data)
})

test_that("spec without a formula uses all non-ID columns", {
  data <- data.frame(
    id = seq_len(20),
    x = rnorm(20),
    flag = rep(c(TRUE, FALSE), 10),
    category = rep(c("a", "b"), 10)
  )
  spec <- rerand_spec(data, n_treat = 10, id = "id", accept_prob = 0.2)

  expect_true(all(c("x", "flagTRUE", "categoryb") %in% colnames(spec$X)))
  expect_false("id" %in% colnames(spec$X))
})

test_that("spec validates IDs, formulas, and unsupported columns", {
  data <- data.frame(id = c(1, 1, 2, 3), x = 1:4)
  expect_error(
    rerand_spec(data, n_treat = 2, id = "id", accept_prob = 0.2),
    "unique"
  )
  expect_error(
    rerand_spec(data.frame(id = 1:4, x = 1:4), n_treat = 2,
                formula = ~ id + x, id = "id", accept_prob = 0.2),
    "id column"
  )
  expect_error(
    rerand_spec(data.frame(x = 1:4, when = as.Date("2020-01-01") + 0:3),
                n_treat = 2, accept_prob = 0.2),
    "Unsupported"
  )
  expect_error(
    rerand_spec(matrix(1:8, nrow = 4), n_treat = 2,
                formula = ~ x1, accept_prob = 0.2),
    "column names"
  )
})

test_that("spec uses one rank definition for collinear covariates", {
  x <- seq_len(20)
  spec <- rerand_spec(cbind(x = x, twice_x = 2 * x), n_treat = 10,
                        accept_prob = 0.2)

  expect_equal(spec$effective_rank, 1L)
  expect_equal(ncol(spec$whitened), 1L)
  expect_error(
    rerand_spec(matrix(1, nrow = 20, ncol = 2), n_treat = 10,
                accept_prob = 0.2),
    "non-constant"
  )
})
