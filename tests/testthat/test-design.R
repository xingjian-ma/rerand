make_design_data <- function() {
  set.seed(11)
  data.frame(matrix(rnorm(120), nrow = 40, ncol = 3,
                    dimnames = list(NULL, c("x1", "x2", "x3"))))
}

test_that("design preparation validates covariates and criteria", {
  data <- make_design_data()
  expect_error(
    rerand_design(data, n_treat = 20),
    "Exactly one of accept_prob or threshold"
  )
  expect_error(
    rerand_design(data, n_treat = 20, accept_prob = 0.5, threshold = 1),
    "mutually exclusive"
  )
  expect_error(
    rerand_design(data, n_treat = 20, accept_prob = 0),
    "accept_prob"
  )
  expect_error(
    rerand_design(data, n_treat = 20, threshold = 0),
    "threshold"
  )

  design <- rerand_design(data, n_treat = 20, accept_prob = 1)
  expect_s3_class(design, "rerand_design")
  expect_equal(design$n_treat, 20)
  expect_equal(design$design_method, "cre")
  expect_equal(unname(crossprod(design$whitened) / 39), diag(3),
               tolerance = 1e-10)
  expect_output(print(design), "Rerandomization design")
  expect_s3_class(summary(design), "summary.rerand_design")

  threshold_design <- rerand_design(data, n_treat = 20, threshold = 100)
  expect_equal(threshold_design$criterion$type, "threshold")
  expect_null(threshold_design$criterion$accept_prob)
  expect_equal(threshold_design$design_method, "rem")

  factor_data <- data.frame(
    unit = paste0("u", seq_len(20)), x = seq_len(20),
    group = rep(c("a", "b"), 10), unused = rnorm(20)
  )
  factor_design <- rerand_design(
    factor_data, n_treat = 10, formula = ~ x + group, id = "unit",
    accept_prob = 0.2
  )
  expect_equal(factor_design$unit_id, factor_data$unit)
  expect_equal(colnames(factor_design$X), c("x", "groupb"))
  expect_equal(factor_design$id_name, "unit")
  expect_identical(factor_design$data, factor_data)

  all_columns <- rerand_design(
    data.frame(
      id = seq_len(20), x = rnorm(20), flag = rep(c(TRUE, FALSE), 10),
      category = rep(c("a", "b"), 10)
    ),
    n_treat = 10, id = "id", accept_prob = 0.2
  )
  expect_true(all(c("x", "flagTRUE", "categoryb") %in% colnames(all_columns$X)))
  expect_false("id" %in% colnames(all_columns$X))

  expect_error(
    rerand_design(data.frame(id = c(1, 1, 2, 3), x = 1:4),
                  n_treat = 2, id = "id", accept_prob = 0.2),
    "unique"
  )
  expect_error(
    rerand_design(data.frame(id = 1:4, x = 1:4), n_treat = 2,
                  formula = ~ id + x, id = "id", accept_prob = 0.2),
    "id column"
  )
  expect_error(
    rerand_design(data.frame(x = 1:4,
                             when = as.Date("2020-01-01") + 0:3),
                  n_treat = 2, accept_prob = 0.2),
    "Unsupported"
  )
  expect_error(
    rerand_design(matrix(1:8, nrow = 4), n_treat = 2, accept_prob = 0.2),
    "data must be a data frame"
  )

  collinear <- rerand_design(
    data.frame(x = seq_len(20), twice_x = 2 * seq_len(20)),
    n_treat = 10, accept_prob = 0.2
  )
  expect_equal(collinear$effective_rank, 1L)
  expect_equal(ncol(collinear$whitened), 1L)
  expect_error(
    rerand_design(data.frame(x = rep(1, 20), y = rep(1, 20)),
                  n_treat = 10, accept_prob = 0.2),
    "non-constant"
  )
})
