rerand
================

# rerand

The `rerand` package provides Mahalanobis rerandomization designs and
randomization-based treatment effect estimators. The package supports
both a readable R reference engine and a C++ engine for larger
simulations.

## Installation

The development version can be installed from the repository with:

``` r
remotes::install_github("YOUR_USERNAME/rerand")
```

## Design

The design API requires exactly one acceptance criterion. An acceptance
probability is converted to a Mahalanobis threshold using the effective
rank of the covariate matrix.

``` r
library(rerand)

set.seed(123)
X <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)

design <- rerand_design(
  X = X,
  n_treat = 50,
  accept_prob = 0.1,
  max_tries = 10000,
  engine = "cpp"
)

design
design$Z
```

The R reference implementation can be selected with `engine = "R"`. A
direct threshold can be supplied instead of `accept_prob`:

``` r
design <- rerand_design(X, n_treat = 50, threshold = 1.5)
```

## Estimation

`rerand_estimate()` supports difference-in-means and Lin’s fully
interacted regression adjustment. Use `accept_prob = 1` for complete
randomization.

``` r
Y0 <- rnorm(100)
Y1 <- Y0 + 1 + X[, 1]
Y_obs <- ifelse(design$Z == 1, Y1, Y0)

estimate <- rerand_estimate(
  Y_obs = Y_obs,
  Z = design$Z,
  X = X,
  method = "dim",
  accept_prob = 0.1
)

estimate
estimate$tau_hat
estimate$se_ding
```

Set `theoretical = TRUE` and provide `Y_full = cbind(Y0, Y1)` to
calculate population-level benchmarks alongside the sample estimate.

## Compatibility

The canonical API uses snake_case names. Historical point-notation and
low-level function names remain available temporarily with deprecation
warnings. New code should use `rerand_design()`, `rerand_estimate()`,
and `get_quantile()`.
