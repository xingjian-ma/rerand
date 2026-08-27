rerand
================

# rerand

The `rerand` package provides Mahalanobis rerandomization designs and
randomization-based treatment-effect estimators. The workflow keeps design,
assignment, estimation, inference, comparison, and evaluation as separate
objects.

## Installation

The development version can be installed from the repository with:

``` r
remotes::install_github("YOUR_USERNAME/rerand")
```

## Design and assignment

Prepare a data-frame design from pretreatment covariates, then generate a
realized assignment. A design uses either `accept_prob` or a direct
Mahalanobis `threshold`.

``` r
library(rerand)

set.seed(123)
covariates <- data.frame(
  id = paste0("unit-", seq_len(100)),
  age = rnorm(100, 45, 12),
  site = sample(c("A", "B", "C"), 100, replace = TRUE)
)

design <- rerand_design(
  covariates,
  n_treat = 50,
  formula = ~ age + site,
  id = "id",
  accept_prob = 0.1
)
assignment <- rerand_assign(design, seed = 123, engine = "cpp")
assigned_data <- as.data.frame(assignment)
```

`rerand_design()` returns a reusable `rerand_design` object. Every estimate
requires the corresponding `rerand_assignment`, including complete
randomization (`accept_prob = 1`).

## Estimation

The public estimate data must be a complete data frame containing the observed
outcome, treatment, and any adjustment covariates. Formula mode derives the
estimator:

``` r
Y0 <- 0.1 * covariates$age + rnorm(100)
Y1 <- Y0 + 1
analysis_data <- covariates
analysis_data$Z <- assignment$Z
analysis_data$Y <- ifelse(analysis_data$Z == 1, Y1, Y0)

dim_fit <- rerand_estimate(analysis_data, assignment, formula = Y ~ Z)
ancova_fit <- rerand_estimate(
  analysis_data, assignment, formula = Y ~ Z + age
)
lin_fit <- rerand_estimate(
  analysis_data, assignment, formula = Y ~ Z * (age + site)
)
```

`Y ~ Z` is difference-in-means (DIM), `Y ~ Z + x` is additive ANCOVA, and
`Y ~ Z * (x1 + x2)` is Lin's fully interacted estimator. Alternatively, use
string selectors with an explicit `estimator = "dim"`, `"ancova"`, or `"lin"`.
Formula mode takes priority and warns once if selector arguments are supplied
at the same time.

## Inference, comparison, and evaluation

Inference inherits the design method. CRE uses a normal reference distribution;
ReM uses the estimator's rerandomization mixture interval.

``` r
dim_inference <- rerand_inference(dim_fit)
lin_inference <- rerand_inference(lin_fit)
confint(dim_inference)

comparison <- rerand_compare(c(dim = dim_inference, lin = lin_inference))
```

To evaluate against finite-population truth, provide a separate two-column
potential-outcomes data frame and map its control and treated columns:

``` r
potential_data <- data.frame(Y0 = Y0, Y1 = Y1)
evaluation <- rerand_evaluate(
  comparison,
  potential_data,
  potential_outcomes = c(control = "Y0", treated = "Y1")
)
evaluation
```

The evaluation reports the true average treatment effect, signed and absolute
errors, squared error, and confidence-interval coverage for each inference.
