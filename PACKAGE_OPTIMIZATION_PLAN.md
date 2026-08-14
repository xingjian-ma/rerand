# rerand Package Optimization Plan

## Goal

This refactor makes the package use a data-and-formula interface for ordinary
workflows while retaining dedicated matrix functions for simulations and
programmatic use. It also introduces reusable design specifications, batch
assignment generation, deterministic rerandomization quantiles, and standard
S3 inference methods.

The work is limited to the `rerand/` repository. It does not add stratified,
clustered, PCA, ridge, sequential, parallel, or randomization-test methods.

## Design API

The primary design interface is:

```r
rerand_design(
  data, n_treat,
  formula = NULL,
  accept_prob = NULL,
  threshold = NULL,
  id = NULL,
  max_tries = NULL,
  failure_prob = 1e-6,
  seed = NULL,
  engine = c("cpp", "R"),
  on_failure = c("error", "warn"),
  tol = 1e-10
)
```

`data` accepts a numeric matrix or data frame. A one-sided `formula` selects
and constructs balance covariates. With no formula, all matrix columns or all
data-frame columns except `id` are used. Data-frame factors and character
columns are encoded with `model.matrix()`; unsupported columns and missing or
non-finite model-matrix values are rejected.

`rerand_spec()` performs data preparation, effective-rank determination,
whitening, and criterion resolution without drawing an assignment.
`rerand_draw()` reuses a specification to generate one or more independent
accepted assignments. `rerand_design()` is the single-draw convenience entry
point.

When `max_tries` is `NULL`, it is derived from the criterion and
`failure_prob`. Failure to reach the criterion is an error by default;
`on_failure = "warn"` retains the final rejected assignment with
`accepted = FALSE`.

Design results retain their input data and specification. `as.data.frame()`
adds the assignment to the saved data, and `balance()` returns covariate-level
balance diagnostics.

## Estimation API

The primary estimation interface is:

```r
rerand_estimate(
  formula, data,
  design = NULL,
  method = c("dim", "lin"),
  treated = NULL,
  accept_prob = NULL,
  threshold = NULL,
  se_type = NULL
)
```

In a formula such as `Y ~ Z + x1 + x2`, the response is the observed outcome,
the first right-hand-side term is the treatment, and all remaining terms are
the estimation covariates. Treatment interactions supplied by the user are
rejected because the Lin estimator constructs centered, fully interacted terms
internally.

For difference-in-means, a supplied design provides the assignment, balance
covariates, and criterion. Additional formula covariates are ignored with a
warning. Without a design, formula covariates define the rerandomization
covariates and require an explicit criterion; `Y ~ Z` requires
`accept_prob = 1`.

For the Lin estimator, the formula must contain at least one non-treatment
covariate even when a design is supplied. Only formula covariates are used.
The estimator uses HC2 inference and does not apply a rerandomization mixture
correction.

`rerand_estimate_matrix()` provides the dedicated numeric-vector and matrix
interface. `rerand_population_stats()` owns population-level calculations that
previously used `theoretical` and `Y_full` in the main estimator.

Estimate results provide `coef()`, `vcov()`, `confint()`, `as.data.frame()`,
`print()`, and `summary()` methods.

## Quantiles and Computation

`get_quantile()` uses deterministic one-dimensional integration by default for
canonical Mahalanobis rerandomization. The existing R and C++ simulation
methods remain available through `method = "simulation"`.

The design engine uses one SVD to determine effective rank and construct a
whitened covariate matrix. Both R and C++ engines compute balance in the
whitened space. The C++ engine reuses total covariate sums and supports batch
assignment generation. Initial implementation remains single-threaded.

An explicit seed is local to the function call: the previous global RNG state
is restored on exit. Calls without a seed advance the global RNG normally.

## Compatibility

This is a breaking development refactor. The design argument is renamed from
`X` to `data`; vector-based estimation moves to `rerand_estimate_matrix()`;
population quantities move to `rerand_population_stats()`; deprecated wrappers
and `p_accept` aliases are removed.

## Commit Sequence

1. `Add package optimization plan`
2. `Add formula-based rerandomization specifications`
3. `Optimize rerandomization assignment generation`
4. `Redesign treatment effect estimation interfaces`
5. `Add deterministic rerandomization inference`
6. `Remove deprecated APIs and document the new workflow`

Each implementation commit receives targeted tests and `git diff --check`.
The final commit is validated with the full test suite and `R CMD check`.
