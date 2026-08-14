# rerand Refactor Plan

## Goal

This refactor stabilizes the existing Mahalanobis rerandomization design,
difference-in-means estimator, Lin estimator, theoretical variance routines,
and quantile simulation. It also makes the R and C++ engines use the same
validation and result contracts.

The work is limited to the `rerand/` repository. It does not add new
rerandomization criteria or new estimators.

## Public API

The canonical user-facing functions are:

```r
rerand_design(
  X, n_treat,
  accept_prob = NULL,
  threshold = NULL,
  max_tries = 10000L,
  seed = NULL,
  engine = c("cpp", "R")
)

rerand_estimate(
  Y_obs, Z, X = NULL,
  method = c("dim", "lin"),
  accept_prob = NULL,
  threshold = NULL,
  theoretical = FALSE,
  Y_full = NULL
)

get_quantile(
  R2, K,
  accept_prob = NULL,
  threshold = NULL,
  alpha = 0.975,
  n_sim = 100000L,
  engine = c("cpp", "R")
)
```

`accept_prob` and `threshold` are mutually exclusive. Both must be supplied
as `NULL` or exactly one must be supplied; design requires one criterion,
while estimation and quantile calculations use `accept_prob = 1` explicitly
for complete randomization.

Design results use class `rerand_design_result` and estimation results use
class `rerand_estimate_result`. Both remain ordinary lists with stable named
fields and provide `print()` and `summary()` methods.

## Implementation batches

1. `docs: add rerand refactor plan`

   Add this file and verify it with `git diff --check`.

2. `refactor: add shared validation and criterion handling`

   Add shared input validation, criterion resolution, covariance handling, and
   the initial `testthat` infrastructure.

3. `refactor: unify rerandomization design engines`

   Refactor the R and C++ design engines, add the design S3 result, and test
   assignment invariants, threshold behavior, and engine parity.

4. `refactor: unify estimators and inference`

   Refactor difference-in-means, Lin adjustment, theoretical statistics, and
   quantile simulation. Fix R/C++ distribution mismatches and add regression
   tests.

5. `refactor: add deprecated API compatibility wrappers`

   Preserve historical exports and parameter names through deprecated wrappers
   that forward to the canonical API.

6. `docs: update API documentation and migration guide`

   Update package metadata, README examples, generated documentation, and
   migration notes. Finish with `testthat` and `R CMD check` validation.

## Validation requirements

Tests cover valid and invalid inputs, treatment allocation, acceptance
thresholds, R/C++ engine invariants, singular or rank-deficient covariates,
both estimators, theoretical statistics, quantiles, deprecated wrappers, and
S3 printing/summary methods.

Each implementation batch is checked with targeted tests, `git diff --check`,
and `git status --short --branch`. Unrelated files and the sibling
`Rerand-Review/` repository remain untouched.
