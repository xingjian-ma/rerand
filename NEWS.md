# rerand 0.0.0.9000

## Usability and API overhaul

- `rerand_design()` now prepares a data-frame design and `rerand_assign()`
  generates a realized assignment with explicit design provenance.
- `rerand_estimate()` supports formula-derived DIM, ANCOVA, and Lin estimators
  or string selectors with an explicit estimator.
- `rerand_inference()` inherits CRE normal or ReM mixture inference from the
  assignment design.
- `rerand_compare()` compares named inference objects and `rerand_evaluate()`
  evaluates them against mapped two-column potential outcomes.
- Legacy matrix shortcuts, standalone quantile/population-statistics helpers,
  balance functions, and deprecated aliases are no longer public.
