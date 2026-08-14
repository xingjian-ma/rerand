# rerand 0.0.0.9000

## Usability and API overhaul

- `rerand_design()` now accepts matrices or data frames, supports one-sided
  balance formulas, ID columns, automatic factor encoding, reusable
  specifications, and diagnostic assignment pools.
- `rerand_estimate()` now uses `formula` plus `data`; its matrix interface is
  available separately as `rerand_estimate_matrix()`.
- `rerand_population_stats()` separates population benchmarks from observed
  outcome estimation.
- `get_quantile()` now uses deterministic numerical integration by default.
- Historical point-notation aliases, deprecated wrappers, and `p_accept` have
  been removed.
