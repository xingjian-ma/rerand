# rerand 0.0.0.9000

## API refactor

- Added `rerand_design()`, `rerand_estimate()`, and `get_quantile()` as the
  canonical APIs.
- Added `rerand_design_result` and `rerand_estimate_result` S3 result classes.
- Unified R and C++ criterion handling and corrected quantile simulation
  consistency.
- Retained historical exported functions as deprecated compatibility wrappers.
