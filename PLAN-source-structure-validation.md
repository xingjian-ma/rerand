---
task: "source-structure-validation"
status: review
created: "2026-09-07"
updated: "2026-09-07"
---

# Source structure and validation refactor plan

## 1. Goal

### Desired outcome

Centralize public-input validation in `R/validate.R`, shared non-validation
helpers in `R/utils.R`, and keep business-specific helpers with their owning
pipeline module. Rename private helpers to concise `.` names without the
`.rerand_` prefix while preserving public behavior.

### Verifiable result

The four obsolete helper files are removed, the six pipeline modules remain,
all private helper references use the new names, and package tests/checks pass.

## 2. Context

### Current state

The package uses the workflow design -> assignment -> estimation -> inference
-> comparison/evaluation. Validation and helper code were mixed across
`data_utils.R`, `criterion.R`, `design_utils.R`, and `validation.R`.

### Key files and modules

| Area | Path | Relevance |
|---|---|---|
| Pipeline | `R/design.R`, `R/assign.R`, `R/estimate.R`, `R/inference.R`, `R/compare.R`, `R/evaluate.R` | Public APIs and owning business helpers. |
| Shared helpers | `R/utils.R` | Shared non-validation preparation, alignment, criterion, and RNG helpers. |
| Validation | `R/validate.R` | `.validate_*` user-input validators. |
| Tests | `tests/testthat/` | Regression checks and private quantile test reference. |

### Worktree and Branch

Use the current working directory; do not create a new worktree.

- **Worktree:** `D:\\xingj\\Documents\\Research\\rerand`
- **Branch:** `main` (existing user checkout)
- **Base branch:** `main`

### Dependencies and assumptions

- Public APIs, S3 registrations, classes, return fields, accepted inputs, and
  primary error/warning behavior remain unchanged.
- Existing user changes to `.gitignore` are preserved and excluded from this
  task.
- C++ code and generated Rcpp bindings are unchanged.

## 3. Constraints

### Must not change or break

- Mahalanobis rerandomization, seeded RNG restoration, R/C++ parity, ID
  alignment, CRE/ReM classification, estimators, and inference calculations.
- Public exports and S3 method names.

### Technical constraints

- Use `.validate_*` for validation helpers.
- Keep shared non-validation helpers in `utils.R` and module-specific helpers
  in their pipeline files.
- Do not retain `.rerand_*` compatibility aliases.

### Scope exclusions

- No new estimators, options, public aliases, C++ changes, or statistical
  formula changes.

### Safety and handling constraints

- Preserve existing user changes.
- Do not expose secrets or credentials.
- Do not run destructive commands without explicit approval.

## 4. Definition of Done

### Batch 1: Centralize public validation

- **Commit:** `refactor: centralize public input validation`
- **Scope:** Create `validate.R`, rename generic validators, move public-input
  validation there, and update all callers.

#### Acceptance

- [x] Public APIs preserve existing rejection and warning behavior.
- [x] Validation helpers are named `.validate_*` and live in `validate.R`.
- [x] Targeted tests pass.

#### Batch Review

- **Result:** `passed`
- **Timestamp:** `2026-09-07T20:45:00-04:00`
- **Checks:** `testthat::test_local(".")` passed for all six test files.

### Batch 2: Consolidate internal helpers

- **Commit:** `refactor: consolidate internal R helpers`
- **Scope:** Add `utils.R`, absorb shared helpers, move module-specific helpers
  into pipeline files, remove obsolete helper files, and update private test
  references.

#### Acceptance

- [x] `data_utils.R`, `criterion.R`, `design_utils.R`, and `validation.R` are removed.
- [x] No private `.rerand_*` helper definitions or calls remain.
- [x] Full testthat suite passes.

#### Batch Review

- **Result:** `passed`
- **Timestamp:** `2026-09-07T20:55:00-04:00`
- **Checks:** `testthat::test_local(".")` passed for all six test files; `git diff --check` passed.

### Final acceptance

- [x] `devtools::document(".")` completes.
- [x] `Rcpp::compileAttributes(".")` produces no unintended changes.
- [x] `R CMD build .`, `R CMD check .`, and `git diff --check` pass.
- [x] Automatic review is recorded; user review remains pending.

#### Automatic Review

- **Result:** `passed`
- **Timestamp:** `2026-09-07T21:05:00-04:00`
- **Checks:** `testthat::test_local(".")`, documentation generation, Rcpp
  attribute generation, clean-copy `R CMD build`, clean-copy `R CMD check
  --no-manual` under `LC_ALL=C`, and `git diff --check` passed. Direct build
  from the working tree was affected only by recursive `.git` metadata paths.

#### User Review

- **Result:** `pending`
- **Timestamp:** `2026-09-07T21:05:00-04:00`
- **Approved Commit SHA:** `pending`
