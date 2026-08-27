---
task: "rerand-api-inference-evaluation"
status: planned
created: "2026-08-28"
updated: "2026-08-28"
---

# Rerand API, inference, comparison, and evaluation plan

## 1. Goal

### Desired outcome

Refactor the `rerand` package into an explicit workflow:

```r
design <- rerand_design(data, ...)
assignment <- rerand_assign(design, ...)
fit <- rerand_estimate(data, assignment, ...)
inference <- rerand_inference(fit)
comparison <- rerand_compare(c(dim = inference, lin = inference))
evaluation <- rerand_evaluate(
  comparison,
  potential_data,
  potential_outcomes = c(control = "Y0", treated = "Y1")
)
```

### Verifiable result

The package exposes the new public API, supports formula-derived DIM/ANCOVA/Lin and selector-based estimation, makes inference inherit the assignment design method, compares inference objects, and evaluates them against an ATE computed from a separate two-column potential-outcomes data frame.

## 2. Context

### Background

The current API mixes design specification, assignment generation, estimation, confidence intervals, quantiles, and population statistics. The refactor separates these responsibilities and makes design provenance determine the inference reference distribution.

### Current state

- Repository: `rerand/` on branch `major-package-refactor`.
- Existing exports include `rerand_spec()`, `rerand_draw()`, single-step `rerand_design()`, `rerand_estimate_matrix()`, `rerand_population_stats()`, `get_quantile()`, and `balance()`.
- The original worktree contains an untracked user file `STATE.md`; preserve it without editing or staging.

### Key files and modules

| Area | Path | Relevance |
|---|---|---|
| Design and assignment | `R/spec.R`, `R/design.R` | Split design configuration from realized assignment. |
| Estimation and inference | `R/estimate.R`, `R/estimate_utils.R` | Add estimator grammar and design-aware inference. |
| Data and validation | `R/data_utils.R`, `R/estimate_data.R` | Enforce complete data-frame and selector contracts. |
| Public API and tests | `NAMESPACE`, `tests/testthat/` | Replace exports and verify the breaking interface. |

### Worktree and Branch

- **Worktree:** `/Users/m/Documents/Research/Rerand-Workspace/.worktrees/rerand-api-inference-evaluation`
- **Branch:** `codex/rerand-api-inference-evaluation`
- **Base branch:** `major-package-refactor`

### Dependencies and assumptions

- Public data-source APIs require complete `data.frame` inputs; `rerand_evaluate()` additionally accepts a separate two-column potential-outcomes data frame.
- Every estimate requires a `rerand_assignment`, including complete randomization.
- A design is classified as `cre` when all assignments are accepted and `rem` otherwise.
- Formula mode takes priority and warns once if selector arguments or `estimator` are also supplied.
- No compatibility wrappers or deprecated aliases are retained.

## 3. Constraints

### Must not change or break

- Preserve Mahalanobis rerandomization, criterion validation, seeded RNG behavior, and R/C++ assignment parity.
- Preserve correct unit/ID alignment between design, assignment, estimation data, and comparison inputs.
- Keep all new or edited repository content in English.
- Do not modify `Rerand-Review/` or the original worktree's `STATE.md`.

### Technical constraints

- `rerand_design()` returns `rerand_design`; `rerand_assign()` returns `rerand_assignment`.
- Formula syntax derives the estimator: `Y ~ Z` is DIM, `Y ~ Z + x` is ANCOVA, and `Y ~ Z * (x1 + x2)` is Lin with internal centering and interactions.
- Selector mode uses string column names (`outcome`, `treatment`, `covariates`) and explicit `estimator`.
- `rerand_inference()` has no user-selectable CI method: CRE assignments use normal intervals, and ReM assignments use estimator-specific ReM intervals.
- `rerand_evaluate()` accepts a single inference or a comparison result, plus a two-column finite numeric potential-outcomes data frame and an explicit control/treated column mapping.
- Quantile, balance, and population-statistics helpers remain internal; they are not public `rerand_` functions.

### Scope exclusions

- Do not add estimators beyond DIM, ANCOVA, and Lin.
- Do not add automatic estimator selection or recommendation.
- Do not preserve matrix/vector public shortcuts or old function aliases.

### Safety and handling constraints

- Preserve existing user changes.
- Do not expose secrets or credentials.
- Do not run destructive commands without explicit approval.

## 4. Definition of Done

The task is complete when the Commit Flow below has been completed successfully.

### First Commit

- Commit: `docs: add PLAN-rerand-api-inference-evaluation.md`
- Include only `PLAN-rerand-api-inference-evaluation.md`.

### Batch 1: Split design and assignment APIs

- **Commit:** `refactor: split design and assignment APIs`
- **Scope:** Replace spec/draw/single-step design APIs with `rerand_design()` and `rerand_assign()`; add design method metadata and new result classes; remove matrix public paths and `balance` export.

#### Acceptance

- [ ] Design and assignment preserve IDs, treatment counts, diagnostics, seeded RNG behavior, and R/C++ parity.
- [ ] CRE and ReM classifications are correct for supported criteria.
- [ ] Legacy design, draw, spec, and balance exports are absent.

#### Batch Review

- **Result:** `passed`
- **Timestamp:** `2026-08-28T00:11:17+08:00`

### Batch 2: Add unified estimation modes

- **Commit:** `feat: add formula and selector estimators`
- **Scope:** Implement formula-derived DIM/ANCOVA/Lin, selector-based explicit estimation, mandatory assignment provenance, and canonical formula validation.

#### Acceptance

- [ ] Formula grammar selects the intended estimator and rejects noncanonical forms.
- [ ] Selector and formula modes agree for equivalent analyses.
- [ ] Treatment, assignment, ID, missing-value, and covariate validation failures are clear.

#### Batch Review

- **Result:** `passed`
- **Timestamp:** `2026-08-28T00:28:00+08:00`

### Batch 3: Add design-aware inference

- **Commit:** `feat: add design-aware inference`
- **Scope:** Add `rerand_inference()`, internalize quantile calculations, and move confidence intervals to inference objects.

#### Acceptance

- [ ] CRE assignments produce normal-reference intervals for DIM, ANCOVA, and Lin.
- [ ] ReM assignments produce estimator-specific ReM-reference intervals for DIM, ANCOVA, and Lin.
- [ ] Inference records inherited design method, reference distribution, critical value, standard error, and CI.
- [ ] Repeated accepted-assignment checks validate ReM critical values and coverage behavior against the implemented theory.

#### Batch Review

- **Result:** `passed`
- **Timestamp:** `2026-08-28T00:42:00+08:00`

### Batch 4: Add comparison and evaluation

- **Commit:** `feat: compare and evaluate inference results`
- **Scope:** Implement `c.rerand_inference()`, `rerand_compare()`, internal `.population_stats()`, and `rerand_evaluate()`.

#### Acceptance

- [ ] Comparison accepts only named, same-level, provenance-compatible inference objects.
- [ ] Evaluation accepts one inference or a comparison, not estimates or raw assignments.
- [ ] Evaluation derives `mean(Y1 - Y0)` from the explicitly mapped potential-outcome columns.
- [ ] Evaluation reports true effect, signed error, absolute error, squared error, and CI coverage.
- [ ] Invalid mappings, non-finite values, incompatible unit counts, and raw estimate inputs fail clearly.

#### Batch Review

- **Result:** `passed`
- **Timestamp:** `2026-08-28T01:03:00+08:00`

### Batch 5: Document migration and validate package

- **Commit:** `docs: document rerand inference workflow`
- **Scope:** Update roxygen, generated manuals, README, NEWS, and migration guidance; remove obsolete documentation.

#### Acceptance

- [ ] Documentation uses only the new design-aware pipeline and explains both estimation modes.
- [ ] Targeted tests and full `testthat` pass.
- [ ] `R CMD check` and `git diff --check` pass.

#### Batch Review

- **Result:** `pending`
- **Timestamp:** to be recorded after passing acceptance checks.

### Final Commit

#### Total Acceptance

- [ ] All batches are complete, reviewed, and recorded.
- [ ] Full tests, `R CMD check`, and documentation checks pass.
- [ ] Automatic review is completed and recorded.
- [ ] User approval includes the latest functional commit SHA.

#### Automatic Review

- **Result:** `pending`
- **Timestamp:** recorded after final checks.

#### User Review

- **Result:** `pending`
- **Timestamp:** recorded after manual approval.
- **Approved Commit SHA:** pending
