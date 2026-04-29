# Slim fitted objects via `detach_*` arguments

## Why

`mfsusie()` returns a model object that retains buffers used
internally by the IBSS loop but unused by the post-fit public
API. At realistic dimensions (n=84, p=3500, M=6, T=128, L=10):

| field | shape | size |
|---|---|---:|
| `model$mu[[l]][[m]]` | p × T_basis[m] | L * sum_m(p*T_m) doubles = ~215 MB |
| `model$mu2[[l]][[m]]` | p × T_basis[m] | same as `mu` |
| `model$residuals[[m]]` | p × T_basis[m] | sum_m(p*T_m) = ~21 MB |
| `model$raw_residuals[[m]]` | n × T_basis[m] | small (per-outcome n×T) |
| `model$fitted_without_l[[m]]` | n × T_basis[m] | same shape |
| `model$fitted[[m]]` | n × T_basis[m] | small but USED by predict |
| `model$alpha` | L × p | small |
| `model$lbf_variable` | L × p | small |
| `model$lbf_variable_outcome` | L × p × M | medium |

Total fit object can exceed 250 MB on this fixture.

Audit of post-fit consumers:

| field | needed by |
|---|---|
| `mu` | `coef`, `predict`, `mf_post_smooth` (TI/HMM/scalewise/smash) |
| `mu2` | `mf_post_smooth` (HMM, scalewise — for posterior sd) |
| `alpha` | `coef`, `summary`, CS extraction, plotting |
| `sigma2` | `mf_post_smooth` (smash, smashr) |
| `fitted` | `predict` (intercept + baseline) |
| `pip`, `sets` | summary, plot, `coef` |
| `dwt_meta` | `coef`, `predict`, all smoothers |
| `residuals`, `raw_residuals`, `fitted_without_l` | NONE — IBSS-internal buffers |
| `lbf_variable`, `lbf` | summary diagnostics, `susie_post_outcome_configuration` |

**Working buffers (`residuals`, `raw_residuals`, `fitted_without_l`)
are unused after the fit returns.** They linger because the IBSS
loop never deletes them.

`mu2` is the largest field in the fit (matches `mu` in size).
It is consumed by `mf_post_smooth` only for posterior-variance
estimation in the HMM and scalewise paths. The TI smoother does
not read `mu2`. The smash smoother does not read `mu2` (it
re-derives variance from residuals). For HMM/scalewise, the
relevant quantity is the alpha-weighted aggregate
`sum_j alpha[l, j] * mu2[l, j, t]` (length T_basis), not the
full `(p × T)` matrix. We can store the aggregate (1/p smaller)
when the user opts in.

## What changes

Three new arguments to `mfsusie()` (and `fsusie()` via the
shared workhorse), each defaulting to a value that preserves
existing semantics so existing scripts keep working.

### `detach_residuals = FALSE` (default)

When `TRUE`, the fit removes:
* `model$residuals`
* `model$raw_residuals`
* `model$fitted_without_l`

after the IBSS loop converges. `model$fitted` is kept (used by
`predict.mfsusie`).

Default is `FALSE` to preserve backward compatibility, but the
roxygen recommends setting it to `TRUE` for production fits
where the fit object will be persisted (saveRDS / large pipeline
storage).

### `detach_mu2 = FALSE` (default)

When `TRUE`, the fit replaces each `model$mu2[[l]][[m]]`
(p × T_basis[m] matrix) with the alpha-weighted aggregate
`sum_j alpha[l, j] * mu2[l, j, t]` (length-T_basis vector).
Stored as `model$mu2_agg[[l]][[m]]`. The original p × T matrix
is dropped.

`mf_post_smooth` paths that consume `mu2` (HMM, scalewise) check
for the aggregate first, fall back to recomputing from `mu2` if
present, and error with a clear message if neither is present
(shouldn't happen in normal flow).

### `detach_per_effect_diagnostics = FALSE` (default)

When `TRUE`, drops `model$lbf` (length-L vector) and
`model$lbf_variable` (L × p matrix). These are diagnostic
fields used by `susie_post_outcome_configuration` and `summary`.
Setting to `TRUE` is a small win (a few MB at p=3500) but
disables `susie_post_outcome_configuration` and reduces
`summary` output.

## Acceptance criteria

* Default fit object size is unchanged from current behavior
  (all three `detach_*` default `FALSE`).
* `detach_residuals = TRUE`: fit object size drops by the size
  of `residuals` + `raw_residuals` + `fitted_without_l`. On the
  reference fixture, ≥ 20 MB reduction.
* `detach_mu2 = TRUE`: fit object size drops by approximately
  half (mu2 was ~half of mu+mu2 footprint). On the reference
  fixture, ≥ 100 MB reduction.
* `detach_per_effect_diagnostics = TRUE`: fit object size
  drops by ~L*p doubles (~5 MB at p=3500, L=10).
* `coef(fit)`, `predict(fit, X_new)`, `mf_post_smooth(fit, ...)`
  produce identical results across all combinations of
  `detach_*` flags (using only the remaining fields). Tests
  cover all 8 combinations.
* `susie_post_outcome_configuration(fit)` works when
  `detach_per_effect_diagnostics = FALSE`; errors with a clear
  message when `TRUE`.
* `summary(fit)` works in all cases (degrades gracefully when
  diagnostics are detached).

## Impact

* `mfsusie()` and `fsusie()` signatures gain 3 args (default-
  preserving — no existing user is affected unless they opt in).
* `mf_post_smooth` gets a small adapter to read the aggregate
  `mu2_agg` form when present.
* `summary.mfsusie` print method gracefully omits per-effect
  diagnostics when detached.
* `NEWS.md` documents the new flags.

## Out of scope

* Changing the defaults of any existing `attach_*` flag
  (`attach_smoothing_inputs`, `attach_lbf_variable_outcome`).
* A generic `slim()` post-hoc method (could be added separately
  but adds API surface; the `detach_*` arguments handle the
  common case at fit time).
* Compressing the retained fields (e.g. via `qs::qserialize`).
* Memory-budget-aware adaptive detach.
* `mu` itself — `mu` is needed by every post-fit consumer
  (coef, predict, all smoothers). Cannot be detached without
  losing the public API.

## Risk + mitigation

* User scripts that read `fit$residuals` directly will see
  `NULL` when `detach_residuals = TRUE`. Mitigation: default is
  `FALSE` (existing behavior preserved); document the flag.
* `detach_mu2 = TRUE` changes the fit object shape such that
  some older HMM/scalewise smoother code paths may not handle
  the aggregated form. Mitigation: small adapter in
  `mf_post_smooth` checks for `mu2_agg` first; integration tests
  exercise both paths.
* Tests in `test_post_smooth_*.R` may compare fit objects
  field-by-field; those tests should not run on slimmed
  objects. Mitigation: gate them on the `detach_*` flags and
  only run the slim variants on a focused slim-test file.
