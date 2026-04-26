# Broaden reference tests to reach 95% / 90% coverage

## Why

The current package coverage measured by
`covr::package_coverage` is 90.83% overall. Per-file lows:

- `R/zzz.R` (0%, untestable: `.onLoad` hooks).
- `R/mfsusie_plot.R` (73.19%, 74 lines uncovered):
  multi-modality, errorbar, stack-facet, lfsr-curve overlay
  paths.
- `R/ibss_methods.R` (83.33%, 8 lines uncovered): edge
  branches in `configure_data`, `get_intercept`, `get_cs`.
- `R/mfsusie_methods.R` (92.64%, 19 lines uncovered):
  `predict.mfsusie`, `coef.mfsusie`, `mf_post_smooth`
  early-return / error-handling branches.
- `R/post_smooth_hmm.R` (94.97%, 9 lines uncovered): edge
  branches in `mf_fit_hmm` and `mf_univariate_hmm_regression`.

Target: overall `>= 95%` and per-file `>= 90%` (zzz.R
exempt). The highest-leverage additions are reference tests
against `fsusieR` and `mvf.susie.alpha` covering public-API
surfaces that lack such tests today.

## What changes

### 1. New test files

`tests/testthat/test_predict_coef_fidelity.R`. For the M=1
case, assert
`predict.mfsusie(fit, newx = X_new)` and
`coef.mfsusie(fit)` match `fsusieR::predict.susiF(fit_ref,
newx = X_new)` and `coef.susiF(fit_ref)` at tolerance
`<= 1e-12`.

`tests/testthat/test_plot_layouts.R`. Smoke tests for the
multi-modality (`M > 1`), errorbar, stack-facet, and lfsr-
curve overlay paths in `mfsusie_plot()` and
`mfsusie_plot_lfsr()`. Assertions are layout-only (panel
counts, return-value invisibility); no numerical assertions.

`tests/testthat/test_post_smooth_hmm_branches.R`. Parameter-
sweep edge cases: `halfK in {3, 5, 10, 20}`, `maxiter in
{1, 5, 50}`, `prefilter in {TRUE, FALSE}`, both with and
without the small-data prefilter triggering. Reference
comparison vs `fsusieR::fit_hmm` on each combination.

### 2. Existing test extensions

`tests/testthat/test_mfsusie_methods.R`. Add coverage for the
early-return and error-handling branches of `predict.mfsusie`,
`coef.mfsusie`, and `mf_post_smooth` flagged at lines 39, 275,
343, 346, 350, 419-420, 464, 480, 595-598, 646, 648.

`tests/testthat/test_ibss_methods.R` (new or extend
existing). Cover the edge branches in `configure_data`,
`get_intercept`, `get_cs`, and `get_variable_names`.

### 3. Coverage CI

Confirm the existing coverage CI step in
`.github/workflows/dispatch_pkgdown_build.yml` (the
`Compute coverage` step) reports the overall percentage to
the shields.io endpoint after this change. No workflow
changes required.

## Impact

- New: three test files plus extensions to existing.
- No source-code changes.
- Coverage badge reflects the new percentage automatically.
- Specs: no new spec capability; this change is a pure test
  expansion.

## Out of scope

- Plot snapshot testing via `vdiffr`. Layout assertions are
  shape-only.
- Coverage of `R/zzz.R` (`.onLoad` is untestable through covr).
