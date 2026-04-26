# Tasks

## 1. Predict / coef fidelity

- [ ] 1.1 `tests/testthat/test_predict_coef_fidelity.R` (new).
      M=1 case bit-identity vs
      `fsusieR::predict.susiF` / `coef.susiF` across
      `n in {100, 200}`, `T in {64, 128}` at tolerance
      `<= 1e-12`. Gated `skip_if_not_installed("fsusieR")`.

## 2. Plot layouts

- [ ] 2.1 `tests/testthat/test_plot_layouts.R` (new).
      Smoke + layout assertions on the multi-modality,
      errorbar, stack-facet, and lfsr-curve overlay paths.

## 3. HMM branches

- [ ] 3.1 `tests/testthat/test_post_smooth_hmm_branches.R`
      (new). Parameter-sweep reference tests vs
      `fsusieR::fit_hmm` exercising the prefilter and
      maxiter edge branches.

## 4. mfsusie_methods branches

- [ ] 4.1 Extend `tests/testthat/test_mfsusie_methods.R` (or
      create if absent) covering the early-return / error-
      handling paths in `predict.mfsusie`, `coef.mfsusie`,
      `mf_post_smooth`.

## 5. ibss_methods branches

- [ ] 5.1 Add or extend tests for `configure_data`,
      `get_intercept`, `get_cs`, `get_variable_names` edge
      branches.

## 6. Verification

- [ ] 6.1 Run `covr::package_coverage()` locally and
      confirm overall `>= 95%`, per-file `>= 90%`
      (zzz.R exempt).
- [ ] 6.2 Push, confirm coverage badge updates after the
      next pkgdown CI run.

## 7. Archive

- [ ] 7.1 `openspec archive broaden-reference-tests-coverage`.
