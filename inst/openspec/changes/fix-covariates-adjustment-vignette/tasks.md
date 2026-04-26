# Tasks

## 1. Port `simu_IBSS_per_level`

- [ ] 1.1 Add `R/simulation.R` with
      `mf_simu_ibss_per_level(lev_res = 7, length_grid = 10,
      pi0, alpha = 0.8, prop_decay = 0.1)`. Body mirrors
      `fsusieR::simu_IBSS_per_level` (per-scale ash mixture
      sampler, inverse DWT to position space).
- [ ] 1.2 Argument validation: `lev_res >= 2`,
      `length_grid >= 2`, `0 <= alpha`, `0 <= prop_decay <= 1`,
      `pi0` length matches `lev_res` when supplied.
- [ ] 1.3 Roxygen with full parameter contract, return-value
      enumeration, and a runnable example. `@export`.
- [ ] 1.4 `@importFrom wavethresh wd accessD wr` and
      `@importFrom ashr normalmix` added to
      `R/mfsusieR-package.R` if not already present.

## 2. Vignette rewrite

- [ ] 2.1 Replace `vignettes/fsusie_covariates_adjustment.Rmd`
      content. Keep the existing intro paragraph verbatim.
      Remove every other section.
- [ ] 2.2 Setup chunk: load `mfsusieR`, `susieR`,
      `wavethresh`; `set.seed(2)`.
- [ ] 2.3 Generating data: port upstream simulation faithfully
      with `mf_simu_ibss_per_level`. Variables renamed to
      `genetic_signal` and `covariate_signal`. Drop unused
      `svd(...)` line. Use explicit `N3finemapping$X` access.
- [ ] 2.4 Account for covariates: call
      `mf_adjust_for_covariates(Y, Cov, method = "wavelet_eb")`,
      plot the first fitted covariate effect overlay vs
      truth, sanity-check `Y_adjusted == Y - Cov %*%
      fitted_func`.
- [ ] 2.5 Recovery scatter: two-panel
      `(Y vs genetic_signal)` and
      `(Y_corrected vs genetic_signal)`.
- [ ] 2.6 Fine-map adjusted: `fsusie(Y_corrected, X, L = 10)`
      with `mfsusie_plot(fit)`.
- [ ] 2.7 Compare with no-adjust fit. Match recovered effect
      curves to ground truth `f1`, `f2` by lead-SNP not by
      hard-coded `l` index.
- [ ] 2.8 New OLS-for-scalar section: simulate a scalar
      response, call `mf_adjust_for_covariates(Y_scalar, Cov,
      X = X, method = "ols")`, show that the function returns
      both `Y_adjusted` and `X_adjusted` (FWL).

## 3. Spec + ledger

- [ ] 3.1 `inst/openspec/changes/fix-covariates-adjustment-vignette/specs/mf-simulation-helpers/spec.md`
      with the requirement that
      `mf_simu_ibss_per_level` exists, is exported, and
      mirrors the upstream sampler.

## 4. Build + archive

- [ ] 4.1 `devtools::test()` passes.
- [ ] 4.2 `devtools::document()` regenerates man + NAMESPACE.
- [ ] 4.3 Render the rewritten vignette; PIP at the causal
      SNP improves visibly under adjustment vs no adjustment.
- [ ] 4.4 Push, CI green.
- [ ] 4.5 `openspec archive fix-covariates-adjustment-vignette`.
