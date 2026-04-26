# Tasks

## 1. Helper additions

- [ ] 1.1 Add `mf_low_count_indices(Y_wd, threshold = 0)` to
      `R/utils_wavelet.R`. Body: `which(apply(abs(Y_wd), 2,
      median) <= threshold)`. Roxygen with manuscript
      reference.
- [ ] 1.2 Add `mf_quantile_normalize(Y_wd)` to
      `R/utils_wavelet.R`. Body: column-wise
      `qqnorm(rank(.x, ties.method = "random"))$x` with
      `set.seed(1)` for tie reproducibility. Roxygen with
      manuscript reference.
- [ ] 1.3 Add internal
      `mf_compute_bhat_shat_with_low_count_mask(X, Y_wd,
      lowc_idx)` to `R/individual_data_methods.R`. Body:
      single call to `susieR::compute_marginal_bhat_shat`
      followed by the post-compute mask
      `Bhat[, lowc_idx] <- 0; Shat[, lowc_idx] <- 1`.
- [ ] 1.4 Unit tests in
      `tests/testthat/test_preprocessing_helpers.R`:
      bit-identity of each helper against the upstream
      equivalent on randomized inputs.

## 2. Pipeline integration

- [ ] 2.1 Extend `create_mf_individual` in
      `R/individual_data_class.R` to accept
      `low_count_filter` and `quantile_norm` arguments.
      After the per-outcome DWT, compute
      `lowc_idx[[m]]`; if `quantile_norm` is `TRUE`, apply
      `mf_quantile_normalize` to `Y_wd[[m]]`. Store
      `lowc_idx[[m]]` on the data object.
- [ ] 2.2 Update `R/individual_data_methods.R` SER step to
      use the masked Bhat/Shat helper.
- [ ] 2.3 Verify `init_scale_mixture_prior_default` accepts
      `lowc_idx`; if not, extend its signature to honor the
      mask in the prior init.
- [ ] 2.4 Verify mixture-weight EM
      (`R/em_helpers.R::mf_em_likelihood_per_scale` and
      `mf_em_m_step_per_scale`) zeros out the contribution
      of masked columns. If the existing helpers already
      ignore zero Bhat / one Shat columns, no change; if
      not, add a `lowc_idx` argument.

## 3. Public API

- [ ] 3.1 `R/mfsusie.R`: add `low_count_filter = 0` and
      `quantile_norm = FALSE` to the public signature with
      roxygen.
- [ ] 3.2 `R/fsusie.R`: forward both arguments to
      `mfsusie()`.
- [ ] 3.3 `R/adjust_covariates.R`: lift the v1 rejects on
      `thresh_lowcount > 0` and `quantile_trans = TRUE`;
      replace with full support, plumbed through
      `mf_residualize_wavelet_eb`. When `quantile_norm =
      TRUE`, invert the transform on the returned
      `Y_adjusted` so the downstream `fsusie(Y_adjusted,
      X)` call works on the original units.
- [ ] 3.4 Update `man/*.Rd` via `devtools::document()`.

## 4. Tests

- [ ] 4.1 `tests/testthat/test_preprocessing.R` (new): bit-
      identity sweep across `n in {84, 200}`,
      `T in {64, 128, 256}`,
      `low_count_filter in {0, 0.1, 0.5}`,
      `quantile_norm in {FALSE, TRUE}` vs the upstream
      functional fine-mapping result. Tolerance
      `<= 1e-12`. Skipped when fsusieR is not installed.
- [ ] 4.2 Targeted test for the zero-median-column case at
      `low_count_filter = 0` that the pre-port mfsusieR
      misses.
- [ ] 4.3 Targeted test for `mf_adjust_for_covariates(
      quantile_norm = TRUE)`: returned `Y_adjusted` on the
      original scale; downstream `fsusie(Y_adjusted, X)`
      converges and recovers a known causal SNP.

## 5. Vignettes

- [ ] 5.1 `vignettes/mfsusie_intro.Rmd`: one-paragraph
      mention of `low_count_filter` and `quantile_norm`
      with a pointer to
      `vignette("practical_data_applications")`.
- [ ] 5.2 `vignettes/fsusie_intro.Rmd`: same.

## 6. Spec delta

- [ ] 6.1 Create
      `inst/openspec/changes/add-preprocessing-low-count-and-quantile-norm/specs/mf-preprocessing/spec.md`
      with requirements + scenarios for the two arguments
      and the masked-Bhat helper.

## 7. Build + archive

- [ ] 7.1 `devtools::test()` passes.
- [ ] 7.2 `devtools::document()` regenerates man + NAMESPACE.
- [ ] 7.3 Push, CI green.
- [ ] 7.4 `openspec archive add-preprocessing-low-count-and-quantile-norm`.
