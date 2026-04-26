# Tasks

## 1. Scalewise SD fix

- [ ] 1.1 Add `.invert_variance_curve(var_w, filter_number,
      family)` helper in `R/utils_wavelet.R` returning the
      position-space variance vector via the squared-filter-
      coefficient inverse-DWT pattern.
- [ ] 1.2 Replace the per-position SD line in
      `R/mfsusie_methods.R::.post_smooth_scalewise` with
      `sd_pos <- sqrt(.invert_variance_curve(var_w, ...))`.
- [ ] 1.3 Closed-form sanity test in
      `tests/testthat/test_post_smooth_scalewise.R`: on a
      small fixture with a known input variance vector,
      assert the returned `sd_pos` matches the exact
      `sqrt(W^2 %*% var_w)` reference.

## 2. Smash method

- [ ] 2.1 Add `R/post_smooth_smash.R` with
      `.post_smooth_smash(fit, ...)` mirroring
      `fsusieR::univariate_smash_regression`. Gates on
      `requireNamespace("smashr", quietly = TRUE)`.
- [ ] 2.2 Extend `mf_post_smooth` enum to
      `c("TI", "scalewise", "HMM", "smash")` and route
      `"smash"` to the new internal.
- [ ] 2.3 `tests/testthat/test_post_smooth_smash.R` (new):
      bit-identity vs upstream at tolerance `<= 1e-12`.

## 3. Default change + method-aware precondition

- [ ] 3.1 Change `mf_post_smooth` default from
      `"scalewise"` to `"TI"`.
- [ ] 3.2 Method-aware precondition check: `"scalewise"`
      does not require `fit$residuals` / `fit$lead_X`;
      the other three do.

## 4. Docs

- [ ] 4.1 Roxygen for `mf_post_smooth` lists all four
      methods with their dependency requirements.
- [ ] 4.2 Rewrite `vignettes/post_processing.Rmd`:
      present TI as the recommended default; remove the
      "Parseval" justification; provide the corrected
      linear-combination derivation for scalewise.
- [ ] 4.3 `NEWS.md` 0.0.1 entry mentions the default
      change.

## 5. Spec delta

- [ ] 5.1 Update or extend
      `inst/openspec/specs/mf-post-processing/spec.md`
      with the four-method enum and the corrected
      scalewise SD requirement.

## 6. Build + archive

- [ ] 6.1 `devtools::test()` passes.
- [ ] 6.2 `devtools::document()`.
- [ ] 6.3 Push, CI green.
- [ ] 6.4 `openspec archive fix-scalewise-sd-and-add-smash-method`.
