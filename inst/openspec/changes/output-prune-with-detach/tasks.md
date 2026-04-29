# Tasks

## 1. Audit consumers of fit-object fields

- [ ] 1.1 In `inst/notes/sessions/<date>-fit-fields-audit.md`,
  document for every field on the returned `mfsusie` /
  `fsusie` fit:
  - which public-facing functions read it (`coef`, `predict`,
    `mf_post_smooth`, `summary`, `susie_post_outcome_configuration`,
    `mfsusie_plot`, `mfsusie_plot_lfsr`)
  - approximate size at the reference fixture
  - whether it can be safely dropped after the IBSS loop
- [ ] 1.2 Identify any dependencies in vignettes / tests that
  read fields beyond the documented "public" set; flag any
  hidden consumers.

## 2. Add `detach_residuals` flag

- [ ] 2.1 Add `detach_residuals = FALSE` to `mfsusie()` and
  `fsusie()` signatures.
- [ ] 2.2 Plumb to `params$detach_residuals` in the
  `susie_workhorse` call.
- [ ] 2.3 In `mfsusie()` post-fit (after `susie_workhorse`
  returns), if `detach_residuals = TRUE`, set
  `fit$residuals <- NULL`, `fit$raw_residuals <- NULL`,
  `fit$fitted_without_l <- NULL`.
- [ ] 2.4 Test: `coef(fit)`, `predict(fit, X_new)`, and
  `mf_post_smooth(fit, method = "TI")` work identically with
  and without the detach.
- [ ] 2.5 Test: object size reduction matches expectation on
  the reference fixture.

## 3. Add `detach_mu2` flag with mu2_agg fallback

- [ ] 3.1 Add `detach_mu2 = FALSE` arg to `mfsusie()` and
  `fsusie()`.
- [ ] 3.2 In `mfsusie()` post-fit, if `detach_mu2 = TRUE`,
  iterate over (l, m) and replace
  `fit$mu2[[l]][[m]]` with the alpha-weighted aggregate
  `colSums(fit$alpha[l, ] * fit$mu2[[l]][[m]])` (length T_basis).
  Store at `fit$mu2_agg[[l]][[m]]`. Set `fit$mu2 <- NULL`.
- [ ] 3.3 Modify `mf_post_smooth` paths that consume `mu2`
  (HMM and scalewise smoothers; trace via grep) to check for
  `fit$mu2_agg` and fall back to the aggregate form.
- [ ] 3.4 Test: each smoother method produces identical
  results with `detach_mu2 = TRUE` vs `FALSE` (on a
  fixed-seed fixture).

## 4. Add `detach_per_effect_diagnostics` flag

- [ ] 4.1 Add `detach_per_effect_diagnostics = FALSE` arg.
- [ ] 4.2 If TRUE, set `fit$lbf <- NULL`,
  `fit$lbf_variable <- NULL`, and
  `fit$lbf_variable_outcome <- NULL`.
- [ ] 4.3 Modify `summary.mfsusie` print method to omit
  the lbf summary when these are NULL.
- [ ] 4.4 Modify `susie_post_outcome_configuration` to
  raise a clear error when `lbf_variable_outcome` is NULL
  (already does; just verify the message names the
  detach flag).

## 5. Documentation

- [ ] 5.1 Roxygen for each new arg, including a recommendation
  on when to set it `TRUE` (production fits / large datasets /
  saveRDS persistence).
- [ ] 5.2 `NEWS.md` entry.
- [ ] 5.3 Vignette note in `mfsusie_intro.Rmd` about object
  size and the `detach_*` flags.

## 6. Tests

- [ ] 6.1 New file
  `tests/testthat/test_detach_options.R` covering all 8
  combinations of the three flags. Each test:
  - fits the same fixture twice (with and without the flag)
  - asserts identical PIPs / alpha / fitted
  - asserts which fields are NULL
  - asserts post-fit operations still work
- [ ] 6.2 Object-size assertion: with all three flags TRUE,
  fit size on the reference fixture is < 50% of with all
  three flags FALSE.

## 7. Commit

- [ ] 7.1 Single focused commit
  `feat(output): detach_* arguments to slim fit objects`.
