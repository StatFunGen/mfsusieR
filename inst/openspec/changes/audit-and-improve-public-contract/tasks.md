# Tasks

Each section below produces its own commit. Sections are
ordered to land independent / read-only audits first, then
commit-shaped fixes.

## 1. Fit-field trim and documentation -- NO-OP

Decision: fit object stays as is. No tasks.

## 2. Feature gap with `fsusieR` / `mvf.susie.alpha`

- [ ] 2.1 List every exported function in `fsusieR`. For each,
  identify whether mfsusieR has an equivalent and where.
- [ ] 2.2 Same for `mvf.susie.alpha`.
- [ ] 2.3 List every public parameter in `fsusieR::susiF` and
  `mvf.susie.alpha::multfsusie`. Identify mfsusieR equivalent
  parameter or note as gap.
- [ ] 2.4 Write `inst/notes/sessions/<date>-feature-parity.md`
  with the table; classify each gap as `port-now` /
  `port-if-asked` / `out-of-scope`. No code change here.

## 3. Trim trivial S3 overrides

- [ ] 3.1 List every `*.mf_individual` and `*.mfsusie` S3
  method registered in `R/zzz.R`. For each, identify whether
  the body diverges from `susieR`'s default or just calls
  super-equivalent code.
- [ ] 3.2 For trivial overrides, drop the override. Re-test;
  expect 1188 PASS unchanged.
- [ ] 3.3 List anything in `susieR` that looks miscalibrated
  or redundant during the review; surface in the audit note
  for separate upstream-PR discussion.

## 4. Verbose output, susieR-arg parity, parameter renames

### 4a. Inherit susieR's check_convergence machinery

- [ ] 4a.1 Delete `check_convergence.mf_individual` from
  `R/ibss_methods.R`. The susieR `check_convergence.default`
  fires for `mf_individual` data.
- [ ] 4a.2 Run `devtools::test()` to verify nothing depends on
  the removed override. Fix any test that relied on the
  monotone-ELBO assumption (susieR also gives that, but with
  a warning rather than silent stop).
- [ ] 4a.3 Drive a small fixture with `verbose = TRUE` and
  capture the output. Confirm it shows the susieR tabular
  format including `mem`, `V`, `delta` columns.

### 4b. Expose susieR args on the public API

- [ ] 4b.1 Add `convergence_method = c("elbo", "pip")` to
  `mfsusie()`. Default `"elbo"` (matches susieR). Forward to
  `params$convergence_method`.
- [ ] 4b.2 Add `pip_stall_window = 5` to `mfsusie()`. Forward
  to `params$pip_stall_window`.
- [ ] 4b.3 Add `estimate_residual_variance = TRUE`. Forward
  to `params$estimate_residual_variance`. Verify
  `update_model_variance.mf_individual` honours the flag
  (currently always estimates).
- [ ] 4b.4 Roxygen: document each of the three new args.

### 4c. Convergence-metric agreement test

- [ ] 4c.1 Add a unit test on a standard fixture that runs
  the same fit twice (`convergence_method = "elbo"` and
  `"pip"`) and asserts they reach the same posterior within
  `tol = 1e-6` on `pip` and `tol = 1e-3` on `elbo`, and that
  the iteration counts differ by at most 2.

### 4d. Rename: mixture_weight_method -> estimate_prior_variance

- [ ] 4d.1 Add `estimate_prior_variance = TRUE` to
  `mfsusie()`. `TRUE` enables mixsqp (current default);
  `FALSE` skips the M-step and keeps pi at init values.
- [ ] 4d.2 Wire to existing internal logic: when
  `estimate_prior_variance = FALSE` set
  `params$estimate_prior_method = "none"` (skip the M-step).
- [ ] 4d.3 Deprecate `mixture_weight_method` via
  `lifecycle::deprecate_warn()`. Map old values:
  `"mixsqp"` -> `TRUE`, `"none"` -> `FALSE`. Keep working for
  one minor version.
- [ ] 4d.4 Update vignettes that reference the old name.

### 4e. Rename: lbf_min -> greed_lbf_cutoff

- [ ] 4e.1 Add `greed_lbf_cutoff = 0.1` to `mfsusie()`. Forward
  to the existing internal logic that consumed `lbf_min`.
- [ ] 4e.2 Deprecate `lbf_min` via `lifecycle::deprecate_warn()`.
- [ ] 4e.3 Update vignettes that reference `lbf_min`.

### 4f. Audit get_cs wrapper

- [ ] 4f.1 Read `get_cs.mf_individual` (R/ibss_methods.R:337).
  If the body is just `susie_get_cs(...)` with no
  customization, drop the override and let
  `get_cs.default` fire.
- [ ] 4f.2 Re-test; expect 0 failures.

## 5. HMM credible band

- [ ] 5.1 Read the manuscript section on HMM smoothing
  (fsusie methods). Derive the per-position credible band
  formula from the mixture posterior.
- [ ] 5.2 Compare derived formula to `fsusieR::fit_hmm`'s band
  output and to `mfsusieR::mf_post_smooth(method = "HMM")`'s
  current band output on a fixed seed.
- [ ] 5.3 If our formula diverges from the derivation, fix.
  If both upstream and our formula match the derivation but
  the band is genuinely wide, document the band semantics in
  the roxygen + vignette and stop suppressing display.
- [ ] 5.4 Update any visualization gating that hides the band
  unconditionally.

## 6. Performance and convergence

- [ ] 6.1 Reproduce the heavy fixture (n=84, p≈3500, M=6,
  T=128). Measure clean wall-clock with current HEAD. Skip
  if not relevant after recent perf work.
- [ ] 6.2 Profile if still slow; identify hot paths.
- [ ] 6.3 Targeted fixes only: Rfast where it cleanly wins,
  cpp11 only if profile justifies. No speculative caching.
- [ ] 6.4 Investigate the 50-iteration warning case: is it a
  pathological fixture or a real algorithmic issue?
- [ ] 6.5 Re-measure; record before/after in a session note.

## 7. Comment polish + helper promotion

- [ ] 7.1 Audit British-English spellings under `R/`,
  `tests/testthat/`, `vignettes/`, `inst/notes/`. Convert
  to American (color, behavior, optimize, …). Avoid
  algorithm names (e.g., `colour_palette` if it's a
  function name we already export).
- [ ] 7.2 Audit comments and roxygen for "port", "upstream",
  "fsusieR", "mvf.susie.alpha", "fidelity", references to
  line numbers in upstream packages. Rewrite to be
  self-contained.
- [ ] 7.3 List internal `.helper` functions in `R/`. For each,
  judge generality. Promote (drop `.` prefix) those that look
  reusable, don't export.
- [ ] 7.4 Run `devtools::test()` after the rename; expect 0
  failures.

## 8. Final integration

- [ ] 8.1 Run full `devtools::test()` on the bundle.
- [ ] 8.2 Update `inst/notes/sessions/<date>-audit-summary.md`
  with what landed and what was deferred.
- [ ] 8.3 Archive this OpenSpec change.
