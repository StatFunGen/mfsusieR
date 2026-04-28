# Tasks

Each section below produces its own commit. Sections are
ordered to land independent / read-only audits first, then
commit-shaped fixes.

## 1. Fit-field trim and documentation

- [ ] 1.1 Enumerate every field on a fresh `mfsusie()` fit on
  the `practical_dataset` fixture; record names, classes, and
  size.
- [ ] 1.2 Cross-reference each field with `susieR::susie()`
  output and `mvsusieR::mvsusie()` output.
- [ ] 1.3 Classify each field: keep / drop / move-to-summary /
  rename. Record decisions in a table inside
  `inst/notes/sessions/<date>-fit-field-audit.md`.
- [ ] 1.4 Apply the decisions: edit `R/mfsusie.R` (return-value
  block), `R/ibss_methods.R` (`cleanup_model.mf_individual`),
  any S3 method that constructs the fit. Update tests that
  read removed fields.
- [ ] 1.5 Roxygen: every kept field gets a one-line description
  in the `@return` block of `mfsusie()`.

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

## 4. Useful verbose output

- [ ] 4.1 Read `susieR`'s per-iter verbose output. Read
  `fsusieR::susiF` and `mvf.susie.alpha::multfsusie`'s
  per-iter prints.
- [ ] 4.2 Distill into a per-iter line for mfsusieR. Suggested
  format: `iter K | elbo=… | max_alpha_diff=… | sigma2_range=… [| MB=…]`.
- [ ] 4.3 Add `convergence_metric = c("pip_diff", "elbo")`
  argument to `mfsusie()`. Default `"pip_diff"` (current).
  Plumb through to `check_convergence`.
- [ ] 4.4 Add a unit test on a standard fixture: PIP-based
  and ELBO-based convergence agree within 2 iterations.

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
