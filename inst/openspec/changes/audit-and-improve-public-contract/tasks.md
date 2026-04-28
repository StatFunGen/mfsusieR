# Tasks

Each section below produces its own commit. Sections are
ordered to land independent / read-only audits first, then
commit-shaped fixes.

## 0. Workflow rule -- 3-auditor review per code-change subsection

After every code-change subsection commits (any subsection
under sections 1, 3, 4, 5, 6, 7), spawn three fresh-context
Agent reviewers in parallel:

- **Auditor 1 (numerical correctness)**: read the modified
  R code and the corresponding susieR default. Confirm
  byte-equivalence where the change is delete-and-inherit;
  for non-trivial edits, confirm the new code matches the
  manuscript or audit-derivation that justifies it. Report
  findings under 200 words.
- **Auditor 2 (API surface)**: read the public-facing
  roxygen, vignettes, and any reference test affected.
  Confirm no silent semantic shift in user-visible behavior;
  flag any rename or default change that is not documented.
  Report findings under 200 words.
- **Auditor 3 (test-suite integrity)**: read the modified
  tests. Confirm tolerances and assertions still mean what
  they claim; flag any test that relaxes a contract without
  documentation. Report findings under 200 words.

Address every finding before proceeding to the next
subsection. If a finding requires an OpenSpec scope change,
update the proposal/spec/tasks before continuing.

The audit-only sections (§2 feature parity, §3a S3-override
audit, §6.1 perf re-measurement) are read-only and skip the
3-auditor step.

## 1. Diagnosis-field cleanup

- [ ] 1.1 Retire `V` from the fit object (held at 1, not
  informative). Update `cleanup_model.mf_individual` to
  not carry it, and update any test that asserts on
  `fit$V`.
- [ ] 1.2 Add a `summary.mfsusie()` line that summarizes
  `pi_V` (e.g., per-(m, s) null mass min/median/max). One
  line of output, useful for diagnosing prior collapse.
- [ ] 1.3 Document `sigma2` shape (`list[M]` of
  length-`S_m`, or scalar when `residual_variance_scope =
  "per_outcome"`) in the `mfsusie()` `@return` block.

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

## 3. Maximize susieR backbone (delete-or-patch)

### 3a. Audit S3 overrides

- [ ] 3a.1 List every `*.mf_individual` and `*.mfsusie` S3
  method registered in `R/zzz.R`.
- [ ] 3a.2 For each, classify as `delete-and-inherit`
  (trivial wrapper / reproduces default),
  `patch-susieR-and-delete` (small divergence that susieR
  could parameterize via a hook), or `keep` (real
  divergence).
- [ ] 3a.3 Record decisions in
  `inst/notes/sessions/<date>-s3-override-audit.md` with
  one-line rationales.

### 3b. Apply delete-and-inherit

- [ ] 3b.1 For each override classified as
  `delete-and-inherit`, drop the override. Confirm
  `devtools::test()` passes.
- [ ] 3b.2 Verify on a fixture that the susieR default
  produces byte-identical output to the pre-deletion
  override.

### 3c. Plan patch-susieR-and-delete

- [ ] 3c.1 For each override classified as
  `patch-susieR-and-delete`, write a one-paragraph susieR
  patch sketch (what generic / option to add upstream).
- [ ] 3c.2 Attach to a tracking note for separate upstream
  coordination per the CLAUDE.md `feature/L_greedy`
  pattern. Do NOT edit susieR in this change without
  explicit user approval.

### 3d. Document survivors

- [ ] 3d.1 For each override classified as `keep`, add a
  one-line comment at the top of the body stating the
  divergence reason.

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

### 4g. Per-iter extra-diag columns (depends on susieR patch)

- [ ] 4g.1 Sketch the susieR patch: add a
  `format_extra_diag(model)` generic called from
  `check_convergence.default` that returns extra columns
  (default empty string).
- [ ] 4g.2 Coordinate with user on susieR PR.
- [ ] 4g.3 Once the susieR patch lands, override
  `format_extra_diag.mfsusie()` to return columns:
  `max_pi_null` (largest null-component mass across (m, s)),
  `max_KL_l`, and `n_eff` from alpha entropy. Each one
  number per iter.
- [ ] 4g.4 If the susieR patch is rejected or delayed, drop
  4g.3 and document the gap in
  `inst/notes/refactor-exceptions.md`.

## 5. HMM credible band -- drop gating, document

Per the audit, mfsusieR's HMM band formula is correct (proper
law-of-total-variance); the suppression was unjustified.

- [ ] 5.1 Identify the visualization gating that hides the
  HMM credible band display.
- [ ] 5.2 Remove the gating; populate the band normally.
- [ ] 5.3 Add a one-paragraph roxygen note in
  `mf_post_smooth()` explaining the band derivation
  (`var_w = mu2_w - mean_w^2`, inverse-DWT via `W_inv^2`,
  `qnorm((1+level)/2)` multiplier for `1-α` coverage).
- [ ] 5.4 Add a unit test that asserts the HMM band from a
  small fixture has `lower <= mean <= upper` everywhere and
  approximate `1-α` coverage on the manuscript's reference
  case.

## 6. Performance and convergence

Heavy fixture (n=84, p≈3500, M=6, T=128) **still exceeds
10 minutes** on current HEAD after cache + subsetting +
cpp11. The SER-step matrix multiplies (X^T R, ~2.3B FP ops
per IBSS iter) were untouched by prior perf work.

- [x] 6.1 Reproduce heavy fixture wall-clock (timed out at
  10 min on current HEAD; documented in
  `2026-04-28-audit-findings.md` section C).
- [ ] 6.2 Run `profvis` on a smaller-but-similar fixture
  (`n=84, p=1000, M=6, T=128`) so the profile completes;
  confirm `compute_residuals.mf_individual` /
  `compute_ser_statistics.mf_individual` matrix multiplies
  dominate. Save flamegraph under
  `bench/profiling/flamegraphs/`.
- [ ] 6.3 cpp11 port of `mf_per_outcome_bhat_shat()` core
  loop (X^T R divide-by-xtx_diag plus per-(scale, outcome)
  sigma broadcast). Pure-R reference comparison at
  `tol = 1e-12` per CLAUDE.md Phase 7 rule.
- [ ] 6.4 Re-measure heavy fixture after the cpp11 port.
  Target: under 5 minutes; ideally under 2.
- [ ] 6.5 Investigate the 50-iteration warning case: is it a
  pathological fixture or a real algorithmic issue? Address
  if real; document if pathological.
- [ ] 6.6 If still over target after (6.3), evaluate whether
  exposing `convergence_method = "pip"` from Section 4
  reduces iteration count enough; if so, document as the
  recommended setting for `p >> n` fixtures.
- [ ] 6.7 Record before/after in a session note.

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
