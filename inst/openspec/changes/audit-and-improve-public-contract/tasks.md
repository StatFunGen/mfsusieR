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

## 2. Feature gap port-now items

The audit (2026-04-28-audit-findings.md section B)
identified the gaps. This section ports them.

### 2a. Post-processing summary helper

- [ ] 2a.1 New function `mf_summarize_effects(fit, ...)` (or
  similar name) that takes a fit (with credible bands
  populated by `mf_post_smooth()`) and returns a list of
  per-(CS, outcome) position ranges where the band excludes
  zero. Reads existing `fit$smoothed[[method]]$credible_bands`
  and `fit$sets$cs`.
- [ ] 2a.2 Roxygen + a unit test on a fixture with known
  causal SNPs at known positions.

### 2b. Yuan 2024 post-hoc causal configurations

- [ ] 2b.1 New function `mf_posthoc_configurations(fit, ...)`
  that takes a fit and returns posterior probability over
  the 2^L causal configurations.
- [ ] 2b.2 Implement in pure R first.
- [ ] 2b.3 If a real fixture measures slow at typical L
  (say > 100 ms), evaluate cpp11 port. For L <= 15 expect
  pure R sufficient.

### 2c. `simu_IBSS_ash_vanilla` simulation helper

- [ ] 2c.1 Port the body from
  `~/GIT/fsusieR/R/Simulations_functions.R` into a new
  `data-raw/` script following the existing simulation
  helper conventions.
- [ ] 2c.2 If a fixture is shipped, add a `data/` entry
  with documentation.

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

### 3b. Apply delete-and-inherit -- mfsusieR (commit)

- [ ] 3b.1 Drop the 9 overrides identified in the audit
  (`validate_prior`, `check_convergence`, `configure_data`,
  `get_cs`, `neg_loglik`, `get_alpha_l.mfsusie`,
  `set_prior_variance_l.mfsusie`,
  `get_prior_variance_l.mfsusie`, `get_variable_names`).
  Verify `get_variable_names` falls through to `.individual`
  byte-equivalently before deletion.
- [ ] 3b.2 Run `devtools::test()`; expect 1188 PASS
  unchanged (or +/- a few that read deleted fields).
- [ ] 3b.3 Run the 3-auditor humanize step.
- [ ] 3b.4 Commit.

### 3c. Apply delete-and-inherit -- mvsusieR (stage, do not commit)

For ~5 mvsusieR overrides identified by the sidebar audit
(`neg_loglik.mv_individual`, `get_cs.mv_individual`,
`get_cs.mv_ss` -- a safety upgrade since `.ss` uses
`safe_cov2cor` instead of `cov2cor`, `get_alpha_l.mvsusie`,
`get_prior_variance_l.mvsusie`, `set_prior_variance_l.mvsusie`):

- [ ] 3c.1 Apply the deletions in `~/GIT/mvsusieR/R/`
  (working tree changes only).
- [ ] 3c.2 Run mvsusieR's test suite; verify no regression.
- [ ] 3c.3 Stop. Do not commit. User reviews and commits
  manually (matching the user's stated workflow for
  cross-package changes).

### 3d. susieR patches -- stage, do not commit

For each `patch-susieR-and-delete` candidate (P1, P2, P3,
optionally P4-P6), apply the patch to `~/GIT/susieR/R/` on
the master branch:

- [ ] 3d.1 P1: `format_sigma2_summary(model)` and
  `format_extra_diag(model)` generics in
  `check_convergence.default`. Default
  `format_sigma2_summary` returns `sprintf("%.4f", model$sigma2)`;
  default `format_extra_diag` returns `""`.
- [ ] 3d.2 P2: change `sum(model$KL)` to
  `sum(model$KL, na.rm = TRUE)` in `get_objective.default`.
- [ ] 3d.3 P3: add `cleanup_extra_fields(data)` generic
  called from `cleanup_model.default`; default returns
  `character(0)`.
- [ ] 3d.4 Run susieR's test suite; verify no regression.
- [ ] 3d.5 Stop. Do not commit. User reviews and commits
  manually.

### 3e. After P1 lands -- delete check_convergence override

- [ ] 3e.1 Once user commits the susieR P1 patch, delete
  `check_convergence.mf_individual`.
- [ ] 3e.2 Override `format_sigma2_summary.mf_individual`
  to return the compact list-of-vectors summary.
- [ ] 3e.3 Override `format_extra_diag.mf_individual` to
  return `max_pi_null` / `max_KL` / `n_eff`.
- [ ] 3e.4 Run mfsusieR test suite + 3-auditor humanize
  step.
- [ ] 3e.5 Commit.

### 3f. After P2/P3 land -- delete more overrides

- [ ] 3f.1 Once user commits the susieR P2/P3 patches,
  delete `get_objective.mf_individual` and
  `cleanup_model.mf_individual`.
- [ ] 3f.2 Run tests + 3-auditor humanize step.
- [ ] 3f.3 Commit.

### 3g. Document survivors

- [ ] 3g.1 For each override classified as `keep`, add a
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

## 5. HMM credible band -- remove suppression

- [ ] 5.1 Identify the gating that hides the HMM band.
- [ ] 5.2 Remove the gating; populate `credible_bands` the
  same way TI / smash / scalewise smoothers do.
- [ ] 5.3 Run tests + 3-auditor humanize step.
- [ ] 5.4 Commit. No roxygen note (pre-alpha; bug fix).

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
