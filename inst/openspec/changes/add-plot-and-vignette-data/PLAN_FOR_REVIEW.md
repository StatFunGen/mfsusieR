# PLAN — please approve before I implement anything

You said:
- vignette must show all three methods,
- cpp11 must accelerate the slow real-data loops,
- unit tests must verify numerical identity at machine precision
  for every ported reference implementation,
- both `fsusie()` and `mfsusie()` fits must work,
- no approximations.

This plan honors all five.

## What ships in this change

### A. Faithful port of the wavelet-variance helpers

Three helpers move from upstream wavelet-utils into
`R/utils_wavelet.R`:

- `wd.var(data, filter.number, family, type = "station", ...)`
  — like `wavethresh::wd` but with squared filter coefficients,
  yielding per-coefficient variances.
- `convert.var(wd, ...)` — like `wavethresh::convert` but on a
  variance object.
- `AvBasis.var(wst, ...)` — like `wavethresh::AvBasis` but on a
  variance object; reaches into `.C("av_basisWRAP", ...,
  PACKAGE = "wavethresh")` via `R_GetCCallable` (or direct .C
  with `PACKAGE`).

These three give the exact pointwise credible band for TI. **No
approximation anywhere.** Ports are line-for-line; no algorithmic
edits.

### B. `mf_post_smooth(fit, method = c("scalewise", "TI", "HMM"))`

Single function, fit-only interface (no X / Y). Already in tree;
**TI's band swaps to the exact `wd.var`/`convert.var`/`AvBasis.var`
path** (no Parseval approximation).

Both `fsusie()` (M = 1) and `mfsusie()` (M >= 1) fits supported by
the same implementation: the dispatch loop iterates over `m in
1:M` and `l in 1:L` regardless of M.

### C. cpp11 acceleration of the hot inner loops

I will profile the three smoothers on a realistic-sized fit
(`n = 500`, `T_basis = 256`, `L = 5`, `M = 5`) **before** writing
any cpp11. Target: identify the loops dominating wall time and
port only those.

Likely candidates (to confirm by profiling):
- TI: per-row stationary-wavelet transform — already C-backed in
  wavethresh; no cpp11 win.
- TI: regression of wavelet coefficients on the lead variable
  (`crossprod` + `colMeans(resid^2)`) — already vectorised.
- TI: scalewise `ashr::ash` — already C-backed.
- **HMM forward-backward**: pure-R loop over `T_basis`
  positions × `2*halfK+1` grid. `O(T * K)` per effect per
  outcome. **This is the most likely cpp11 win**; expect 5–20x on
  large `T`.
- **HMM emission matrix**: pure-R `dnorm` over `T x K` cells. cpp11
  port halves memory + runtime.

Ports go in `src/post_smooth.cpp` and are exposed via R wrappers
in `R/mfsusie_methods.R`. Pure-R reference oracles stay in
`R/reference_implementations.R` (per the established cpp11/oracle
pattern in mfsusieR; see `mixture_log_bf_per_scale_R`).

If profiling shows TI is also a bottleneck, the inner `wd` call
(stationary transform) is wavethresh's C internal, so the cpp11
window is limited; we'd port the surrounding R loop instead.

### D. Reference-test fidelity at machine precision

Three new test files:

- `test_post_smooth_wavelet_var.R` — bit-identity at `<= 1e-14`
  for `wd.var`/`convert.var`/`AvBasis.var` vs the upstream
  reference (skip-if-not-installed).
- `test_post_smooth_TI.R` — bit-identity at `<= 1e-14` for
  `mf_post_smooth(fit, method = "TI")$effect_curves[[1]][[l]]`
  and `$credible_bands[[1]][[l]]` vs the upstream
  `univariate_TI_regression` output, on a single-effect M=1 fit.
  Skip-if-not-installed.
- `test_post_smooth_HMM.R` — bit-identity at `<= 1e-14` for
  `mf_post_smooth(fit, method = "HMM")$effect_curves[[1]][[l]]`
  and `$lfsr_curves[[1]][[l]]` vs upstream
  `univariate_HMM_regression`. Skip-if-not-installed.

Plus `test_post_smooth_cpp_oracle.R` — bit-identity (`<= 1e-14`)
between the cpp11 kernel(s) and the pure-R oracle for the HMM
forward-backward and any other cpp11 port. Always-on test
(no skip).

Plus `test_post_smooth_smoke.R` — runs `mf_post_smooth(...)` on
both an `fsusie()` (M=1) and an `mfsusie()` (M=2) fit; checks
shapes, no errors, both supported.

### E. Vignette `post_processing.Rmd` shows all three

Side-by-side `mfsusie_plot()` panels:
- raw (`coef()` only, no smoothing) — wiggly, no bands
- `mf_post_smooth(method = "scalewise")` — fast smoothing, bands
- `mf_post_smooth(method = "TI")` — exact TI bands
- `mf_post_smooth(method = "HMM")` — bands + lfsr overlay

Demonstrates both single-outcome (`fsusie()`) and multi-outcome
(`mfsusie()`) fits. Uses packaged data
(`data(fsusie_methyl_example)`, `data(mfsusie_joint_example)`) so
the figures are stable across knit runs.

### F. Always-save residuals + lead_X

- `fit$residuals[[m]]` and `fit$lead_X[[l]]` populated on every
  fit. `save_residuals` argument removed (already in tree, tests
  pass).

### G. Audit agents

After every implementation step, I run an audit agent:
- "Verify that `mf_post_smooth(method = 'TI')` gives results
  bit-identical to the upstream reference at `<= 1e-14` on the
  test fixtures"
- "Verify that the cpp11 kernel and the pure-R oracle agree at
  `<= 1e-14` on randomized inputs (1000 trials)"

Findings go in this PLAN as commit notes.

## Order of execution

1. Port `wd.var` / `convert.var` / `AvBasis.var` into
   `R/utils_wavelet.R`. Test at `<= 1e-14` vs upstream.
2. Replace TI's approximate band with the exact path. Test TI's
   point estimate AND band at `<= 1e-14`.
3. Implement HMM forward-backward + emission matrix in cpp11.
   Pure-R oracle stays in `R/reference_implementations.R`.
4. Profile TI and HMM on the realistic fixture; port any
   additional hot loop only if it dominates wall time.
5. Bit-identity tests for all three smoothers vs upstream at
   machine precision.
6. Smoke tests for both `fsusie()` and `mfsusie()` paths.
7. Build packaged data fixtures (`data/fsusie_methyl_example.rda`
   + `data/mfsusie_joint_example.rda`) — trimmed,
   de-identified.
8. Refresh `post_processing.Rmd` to show all 3 methods on both
   single and multi outcome fixtures, with before/after panels.
9. Run audit agents.

## What is NOT in this change

Nothing. Everything in the list above ships in this change.
No "future PR" anywhere in this plan.

## I will start when you say "approved"

Reply "approved" or list edits. I will not touch code until you
approve.
