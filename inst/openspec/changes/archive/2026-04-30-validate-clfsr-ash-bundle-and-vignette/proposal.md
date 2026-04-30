# Validate clfsr / ash / TI-uniform bundle and surface in the vignette

## Why

The post-smooth augmentation bundle (`add-clfsr-smash-lw-ti-uniform`)
introduced three new code paths that ship without
upstream-fidelity tests at the C2 standard the rest of the package
holds:

- `mf_smash_ash` (cycle-spinning + per-coefficient `ashr::ash` on
  the wavelet decomposition). Direct port of `fsusieR:::smash_lw`.
- `univariate_ti_regression(scaling = "uniform")`. Direct port of
  `fsusieR:::univariate_TI_regression_IS`.
- `clfsr_curves` per smoother payload. Per-variant lfsr from the
  SuSiE posterior moments; no upstream equivalent (mfsusieR-only).

Three corresponding gaps:

1. The two ported kernels have no bit-identity test against their
   fsusieR sources. Port-source-bug-fix discipline (per
   `inst/notes/refactor-discipline.md`) requires apple-to-apple
   `tolerance <= 1e-8` (or 1e-12 when explicit) tests on every
   port.
2. The `clfsr_curves` slot is defensible by formula but its
   plot-overlay path (`mfsusie_plot(lfsr_source = "clfsr")`)
   has no smoke test.
3. The vignette `post_processing.Rmd` does not expose any of the
   new features. New users see only the original four methods +
   lfsr.

## What changes

### 1. Reference fidelity tests (machine precision)

Two new test files, each a bit-identity assertion against the
fsusieR port source on a fixed-seed fixture:

- `tests/testthat/test_smash_ash_fidelity.R`. Asserts
  `mf_smash_ash(noisy, noise_level)` equals
  `fsusieR:::smash_lw(noisy, noise_level)` at
  `tolerance = 1e-12` for both `mu.est` and `mu.est.var`.
  Multiple fixtures: power-of-2 length, non-power-of-2 length
  (exercises the reflective padding), constant noise level,
  per-position noise level.
- `tests/testthat/test_ti_uniform_fidelity.R`. Asserts the
  position-space output of
  `univariate_ti_regression(Y, X, ..., scaling = "uniform")`
  equals `fsusieR:::univariate_TI_regression_IS(Y, X, ...)` at
  `tolerance = 1e-12` for both `effect_estimate` and
  `cred_band`. Fixed-seed fixtures with M = 1, T_m in {64, 128}.

Both files `skip_if_not_installed("fsusieR")` (port source is a
read-only sibling repo, not in CRAN). Both reach into fsusieR
internals via `:::` because neither helper is exported.

### 2. clfsr plot-overlay smoke test

A new test in `tests/testthat/test_clfsr.R` (extends the
existing file): assert that
`mfsusie_plot(lfsr_source = "clfsr")` produces no error and
draws to a closed PDF device (visual content not inspected).

### 3. Vignette additions

`vignettes/post_processing.Rmd` gains one new short section:

- Title: `Per-variant lfsr (clfsr) and the smashr-free smoother`.
- Three subsections, each ~5-10 lines:
  - `clfsr_curves`: per-variant view of where each effect is
    credibly nonzero, derived from the SuSiE posterior moments.
    Show the plot toggle `mfsusie_plot(lfsr_source = "clfsr")`.
  - `method = "ash"`: smashr-free alternative to
    `method = "smash"`. Same shape of payload; runs without
    the `smashr` Suggests dependency.
  - Forwarding `nullweight` (and other ashr knobs) via `...`
    on `mf_post_smooth(method = "TI" | "HMM" | "ash", ...)`.
    One-line example.

Vignette stays terse. The TI `scaling = "uniform"` mode and
the HMM `nullweight` knob are NOT documented at vignette level
(too niche for new users). Reference-docs roxygen already
covers them.

## Impact

- **Severity**: HIGH for the fidelity tests (these are the C2
  contract for ported routines and have been slipping). LOW for
  the vignette additions.
- **Source-code changes**: zero (all bundled features already
  shipped in `add-clfsr-smash-lw-ti-uniform`).
- **Tests**: 2 new fidelity test files + 1 plot smoke test
  appended to `test_clfsr.R`.
- **Documentation**: ~20-30 lines added to
  `vignettes/post_processing.Rmd`.
- **Dependencies**: no DESCRIPTION change. fsusieR remains a
  read-only sibling for fidelity comparisons (test-only).

## Out of scope

- Upstream-fidelity tests for `clfsr_curves` (no upstream
  equivalent — fsusieR's `cal_clfsr` is wavelet-domain on
  Bhat/Shat, semantically different).
- Performance / Rcpp hot-path optimization.
- Reorganizing the vignette structure beyond appending the new
  section.
- Reference tests for the `...` plumbing through to ash
  internals (already covered by `test_smash_lw.R` and
  `test_ti_uniform.R` at the behavioral level).

## References

- Source bundle:
  `inst/openspec/changes/add-clfsr-smash-lw-ti-uniform/`.
- Refactor discipline (binary tolerance philosophy + Pattern A/B
  port-source-bug fix exception):
  `inst/notes/refactor-discipline.md`.
- Existing fidelity-test exemplar:
  `tests/testthat/test_post_smooth_smash.R`
  (test 1: mfsusieR vs fsusieR at `tolerance = 1e-12`).
