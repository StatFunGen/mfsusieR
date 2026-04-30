# Add per-variant clfsr, smash flavor switch, TI uniform mode

## Why

Three small surface additions surfaced by the 2026-04-30 cross-package
coverage analysis as upstream features mfsusieR did not expose.
Bundled because they are independent and each is small.

## What changes

### 1. clfsr (per-variant credible-level lfsr)

A new top-level slot `fit$clfsr[[l]][[m]]` of shape `p x T_basis[m]`
holding the per-variant lfsr derived from the per-variant SuSiE
posterior moments `mu` and `mu2`. Computed once at the end of
`mfsusie()` from
`pnorm(-|mu[[l]][[m]]| / sqrt(pmax(mu2[[l]][[m]] - mu[[l]][[m]]^2, 0)))`.

A new internal helper `clfsr_from_gaussian(mu, sd)` mirrors
`lfsr_from_gaussian` and is matrix-aware.

No change to the existing `lfsr_curves` slot or to the plot defaults.
Documentation is one sentence: clfsr is the per-variant version of
lfsr.

### 2. smash flavor switch + ash kernel

A new arg `flavor = c("ash", "smashr")` on `mf_post_smooth(method = "smash")`
selects the smoother kernel. Default `"ash"`.

- `flavor = "ash"`: cycle-spinning + per-coefficient `ashr::ash` on
  the wavelet coefficients of the per-position OLS estimate, with
  the per-position OLS Shat used as the noise level (no
  user-provided value). Pure dependencies: `wavethresh` + `ashr`,
  both already in `Imports`. Does not touch `smashr`.
- `flavor = "smashr"`: existing `smashr::smash.gaus` path. Gated on
  `requireNamespace("smashr")` as today.

Internal helper: a new `mf_smash_lw(noisy_signal, noise_level, ...)`
function ports the cycle-spinning + ash kernel.
`univariate_smash_regression` gains the same `flavor` argument and
forwards to one of the two kernels.

`smashr` stays in Suggests (current state). The `lw` default avoids
the smashr dependency on installs that do not have it; the
`requireNamespace` gate moves to the `flavor = "smashr"` path only.

### 3. TI uniform mode

`univariate_ti_regression` gains one arg
`scaling = c("per_scale", "uniform")` (default `"per_scale"`).

- `"per_scale"` (default, current behavior): scale Y via
  `col_scale`, per-scale `ashr::ash` loop with `nullweight = 30`
  per scale, final rescale by `csd_Y / csd_X`.
- `"uniform"`: do not scale Y, single `ashr::ash` across all
  wavelet coefficients with `nullweight = 3`, final rescale by
  `1 / csd_X`.

Selected from `mf_post_smooth(method = "TI", scaling = "uniform")`.
No new method enum.

## Impact

- **Severity**: low. Three independent additions; each preserves the
  existing default behavior except where the user explicitly opts in.
- **Source-code changes**:
  - `R/mfsusie_methods.R` (clfsr_from_gaussian helper, post-fit
    clfsr loop, univariate_ti_regression `scaling` arg, mf_post_smooth
    `flavor` and `scaling` args).
  - `R/post_smooth_smash.R` (mf_smash_lw helper,
    univariate_smash_regression `flavor` arg).
- **Tests**:
  - `tests/testthat/test_clfsr.R`: shape + lfsr-vs-clfsr aggregation
    identity at machine precision.
  - `tests/testthat/test_smash_lw.R`: shape + parity against the
    `gaus` flavor on a controlled fixture (loose tolerance, the two
    smoothers differ by design).
  - `tests/testthat/test_ti_uniform.R`: shape + the three structural
    differences vs `scaling = "per_scale"` (raw Y not scaled, single
    ash call, output differs).
- **Dependencies**: `smashr` stays Suggests (no DESCRIPTION change).
  CI / pixi.toml unchanged.
- **Documentation**: roxygen tweaks per added arg; no NEWS entry
  pre-alpha.

## Out of scope

- Per-variant position-space smoothing (would be p x cost; not
  needed for clfsr because the per-variant version is wavelet-domain).
- Lead-variant-only post-processing (the existing alpha-aggregated
  smoothers already collapse to lead-variant when alpha concentrates).
- Removing `smashr` Suggests entirely (kept for `flavor = "smashr"`).
- Generalizing other smoothers to the same uniform / per-scale
  knob.

## References

- 2026-04-30 cross-package coverage analysis (see chat transcript).
- Existing audit OpenSpec change:
  `inst/openspec/changes/audit-cross-package-post-hooks/`.
- Existing post-smooth code:
  `R/mfsusie_methods.R::mf_post_smooth`,
  `R/post_smooth_smash.R::univariate_smash_regression`.
