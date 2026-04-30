# Tasks

Severity: LOW. Three independent additions. Each task lands as a
standalone PR; no ordering requirement.

## 1. clfsr slot

- [ ] 1.1 Add `clfsr_from_gaussian(mu, sd)` internal helper next to
      `lfsr_from_gaussian` in `R/mfsusie_methods.R`. Same formula,
      matrix-friendly.
- [ ] 1.2 In `R/mfsusie.R` post-fit block (right before the
      `attach_smoothing_inputs` step), populate
      `fit$clfsr[[l]][[m]]` for every (l, m) from `fit$mu` /
      `fit$mu2` via `clfsr_from_gaussian`.
- [ ] 1.3 New test `tests/testthat/test_clfsr.R`: shape
      `p x T_basis[m]` per (l, m); aggregation identity
      `sum_j alpha[l, j] * clfsr[j, t]` reproduces a single-pass
      lfsr at machine precision.
- [ ] 1.4 One-sentence roxygen on `clfsr_from_gaussian`.

## 2. smash flavor switch + ash kernel

- [ ] 2.1 Add `mf_smash_lw(noisy_signal, noise_level, n.shifts = 50,
      filter.number = 1L, family = "DaubExPhase")` in
      `R/post_smooth_smash.R`. Cycle-spinning + per-coefficient
      `ashr::ash` with reflective padding for non-power-of-2
      lengths.
- [ ] 2.2 Add `flavor = c("ash", "smashr")` arg to
      `univariate_smash_regression`. Default `"ash"`. Forward to the
      ported kernel (`mf_smash_lw`) using the OLS Shat as the noise
      level; `"smashr"` keeps the existing `smashr::smash.gaus` call.
- [ ] 2.3 Add `flavor = c("ash", "smashr")` arg to `mf_post_smooth`.
      Forward to `.post_smooth_smash` and from there to
      `univariate_smash_regression`.
- [ ] 2.4 Move the `requireNamespace("smashr")` gate in
      `mf_post_smooth` so it fires only when `flavor = "smashr"`.
- [ ] 2.5 New test `tests/testthat/test_smash_lw.R`: shape;
      structural parity with `flavor = "smashr"` on a small fixture
      (loose tolerance — the two smoothers are different
      algorithms).

## 3. TI uniform mode

- [ ] 3.1 Add `scaling = c("per_scale", "uniform")` arg to
      `univariate_ti_regression` in `R/mfsusie_methods.R`. Default
      `"per_scale"`. Branch internally:
      - `"per_scale"`: current path (scale Y, per-scale ash with
        nullweight = 30, rescale by `csd_Y / csd_X`).
      - `"uniform"`: skip Y scaling, single ash with nullweight = 3,
        rescale by `1 / csd_X`.
- [ ] 3.2 Forward `scaling` from `.post_smooth_ti` and
      `mf_post_smooth(method = "TI", scaling = ...)`.
- [ ] 3.3 New test `tests/testthat/test_ti_uniform.R`: shape;
      `scaling = "uniform"` differs from `scaling = "per_scale"` on
      a controlled fixture; the differences match the three
      documented branch points (Y not scaled, single ash, output
      differs).

## 4. Validate

- [ ] 4.1 `openspec validate add-clfsr-smash-lw-ti-uniform` from
      `inst/`.
- [ ] 4.2 Full `devtools::test()` clean.
