# Tasks

Severity: HIGH for Section 1 (port-source fidelity, currently
missing). LOW for Section 2 (vignette polish).

## 1. Reference fidelity tests

- [ ] 1.1 `tests/testthat/test_smash_ash_fidelity.R`. Three
      `test_that` blocks:
      - `mf_smash_ash` matches `fsusieR:::smash_lw` on a
        power-of-2 fixture, scalar `noise_level`, fixed seed.
        `tolerance = 1e-12` on both `mu.est` and `mu.est.var`.
      - Same on a NON-power-of-2 fixture (exercises the
        reflective padding branch).
      - Same with a per-position `noise_level` vector.
- [ ] 1.2 `tests/testthat/test_ti_uniform_fidelity.R`. Two
      `test_that` blocks:
      - `univariate_ti_regression(scaling = "uniform")` matches
        `fsusieR:::univariate_TI_regression_IS` on T_m = 64
        fixed-seed fixture. `tolerance = 1e-12` on
        `effect_estimate` and `cred_band`.
      - Same at T_m = 128.
- [ ] 1.3 Both test files start with
      `testthat::skip_if_not_installed("fsusieR")`. Both reach
      into fsusieR internals via `:::` (helpers are not
      exported).

## 2. clfsr plot-overlay smoke test

- [ ] 2.1 Append to `tests/testthat/test_clfsr.R`: a test that
      runs `mfsusie_plot(fit, lfsr_source = "clfsr")` to a PDF
      device and confirms it returns invisibly without error.

## 3. Vignette additions

- [ ] 3.1 Append a section
      `## Per-variant lfsr (clfsr) and the smashr-free smoother`
      to `vignettes/post_processing.Rmd`. Three short
      subsections:
      - clfsr (~10 lines + 1 plot).
      - `method = "ash"` (~5 lines + 1 code chunk).
      - Forwarding `nullweight` via `...` (1 chunk, 3 lines of
        prose).
- [ ] 3.2 Build the vignette locally and confirm it renders
      without errors / warnings.

## 4. Validate

- [ ] 4.1 `openspec validate validate-clfsr-ash-bundle-and-vignette`
      from `inst/`.
- [ ] 4.2 Full `devtools::test()` clean.
- [ ] 4.3 `R CMD check --as-cran` returns 0 errors, 0 warnings on
      the new files.
