# Validation and vignette delta for the clfsr / ash / TI-uniform bundle

## ADDED Requirements

### Requirement: ported smoother kernels match upstream at machine precision

Each mfsusieR smoother kernel SHALL have a reference test asserting bit-identity to its `fsusieR` upstream source at `tolerance = 1e-12` on a fixed-seed fixture.

#### Scenario: mf_smash_ash matches fsusieR::smash_lw

- **WHEN** `mfsusieR:::mf_smash_ash(noisy, noise_level)` and
  `fsusieR:::smash_lw(noisy, noise_level)` are called on the
  same input
- **THEN** `mu.est` SHALL agree at `tolerance = 1e-12`
- **AND** `mu.est.var` SHALL agree at `tolerance = 1e-12`

#### Scenario: univariate_ti_regression(uniform) matches fsusieR

- **WHEN**
  `mfsusieR:::univariate_ti_regression(Y, X, ..., scaling = "uniform")`
  and `fsusieR:::univariate_TI_regression_IS(Y, X, ...)` are
  called on the same input
- **THEN** `effect_estimate` SHALL agree at `tolerance = 1e-12`
- **AND** `cred_band` SHALL agree at `tolerance = 1e-12`

### Requirement: clfsr plot toggle exercised end-to-end

`mfsusie_plot(fit, lfsr_source = "clfsr")` SHALL run end-to-end
without error on a smoothed fit, drawing to whichever device is
active.

#### Scenario: smoke test on a smoothed multi-effect fit

- **WHEN** `mfsusie_plot(fit, lfsr_source = "clfsr")` is called
  on a fit that has been through `mf_post_smooth(method = "TI")`
- **THEN** the call SHALL return invisibly with no error and no
  warning

### Requirement: vignette surfaces the user-visible new features

`vignettes/post_processing.Rmd` SHALL contain a section that
introduces (i) `clfsr_curves` and the `lfsr_source = "clfsr"`
plot toggle, (ii) `mf_post_smooth(method = "ash")` as the
smashr-free alternative, and (iii) the `nullweight` pass-through
via `...`. Niche knobs (`scaling = "uniform"`, HMM-specific
ashr args) are NOT documented at vignette level.

#### Scenario: a new user finds the new features in the vignette

- **WHEN** a new user opens `vignettes/post_processing.Rmd`
- **THEN** the rendered vignette SHALL contain a section whose
  title names `clfsr` and `ash`
- **AND** that section SHALL include at least one runnable code
  chunk for each of the three subtopics
