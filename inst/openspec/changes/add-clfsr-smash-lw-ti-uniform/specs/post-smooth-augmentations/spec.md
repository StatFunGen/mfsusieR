# Post-smooth augmentations spec delta

## ADDED Requirements

### Requirement: per-variant clfsr in the smoother payload

Every `mf_post_smooth(method = ...)` call SHALL populate a
`clfsr_curves[[m]][[l]]` slot of shape `p x T_basis[m]` on the
returned smoother payload alongside `lfsr_curves`. Each entry
holds the per-variant credible-level lfsr derived from the
per-variant SuSiE posterior moments `mu` and `mu2`.

#### Scenario: clfsr matches the per-variant Gaussian formula

- **WHEN** `mf_post_smooth(fit, method = ...)` returns a fit with
  populated `clfsr_curves`
- **THEN** `fit$smoothed[[method]]$clfsr_curves[[m]][[l]][j, t]`
  SHALL equal
  `pnorm(-|fit$mu[[l]][[m]][j, t]| / sqrt(pmax(fit$mu2[[l]][[m]][j, t] - fit$mu[[l]][[m]][j, t]^2, 0)))`
  at machine precision for every `(method, l, m, j, t)`
- **AND** every entry of
  `fit$smoothed[[method]]$clfsr_curves[[m]][[l]]` SHALL lie in
  `[0, 0.5]`

#### Scenario: plot consumes either lfsr or clfsr

- **WHEN** `mfsusie_plot(fit, lfsr_source = "clfsr")` is called on
  a fit that has been smoothed
- **THEN** the plot's lfsr overlay SHALL use the alpha-aggregated
  per-variant clfsr (`alpha[l, ] %*% clfsr_curves[[m]][[l]]`)
  truncated to the position-grid length, instead of the
  smoother's position-space lfsr

### Requirement: smash and ash as separate post-smooth methods

`mf_post_smooth(fit, method = ...)` SHALL accept `"smash"` and
`"ash"` as distinct method names instead of a `flavor` argument.
Method-specific options SHALL flow through `...` to the
underlying shrinkage tool: `ashr::ash` for `"ash"`,
`smashr::smash.gaus` for `"smash"`.

- `"smash"` SHALL invoke `smashr::smash.gaus` and gate on
  `requireNamespace("smashr")`.
- `"ash"` SHALL invoke a cycle-spinning + per-coefficient
  `ashr::ash` kernel that uses the per-position OLS Shat as the
  noise level. This path SHALL NOT require the `smashr` package.

The `mf_post_smooth(method = "ash", nullweight = K, ...)`
spelling subsumes fsusieR's orphan `smash_lw` (default
`nullweight = 300`) and `smash_2lw` (mild `nullweight`) helpers
without exposing them as separate functions.

#### Scenario: ash method on a smashr-less install

- **WHEN** `mf_post_smooth(method = "ash")` is called on a system
  where `smashr` is not installed
- **THEN** the call SHALL succeed and populate
  `fit$smoothed[["ash"]]`

#### Scenario: smash method on a smashr-less install

- **WHEN** `mf_post_smooth(method = "smash")` is called on a
  system where `smashr` is not installed
- **THEN** the call SHALL error with a message naming `smashr` and
  pointing at `method = "ash"` as the alternative

### Requirement: TI uniform scaling mode

The TI smoother SHALL expose a `scaling = c("per_scale", "uniform")`
argument on both `univariate_ti_regression` and
`mf_post_smooth(method = "TI", ...)`, defaulting to `"per_scale"`.
Each value MUST select a distinct internal branch:
`"per_scale"` SHALL column-scale Y, run a per-scale
`ashr::ash` loop with `nullweight = 30`, and rescale the output
by `csd_Y / csd_X`. `"uniform"` SHALL skip Y scaling, run a
single `ashr::ash` across all wavelet coefficients with
`nullweight = 3`, and rescale the output by `1 / csd_X`.

#### Scenario: uniform mode preserves high-frequency signal

- **WHEN** `mf_post_smooth(method = "TI", scaling = "uniform")` and
  `mf_post_smooth(method = "TI", scaling = "per_scale")` are run on
  the same fit
- **THEN** the two outputs SHALL differ
- **AND** the `"uniform"` output SHALL apply less aggressive
  shrinkage on detail coefficients than the `"per_scale"` output
  (verifiable by L2 norm of the wavelet-detail residual)
