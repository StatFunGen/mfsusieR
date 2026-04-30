# effective-V-consumer-coherence Specification

## Purpose
TBD - created by archiving change address-effective-V-consumer-coherence. Update Purpose after archive.
## Requirements
### Requirement: trim_null_effects coherent with V-filter

mfsusieR's `trim_null_effects.mf_individual` SHALL be coherent
with `susie_get_pip`'s `V[l] > prior_tol` filter: every effect
dropped from the PIP SHALL also have its
`alpha[l, ]`, `mu[[l]]`, `mu2[[l]]`, `lbf[l]`, `KL[l]` zeroed on
the returned fit.

#### Scenario: collapsed effect appears in both filters

- **WHEN** an IBSS fit contains an effect `l` with effective
  `V[l] < prior_tol` (the mixture pi has collapsed onto the null
  component on every (m, s) group)
- **THEN** the returned fit SHALL have `alpha[l, ] = 0`,
  `mu[[l]][[m]] = 0` for all m, `mu2[[l]][[m]] = 0` for all m,
  `lbf[l] = 0`, `KL[l] = 0`
- **AND** `susie_get_pip(model)` SHALL exclude effect `l` from
  the PIP product (no double-counting)

### Requirement: V-filter exposed in public API

The `mfsusie()` return-value documentation SHALL name the
V-filter and the `prior_tol` argument so users understand why
some IBSS effects do not contribute to the returned `pip`.

#### Scenario: user reads the pip docstring

- **WHEN** a user reads `?mfsusie` or the `pip` field's roxygen
- **THEN** the documentation SHALL describe the
  `V[l] > prior_tol` filter and reference the `prior_tol`
  argument so the user can disable the filter (set
  `prior_tol = 0`) or tighten it

