# `mixsqp_alpha_eps` default change (delta)

## MODIFIED Requirements

### Requirement: `mixsqp_alpha_eps` default SHALL be `5e-5`

The default value of `mixsqp()`'s alpha-thinning threshold,
`mixsqp_alpha_eps`, is changed from `1e-6` to `5e-5`. The
escape hatch `mixsqp_alpha_eps = 0` (use every variant
regardless of alpha) MUST continue to work.

#### Scenario: default value is `5e-5`

- **WHEN** `mfsusie()` is called without supplying
  `mixsqp_alpha_eps`
- **THEN** `formals(mfsusie)$mixsqp_alpha_eps` MUST evaluate
  to `5e-5`.

#### Scenario: bumped default does not change fits at practical tolerance

- **WHEN** the canonical fixtures (`why_functional`,
  `susie_post_outcome_configuration` end-to-end) are fit at
  the new default
- **THEN** `max|alpha_new - alpha_old|` MUST be ≤ `1e-3`
  and `max|pip_new - pip_old|` MUST be ≤ `1e-3` against the
  same fits at `mixsqp_alpha_eps = 1e-6`. Credible-set
  membership and `niter` MUST be identical.

#### Scenario: `mixsqp_alpha_eps = 0` still uses every variant

- **WHEN** `mfsusie(mixsqp_alpha_eps = 0)` is called on a
  fixture
- **THEN** the M-step input MUST include every variant
  (no thinning), reproducing the un-thinned behavior at
  `tol = 1e-12`.
