# mfsusieR public contract

## ADDED Requirements

### Requirement: every S3 override SHALL fall into one of three classes

Every `*.mf_individual` or `*.mfsusie` S3 override SHALL be
classified as `delete-and-inherit` (drop the override; susieR's
default fires), `patch-susieR-and-delete` (open an upstream
hook in susieR, then delete the override once the upstream
patch lands), or `keep` (real, irreducible divergence;
documented with a one-line rationale comment in the body).
Overrides whose body trivially wraps the default SHALL NOT
remain in the codebase.

#### Scenario: a removed override defers to susieR's default

- **WHEN** mfsusieR removes its trivial override of a generic
  `g`
- **THEN** `g(data, params, model, ...)` for an
  `mf_individual` data object SHALL dispatch to
  `g.default()` from susieR
- **AND** the resulting `mfsusie()` fit on the standard test
  fixtures SHALL be byte-identical to the pre-removal fit.

### Requirement: the `V` field SHALL be retired from the fit object

`mfsusie()` SHALL NOT carry a `V` field on the returned object.
Per-effect prior adaptation in mfsusieR is in the mixture
weights `pi_V`; the scalar `V` was held at 1 for all effects
and provided no diagnostic value. `summary.mfsusie()` SHALL
provide a one-line summary of `pi_V` (e.g., null-component
mass quantiles across (m, s) groups) instead.

#### Scenario: V field absent

- **WHEN** the user calls `names(mfsusie(X, Y, ...))`
- **THEN** `"V"` SHALL NOT be in the returned names.

### Requirement: `verbose = TRUE` SHALL emit susieR's per-iteration tabular output

`mfsusie()` SHALL inherit susieR's `check_convergence.default`
verbose output. Per-iteration lines SHALL include at minimum
the iteration number, ELBO, ELBO delta, sigma2 summary, peak
memory, and prior-variance summary, in tabular form. The
existing `check_convergence.mf_individual` override SHALL be
removed.

#### Scenario: tabular per-iter output

- **WHEN** `mfsusie(X, Y, verbose = TRUE)` is called
- **THEN** for each IBSS iteration the printed output SHALL
  include the columns `iter`, `ELBO`, `delta`, `sigma2`,
  `mem`, `V` (in that order)
- **AND** the header row SHALL be printed once before the
  iteration body executes.

### Requirement: convergence behaviour SHALL be selectable

`mfsusie()` SHALL accept `convergence_method = c("elbo", "pip")`
(default `"elbo"`) and `pip_stall_window` (default `5`),
forwarding both to `params`. The susieR scaffold's PIP-based
fallback path with stall detection SHALL be reachable.

#### Scenario: PIP-based convergence reaches same posterior as ELBO

- **WHEN** the user runs `mfsusie(X, Y, convergence_method =
  "pip")` and again with `convergence_method = "elbo"` on a
  standard fixture
- **THEN** the two fits SHALL agree on `pip` within
  `1e-6` and on the final `elbo` within `1e-3`
- **AND** the iteration counts SHALL differ by at most 2.

### Requirement: `estimate_residual_variance` SHALL be user-controllable

`mfsusie()` SHALL accept `estimate_residual_variance` (default
`TRUE`). When `FALSE`, `update_model_variance` SHALL leave
`model$sigma2` at its initial value across all IBSS iterations.

#### Scenario: fixed sigma2 with FALSE

- **WHEN** the user runs `mfsusie(X, Y,
  estimate_residual_variance = FALSE, residual_variance =
  s2_init)`
- **THEN** the final fit's `sigma2` SHALL be byte-identical to
  `s2_init` for every outcome.

### Requirement: `mixture_weight_method` SHALL be replaced by `estimate_prior_variance`

`mfsusie()` SHALL accept `estimate_prior_variance` (default
`TRUE`) replacing the previous `mixture_weight_method`
argument. `TRUE` SHALL enable the mixsqp M-step (former
`"mixsqp"`); `FALSE` SHALL skip the M-step and hold mixture
weights at their init values (former `"none"`). The old
argument SHALL emit a `lifecycle::deprecate_warn()` and remain
functional for one minor version with the obvious value
mapping.

#### Scenario: deprecation warning on old argument

- **WHEN** the user calls `mfsusie(X, Y, mixture_weight_method
  = "mixsqp")`
- **THEN** a `lifecycle::deprecate_warn()` SHALL fire pointing
  at `estimate_prior_variance`
- **AND** the fit SHALL run as if `estimate_prior_variance =
  TRUE` had been passed.

### Requirement: `lbf_min` SHALL be renamed to `greed_lbf_cutoff`

`mfsusie()` SHALL accept `greed_lbf_cutoff` (default `0.1`)
replacing `lbf_min`. The renamed argument SHALL gate the
L-greedy outer loop the same way `lbf_min` did. `lbf_min`
SHALL emit a `lifecycle::deprecate_warn()` and remain
functional for one minor version.

#### Scenario: deprecation warning on old argument

- **WHEN** the user calls `mfsusie(X, Y, lbf_min = 0.3)`
- **THEN** a `lifecycle::deprecate_warn()` SHALL fire pointing
  at `greed_lbf_cutoff`
- **AND** the fit SHALL run as if `greed_lbf_cutoff = 0.3` had
  been passed.

### Requirement: HMM smoother SHALL return a credible band

`mf_post_smooth(method = "HMM")` SHALL populate
`fit$smoothed$HMM$credible_bands[[m]][[l]]` with a per-
position lower/upper bound at the requested coverage level.
The band formula SHALL be derived from the HMM posterior
mixture and SHALL be documented in the function roxygen.

#### Scenario: HMM band populated

- **WHEN** the user runs `mf_post_smooth(fit, method = "HMM")`
  on a fit with at least one credible set
- **THEN** for each `(m, l)` with a non-null effect curve,
  `fit$smoothed$HMM$credible_bands[[m]][[l]]` SHALL be a
  `T_basis[m] x 2` matrix
- **AND** the lower-bound column SHALL not exceed the
  upper-bound column at any position.
