# Fit-object retention contract

## ADDED Requirements

### Requirement: detach_residuals SHALL drop IBSS-internal buffers

`mfsusie()` and `fsusie()` SHALL accept a `detach_residuals`
argument (default `FALSE`). When `detach_residuals = TRUE`, the
returned fit object MUST NOT contain `$residuals`,
`$raw_residuals`, or `$fitted_without_l`. `$fitted` MUST be
retained (used by `predict.mfsusie`).

#### Scenario: residuals dropped on opt-in

- **GIVEN** `mfsusie(..., detach_residuals = TRUE)` returns
  `fit`
- **THEN** `is.null(fit$residuals)` AND
  `is.null(fit$raw_residuals)` AND
  `is.null(fit$fitted_without_l)` MUST all be TRUE
- **AND** `is.null(fit$fitted)` MUST be FALSE
- **AND** `predict(fit, X_new)` MUST produce identical output
  to a fit with `detach_residuals = FALSE`

### Requirement: detach_mu2 SHALL replace mu2 with alpha-weighted aggregate

When `detach_mu2 = TRUE`, the returned fit MUST replace each
`fit$mu2[[l]][[m]]` (p × T_basis[m] matrix) with the
alpha-weighted length-T_basis[m] aggregate
`colSums(fit$alpha[l, ] * fit$mu2[[l]][[m]])`, stored at
`fit$mu2_agg[[l]][[m]]`. `fit$mu2` MUST be set to NULL.

`mf_post_smooth` paths that consume `mu2` MUST check for
`mu2_agg` first and use the aggregate form when present.

#### Scenario: mu2 replaced with aggregate

- **GIVEN** `mfsusie(..., detach_mu2 = TRUE)` returns `fit`
- **THEN** `is.null(fit$mu2)` MUST be TRUE
- **AND** `fit$mu2_agg[[l]][[m]]` MUST be a numeric vector of
  length `data$T_basis[m]` for every (l, m)
- **AND** `mf_post_smooth(fit, method = "HMM")` MUST produce
  identical output to a fit with `detach_mu2 = FALSE`
- **AND** `mf_post_smooth(fit, method = "scalewise")` MUST
  produce identical output to a fit with `detach_mu2 = FALSE`

### Requirement: detach_per_effect_diagnostics SHALL drop diagnostic fields

`mfsusie()` and `fsusie()` SHALL accept a
`detach_per_effect_diagnostics` argument (default `FALSE`).
When `detach_per_effect_diagnostics = TRUE`, the returned fit
MUST drop `$lbf`, `$lbf_variable`, and `$lbf_variable_outcome`.

#### Scenario: diagnostics dropped, downstream errors clearly

- **GIVEN** `mfsusie(..., detach_per_effect_diagnostics = TRUE)`
  returns `fit`
- **THEN** `is.null(fit$lbf)` AND `is.null(fit$lbf_variable)`
  AND `is.null(fit$lbf_variable_outcome)` MUST all be TRUE
- **AND** `summary(fit)` MUST run without error and omit
  per-effect lbf rows
- **AND** `susie_post_outcome_configuration(fit)` MUST raise an
  error whose message names the `detach_per_effect_diagnostics`
  flag (so the user knows what to flip)

### Requirement: backward compatibility for default arguments

All three `detach_*` flags SHALL default to `FALSE`. The fit
object returned by default-argument `mfsusie()` / `fsusie()`
calls MUST be bit-identical to the pre-change implementation
(modulo unrelated perf changes). Existing user scripts MUST
continue to work without modification.

#### Scenario: default-arg fits unchanged

- **GIVEN** a script that calls `mfsusie(X, Y_list, L = 5)`
  before this change
- **WHEN** the same call is made after this change with no
  new arguments
- **THEN** the returned fit object MUST contain the same
  fields with the same shapes as before (`$mu`, `$mu2`,
  `$residuals`, `$raw_residuals`, `$fitted_without_l`,
  `$fitted`, `$alpha`, `$lbf`, `$lbf_variable`, etc.)
- **AND** PIPs / alpha / mu MUST agree at `tol = 1e-12`
  with the pre-change fit on the same fixture / seed
