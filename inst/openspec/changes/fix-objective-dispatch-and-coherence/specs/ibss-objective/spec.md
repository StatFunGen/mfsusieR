# IBSS objective contract

## ADDED Requirements

### Requirement: The IBSS objective SHALL use susieR's default path

mfsusieR SHALL NOT register a `get_objective` method for the
`"mfsusie"` or `"mf_individual"` classes. The per-iteration objective
SHALL be susieR's `get_objective.default`, which evaluates
`Eloglik.mf_individual(data, model) - sum(model$KL, na.rm = TRUE)` under
both `convergence_method` values.

#### Scenario: No objective override is registered

- **GIVEN** the loaded mfsusieR namespace
- **WHEN** `getS3method("get_objective", cls, optional = TRUE,
  envir = asNamespace("susieR"))` is queried for `cls` in
  {`"mfsusie"`, `"mf_individual"`}
- **THEN** both lookups MUST return `NULL`

#### Scenario: ELBO is reported under PIP convergence

- **GIVEN** `mfsusie(..., convergence_method = "pip")`
- **WHEN** the fit returns
- **THEN** `fit$elbo[-1]` MUST be finite
- **AND** convergence MUST still be decided by the PIP/alpha rule

### Requirement: Per-effect state SHALL be consistent with the returned prior

After each effect's update in an IBSS iteration, `alpha[l, ]`,
`mu[[l]]`, `mu2[[l]]`, `lbf[l]` and `KL[l]` SHALL be computed against
the same per-effect prior that is stored in `fitted_g_per_effect[[l]]`
and `pi_V[[l]]`, and `fitted` SHALL equal
`X %*% sum_l alpha_l * mu_l` per outcome.

#### Scenario: Inner EM hits its cap

- **GIVEN** an effect whose inner EM runs `max_inner_em_steps + 1`
  M-steps without meeting the lbf tolerance
- **WHEN** `post_loglik_prior_hook.mf_individual` returns
- **THEN** a final E-step MUST have run against the last M-step's prior
- **AND** `KL[l]` MUST equal the SER identity recomputed from the stored
  state at `<= 1e-8`

#### Scenario: Trimmed effects leave fitted values consistent

- **GIVEN** an effect with `V[l] < prior_tol` at the end of the fit
- **WHEN** `trim_null_effects.mf_individual` zeroes its moments
- **THEN** `fitted[[m]]` MUST equal `X %*% sum_l alpha_l * mu_l[[m]]`
  at `<= 1e-8` for every outcome
