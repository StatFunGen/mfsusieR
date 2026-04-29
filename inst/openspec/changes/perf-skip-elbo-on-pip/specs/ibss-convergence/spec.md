# IBSS convergence behavior contract

## ADDED Requirements

### Requirement: get_objective SHALL short-circuit on PIP convergence path

`get_objective.mfsusie(data, params, model)` SHALL return
`NA_real_` without invoking `refresh_lbf_kl.mf_individual` or
`Eloglik` when `params$convergence_method == "pip"`. The full
coherent-ELBO computation SHALL run only when
`params$convergence_method == "elbo"`.

#### Scenario: PIP convergence path skips refresh

- **GIVEN** a `mfsusie` fit invocation with
  `convergence_method = "pip"` (default)
- **WHEN** the IBSS loop calls `get_objective(data, params, model)`
  via susieR's workhorse
- **THEN** `get_objective.mfsusie` MUST return `NA_real_` and MUST
  NOT call `refresh_lbf_kl.mf_individual` or `Eloglik`
- **AND** `fit$elbo` MUST be `NA_real_` for entries written on
  the PIP path

#### Scenario: ELBO path produces coherent free energy

- **GIVEN** a `mfsusie` fit invocation with
  `convergence_method = "elbo"`
- **WHEN** `get_objective.mfsusie` is called
- **THEN** it MUST call `refresh_lbf_kl.mf_individual` to refresh
  per-effect lbf and KL against iter-final pi_V
- **AND** return `Eloglik(data, model) - sum(model$KL, na.rm = TRUE)`
- **AND** the resulting ELBO trajectory MUST be the coherent
  variational free energy (not a hybrid)
