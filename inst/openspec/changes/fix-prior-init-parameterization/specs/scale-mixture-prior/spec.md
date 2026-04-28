# Scale-mixture prior init

## ADDED Requirements

### Requirement: `init_scale_mixture_prior_default()` SHALL set initial mixture weights from `null_prior_weight`

`init_scale_mixture_prior_default()` SHALL discard `ash`'s fitted
mixture weights after fitting on `(bhat, shat)`, and set the
initial pi vector via `pi_null = null_prior_weight / (K + 1)`,
`pi_k = (1 - pi_null) / (K - 1)` for `k = 2, …, K`. The same
formula is used by the user-supplied-grid path via
`distribute_mixture_weights()`.

#### Scenario: ash path with default `null_prior_weight = 2` and `K = 21` mixture components

- **WHEN** `init_scale_mixture_prior_default()` is called with
  default arguments and `ash` returns a 21-component mixture
- **THEN** the resulting `G_prior[[1]]$fitted_g$pi` SHALL satisfy
  `pi[1] = 2 / 22 ≈ 0.0909` and `pi[k] = 20 / 22 / 20 ≈ 0.0455`
  for `k = 2, …, 21`.

#### Scenario: ash path with `null_prior_weight = 16.8` matches upstream's hardcoded init

- **WHEN** `init_scale_mixture_prior_default()` is called with
  `null_prior_weight = 16.8` and `ash` returns a 21-component
  mixture (so `K = 21`)
- **THEN** the resulting `pi[1] = 16.8 / 22 ≈ 0.764` is close
  (within rounding) to upstream's `0.8`.
