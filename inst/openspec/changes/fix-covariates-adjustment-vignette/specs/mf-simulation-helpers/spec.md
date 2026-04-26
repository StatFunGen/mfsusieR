# mf-simulation-helpers capability

## ADDED Requirements

### Requirement: `mf_simu_ibss_per_level` is exported and matches the upstream IBSS-prior sampler

The package SHALL export a function `mf_simu_ibss_per_level(lev_res, length_grid, pi0, alpha, prop_decay)` that draws a function on a length-`2^lev_res` grid by sampling each wavelet scale's coefficients from a per-scale ash normal mixture and applying the inverse DWT. The implementation MUST mirror `fsusieR::simu_IBSS_per_level` mathematically.

#### Scenario: signature

```
mf_simu_ibss_per_level(
  lev_res     = 7,
  length_grid = 10,
  pi0,
  alpha       = 0.8,
  prop_decay  = 0.1
)
```

`pi0` is optional. When omitted it defaults to `1 - exp(-prop_decay * 1:lev_res)`.

#### Scenario: argument validation

The function SHALL error when `lev_res < 2`, `length_grid < 2`, `alpha < 0`, `prop_decay < 0`, `prop_decay > 1`, or when a supplied `pi0` does not have length `lev_res`.

### Requirement: returned object exposes the sampled function and metadata

The returned list SHALL contain `sim_func` (numeric vector of length `2^lev_res`), `true_coef` (the wavelet-domain detail coefficients used to construct the sample), `mix_per_scale` (length-`lev_res` list of `ashr::normalmix` objects), and `emp_pi0` (length-`lev_res` numeric vector of empirical zero-fractions per scale).

#### Scenario: shape contract

`length(out$sim_func) == 2^lev_res`. `length(out$mix_per_scale) == lev_res`. `length(out$emp_pi0) == lev_res`.
