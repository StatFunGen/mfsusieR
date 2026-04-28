# M-step performance contract

## ADDED Requirements

### Requirement: per-effect M-step SHALL reuse cached per-iter quantities

`optimize_prior_variance.mf_individual()` SHALL read
`sigma2_per_pos`, `shat2_m`, and `sdmat[m, s]` from
`model$em_cache` instead of recomputing them per effect.
`update_variance_components.mf_individual()` SHALL repopulate the
cache after each `sigma2` update. Numerical results MUST be
identical to the pre-cache implementation at `tol = 1e-12` on
the standard test fixtures.

#### Scenario: cached invariants are reused across effects

- **WHEN** `optimize_prior_variance.mf_individual()` is called
  for `l = 1, …, L` in one IBSS iteration with the
  cache populated
- **THEN** the function MUST NOT call `outer(svec^2, sd_grid^2, +)`
  or `outer(1 / xtx_diag, sigma2_per_pos)` more than once per
  `(m, s)` pair within that iteration.

### Requirement: adaptive variant subsetting SHALL exclude near-zero alpha SNPs

`optimize_prior_variance.mf_individual()` SHALL build the L-mat
input only from SNPs `j` with `model$alpha[l, j] > mixsqp_alpha_eps`
when `mixsqp_alpha_eps > 0` (default `1e-6`). The full-set
behaviour is recoverable by setting `mixsqp_alpha_eps = 0`.

#### Scenario: concentrated alpha drops most SNPs

- **WHEN** `model$alpha[l, ]` has 990 entries below `1e-6` and
  10 entries above (sparse posterior)
- **THEN** the L-mat passed to mixsqp MUST have at most
  10 × `length(idx)` rows.

#### Scenario: `mixsqp_alpha_eps = 0` reproduces the unsubsetted M-step

- **WHEN** `optimize_prior_variance.mf_individual()` is called
  with `mixsqp_alpha_eps = 0`
- **THEN** the resulting `model$pi_V` MUST match the
  pre-subsetting implementation at `tol = 1e-12`.
