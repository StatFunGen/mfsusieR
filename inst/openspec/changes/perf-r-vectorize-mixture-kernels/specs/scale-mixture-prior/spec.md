# Mixture-of-normals kernels implementation contract

## ADDED Requirements

### Requirement: per-(j, t, k) log density SHALL be computed once per (l, m, scope group)

`mf_em_log_dens_per_scale` SHALL compute the per-(j, t, k)
mixture-component log density
`log dnorm(bhat[j, t]; 0, sqrt(shat[j, t]^2 + sd_k^2 * V))`
exactly once per (effect l, outcome m, scope group s) within the
SER step. The cache MUST be reused by the
downstream consumers (mixsqp likelihood matrix, mixture log Bayes
factor aggregation, posterior moment closed form) without
recomputation. The cache SHALL NOT persist across (l, m)
iterations -- it is ephemeral to the current SER step. Concrete
consumers that MUST read the cache:

1. The mixsqp likelihood matrix construction
2. The per-(j, t) mixture log Bayes factor aggregation
3. The per-(j, t) posterior mean / second-moment closed form

The cache SHALL NOT persist across (l, m) iterations; it is
ephemeral to the current SER step.

#### Scenario: log_dens computed once per SER step

- **GIVEN** an IBSS iteration entering
  `optimize_prior_variance.mf_individual` for effect l
- **WHEN** `compute_ser_statistics.mf_individual` runs
- **THEN** `mf_em_log_dens_per_scale` MUST be called exactly
  `M` times (once per outcome) per scope group s
- **AND** the resulting `log_dens` matrix MUST be reused without
  recomputation by `mf_em_likelihood_per_scale`,
  `loglik.mf_individual`, and `calculate_posterior_moments.mf_individual`
  within the same SER step

### Requirement: cpp kernels for mixture density SHALL be removed

The repository SHALL NOT contain the cpp implementations
`mixture_log_bf_per_scale_cpp`,
`mixture_posterior_per_scale_cpp`, and
`mf_em_log_likelihood_per_scale_cpp`. The vectorized R
implementations in `R/posterior_mixture.R` and `R/em_helpers.R`
are the sole production path.

#### Scenario: src/ contains no per-cell mixture kernels

- **GIVEN** a clean checkout
- **WHEN** the source tree is inspected
- **THEN** `src/posterior_mixture.cpp` MUST NOT contain the
  three named kernels
- **AND** `R/cpp11.R` MUST NOT export R wrappers for these
  kernels

## MODIFIED Requirements

### Requirement: numerical parity with reference R implementation

The production R implementation MUST match the reference R
oracles (`mixture_log_bf_per_scale_R`,
`mixture_posterior_per_scale_R`) at `tol = 1e-12`. The reference
oracles SHALL remain in the repo. The `tol = 1e-12` is elevated
from the prior cpp-vs-R parity threshold because the production
path is now itself the (vectorized) R reference and benefits from
R's IEEE-754 deterministic ordering.

#### Scenario: production R matches reference R at machine precision

- **WHEN** running the reference parity tests on a fixed-seed
  fixture
- **THEN** PIPs, alpha, mu, mu2, lbf MUST agree at
  `tol = 1e-12`
