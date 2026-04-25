## ADDED Requirements

### Requirement: IBSS delegation to susieR

The package SHALL fit an mfsusie model by calling `susieR::susie_workhorse`
with an `mf_individual` data object and a `params` list, registering
`mf_individual` and `mfsusie` S3 methods at package load via
`zzz.R::.onLoad`. mfsusieR SHALL NOT reimplement the outer IBSS iteration
or the per-effect sweep.

#### Scenario: workhorse invocation

- **WHEN** `mfsusie()` is called with valid X, Y, pos
- **THEN** the call stack from `mfsusie()` to the IBSS body SHALL pass
  through `susieR::susie_workhorse` and dispatch into mfsusieR's S3
  methods (`compute_residuals.mf_individual`,
  `compute_ser_statistics.mf_individual`,
  `calculate_posterior_moments.mf_individual`, etc.)

#### Scenario: method registration on load

- **WHEN** the mfsusieR package is loaded
- **THEN** `.onLoad` SHALL register S3 methods on `mf_individual` and
  `mfsusie` into `susieR`'s namespace, mirroring the pattern at
  `mvsusieR/R/zzz.R:37-126`

### Requirement: per-(scale, modality) residual variance update by default

`update_variance_components.mf_individual` SHALL default to estimating
residual variance per (scale, modality), producing one variance per
scale per modality. The legacy mode `residual_variance_method =
"shared_per_modality"` SHALL be available and SHALL produce one variance
per modality, replicated across scales, matching `mvf.susie.alpha`.

#### Scenario: default produces per-(scale, modality) variances

- **WHEN** `mfsusie()` runs to convergence with default arguments
- **THEN** `fit$sigma2` SHALL be a list of M length-`S_m` numeric
  vectors

#### Scenario: legacy mode produces shared-per-modality variances

- **WHEN** `mfsusie()` runs with `residual_variance_method =
  "shared_per_modality"`
- **THEN** `fit$sigma2` SHALL be a list of M scalars, and the resulting
  posterior summaries (alpha, mu, mu2, lbf, KL, pip) SHALL match
  `mvf.susie.alpha::multfsusie` on the same fixed seed at tolerance
  1e-6 for posterior summaries and 1e-10 for deterministic
  intermediates, where any deviation is documented in the test with a
  cited reason

### Requirement: ELBO matches the manuscript

`get_objective.mf_individual` SHALL compute the ELBO of the mfsusie
model as defined in `methods/derivation.tex`
eq:elbo_frorm_mean_feild and eq:ERSS, including the posterior-variance
term in the expected residual sum of squares.

#### Scenario: ELBO is monotone non-decreasing on a fixed scenario

- **WHEN** the IBSS loop runs on a deterministic toy scenario
- **THEN** the ELBO trajectory `fit$elbo` SHALL be non-decreasing across
  iterations within numerical tolerance 1e-10, and SHALL be
  reproducible across runs given the same seed

### Requirement: convergence criteria

The IBSS loop SHALL converge when either the change in ELBO between
iterations is less than `tol` (default 1e-4), OR the maximum change in
alpha across iterations is less than `tol` when `convergence_method =
"pip"`. The default convergence method SHALL be `"elbo"`.

#### Scenario: ELBO-based convergence

- **WHEN** `mfsusie()` runs with `convergence_method = "elbo"` (default)
- **THEN** the loop SHALL terminate when `0 <= elbo[iter+1] -
  elbo[iter] < tol` and SHALL set `fit$converged <- TRUE`

### Requirement: greedy R-search via L_greedy passthrough

`mfsusie()` SHALL pass `L_greedy` through to `susieR::susie_workhorse`
without modification once the susieR generalization (recorded in
`design.md` "External coordination") is in place. Until that
generalization lands, mfsusieR SHALL accept the argument and SHALL
emit a one-time deprecation-style warning that the value is being
ignored.

#### Scenario: post-generalization passthrough

- **WHEN** the susieR upstream generalization is available and
  `mfsusie()` is called with `L_greedy = 3`
- **THEN** the underlying `susie_workhorse` call SHALL receive
  `L_greedy = 3` and the IBSS SHALL fit with L = 3, check pruning,
  grow L by 1, and re-fit until pruning occurs or `L` is reached

#### Scenario: pre-generalization graceful degradation

- **WHEN** the susieR upstream generalization is NOT yet available and
  `mfsusie()` is called with a non-NULL `L_greedy`
- **THEN** the function SHALL fit with L as a fixed upper bound, ignore
  `L_greedy`, and emit one warning per session referencing the
  upstream coordination plan

### Requirement: degenerate case reduces exactly to susieR::susie

`mfsusie()` SHALL reduce mathematically and numerically to
`susieR::susie()` under the degenerate input combination defined in
`design.md` D11: `M = 1`, `T_1 = 1`, `prior_variance_grid` of length 1,
`null_prior_weight = 0`, `cross_modality_prior = NULL`,
`prior_variance_scope = "per_modality"`, `residual_variance_method =
"shared_per_modality"`, `post_processing = "none"`, `L_greedy = NULL`.

#### Scenario: numerical equivalence to susieR::susie at degenerate inputs

- **WHEN** `mfsusie(X, list(matrix(y, ncol = 1)), pos = list(1), L = L,
  prior_variance_grid = sigma_0_squared, null_prior_weight = 0,
  prior_variance_scope = "per_modality", residual_variance_method =
  "shared_per_modality", post_processing = "none", L_greedy = NULL,
  ...)` and `susieR::susie(X, y, L = L, scaled_prior_variance =
  sigma_0_squared / var(y), null_weight = 0, ...)` are run with the
  same fixed seed
- **THEN** the two fits SHALL agree element-wise on `alpha`, `mu`,
  `mu2`, `lbf`, `lbf_variable`, `KL`, `sigma2`, `elbo`, `niter`, `pip`
  at tolerance 1e-10, and on credible-set membership exactly

#### Scenario: degeneracy holds across L values

- **WHEN** the degeneracy test above is repeated for `L in c(1, 5,
  10)` on the same fixture
- **THEN** the agreement SHALL hold for every L value tested

### Requirement: roxygen tags reference manuscript and original code

Every internal function ported from `mvf.susie.alpha` SHALL carry a
roxygen tag `@references_original mvf.susie.alpha/R/<file>.R#L<lo>-L<hi>`
and every function implementing a manuscript formula SHALL carry
`@manuscript_ref methods/<file>.tex eq:<label>`.

#### Scenario: tags present and accurate

- **WHEN** `roxygen2::roxygenise()` runs on the package
- **THEN** every exported function and every internal function with a
  manuscript-derived formula SHALL have at least one `@manuscript_ref`
  tag, and every ported internal function SHALL have at least one
  `@references_original` tag, both pointing at valid file paths
