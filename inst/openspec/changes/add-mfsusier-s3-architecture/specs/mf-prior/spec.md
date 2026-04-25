## ADDED Requirements

### Requirement: scale-mixture-of-normals prior with per-(scale, modality) parameters

The package SHALL provide a default prior class `mf_prior_scale_mixture`
that stores per-(scale, modality) mixture weights `pi_{k,s,m}` and
per-(scale, modality) prior variances `sigma_{k,s,m}^2`, matching the
manuscript model defined in `methods/derivation.tex` eq:additive_model.
The default null-component weight SHALL be `null_prior_weight = 2`,
chosen as the post-scaling-equivalent of `mvf.susie.alpha`'s
`nullweight = 0.7` * `max(K) = 3` and exposed directly so the
multiplication-by-K step in the original code is removed.

#### Scenario: default null_prior_weight is 2

- **WHEN** `mfsusie()` is called without specifying `null_prior_weight`
- **THEN** the prior object's `null_prior_weight` field SHALL equal 2
  and the EM/mixsqp updates SHALL use this value directly without any
  additional scaling by K

#### Scenario: per-(scale, modality) shapes

- **WHEN** an `mf_prior_scale_mixture` object is constructed with M = 3
  modalities and `T_padded = c(64, 128, 1024)`
- **THEN** the object SHALL store `pi` as a length-3 list of matrices
  with shapes `(log2(64) x K)`, `(log2(128) x K)`, `(log2(1024) x K)`,
  and `V_grid` SHALL be a length-3 list of length-K vectors (or a
  length-3 list of `S_m x K` matrices when `prior_variance_scope =
  "per_scale_modality"` is requested)

#### Scenario: per-modality scope as legacy mode

- **WHEN** the user constructs the prior with `prior_variance_scope =
  "per_modality"`
- **THEN** the prior SHALL collapse the scale dimension and store
  `V_grid` as a length-M list of length-K vectors, replicated per
  scale during SER computation, matching `mvf.susie.alpha`'s
  per-trait treatment

### Requirement: per-modality init via susieR helper + ash fit

The per-modality scale-mixture-of-normals init logic SHALL:

1. Honor `prior_variance_grid` when supplied: use the user-given
   K-vector grid directly, distribute mixture weights with
   `null_prior_weight` taking weight `null_prior_weight / (K + 1)`
   on the null component (V = 0) and the remaining `K` non-null
   components splitting `1 - null_prior_weight / (K + 1)`
   uniformly. (For the susieR-degeneracy contract C1, the user
   supplies `prior_variance_grid` of length 1 and
   `null_prior_weight = 0`, producing a single-component Gaussian
   prior identical to `susieR::susie`'s `scaled_prior_variance`
   form.)
2. Fall back to data-driven ash fit when `prior_variance_grid` is
   `NULL`: per modality, compute `Bhat`, `Shat` via
   `susieR::compute_marginal_bhat_shat(X, Y_m)` (the new shared
   helper landing on `~/GIT/susieR` master per design.md
   "Migration Plan"), draw a sample, fit `ashr::ash` to obtain
   the K-vector grid, then assemble `V_grid` and `pi` with the
   `null_prior_weight` distribution above.

There SHALL be NO `T_m == 1` branch in the init logic. The same
code path produces a valid prior for any `T_m`, including the
univariate degenerate case. The public `mfsusie()` SHALL emit
`warning_message` when any `T_m <= 3` advising the user to
consider `susieR::susie` (for `T_m == 1`) or `fsusieR::susiF`
(for short functional traits) before proceeding.

The routine is internal to mfsusieR. mfsusieR SHALL NOT have
`fsusieR` in `DESCRIPTION` Imports and SHALL NOT call
`fsusieR::*` at runtime.

#### Scenario: user-supplied prior_variance_grid is honored

- **WHEN** `mf_prior_scale_mixture(data, prior_variance_grid =
  v_grid, null_prior_weight = w_null)` is called for any
  `T_m >= 1`
- **THEN** the returned prior's `V_grid` SHALL equal `v_grid`
  per modality (replicated where appropriate by
  `prior_variance_scope`), and `pi` SHALL be the
  `null_prior_weight` distribution; no ash fit SHALL run

#### Scenario: data-driven ash fit when prior_variance_grid is NULL

- **WHEN** `mf_prior_scale_mixture(data, prior_variance_grid =
  NULL)` is called for `M = 1, T_1 > 1`
- **THEN** `compute_marginal_bhat_shat(X, data$D[[1]])` SHALL
  feed `ashr::ash` with the same arguments and seed sequence
  that `fsusieR::init_prior.default` uses, producing
  element-wise identical `pi`, `sd`, `mean` fields at tolerance
  `<= 1e-12` (validated by the C2 contract in `mf-ibss/spec.md`)

#### Scenario: warning at very-short T_m

- **WHEN** `mfsusie()` (the public API) is called with any
  `T_m <= 3` (including `T_m == 1`)
- **THEN** `mfsusie()` SHALL emit `susieR::warning_message`
  advising consideration of `susieR::susie` / `fsusieR::susiF`,
  but SHALL still complete the fit through the standard data-
  driven init path

#### Scenario: port-quality audit logged

- **WHEN** the susieR helper integration and the prior init land
  in a Phase 3 PR
- **THEN** the PR SHALL include a port-quality audit per
  `design.md` D13, and any intentionally omitted lines from
  `fsusieR::init_prior.default` SHALL be entered in
  `inst/notes/refactor-exceptions.md`

### Requirement: empirical-Bayes mixture weight update via mixsqp

`update_model_variance.mf_individual` SHALL update mixture weights by
solving M x S independent weighted maximization problems per the
factorization in `methods/derivation.tex` line 216, using `mixsqp` by
default and a weighted-EM fallback when the user requests
`mixture_weight_method = "em"`.

#### Scenario: factorized solve

- **WHEN** the IBSS loop completes a per-effect sweep and the model
  variance update is invoked
- **THEN** the function SHALL produce one updated `pi[m, s, ]` vector
  per (modality, scale) pair, with each solve being convex and
  independent

### Requirement: cross-modality prior plug-in seam

The package SHALL accept an optional `cross_modality_prior`
argument to `mfsusie()` that, when non-NULL, modifies the per-
modality log-Bayes factors before they are summed into a joint
log-BF. The default seam implementation SHALL be
`cross_modality_prior_independent()`, a no-op that preserves the
manuscript's modality-independence assumption
(`methods/online_method.tex` line 41).

The seam dispatch verb is `combine_modality_lbfs` (an S3 generic),
NOT `apply` (avoids collision with `base::apply`). The default
implementation lives at
`combine_modality_lbfs.mf_prior_cross_modality_independent`.

#### Scenario: independent default

- **WHEN** `mfsusie()` is called without specifying
  `cross_modality_prior`, or with
  `cross_modality_prior = cross_modality_prior_independent()`
- **THEN** the joint log-BF for variant j and effect l SHALL
  equal `sum_m lbf_{l,j,m}` and the fit output SHALL match a
  per-modality independent fit run separately and post-combined

#### Scenario: contract for future implementations

- **WHEN** a non-default `cross_modality_prior` is supplied
- **THEN** an S3 method `combine_modality_lbfs(prior,
  modality_lbfs, model_state)` SHALL exist for the prior's class
  and SHALL return a length-J numeric vector for effect l; the
  IBSS SER step SHALL invoke this generic per effect per
  iteration with no other code changes in the IBSS loop
