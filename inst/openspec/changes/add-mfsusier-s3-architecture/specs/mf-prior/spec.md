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
  `V_grid` as a length-M list of length-K vectors, replicated per scale
  during SER computation, matching `mvf.susie.alpha`'s per-trait
  treatment

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

The package SHALL accept an optional `cross_modality_prior` argument to
`mfsusie()` that, when non-NULL, modifies the per-modality log-Bayes
factors before they are summed into a joint log-BF. The default seam
implementation SHALL be `cross_modality_prior_independent()`, a no-op
that preserves the manuscript's modality-independence assumption
(`methods/online_method.tex` line 41).

#### Scenario: independent default

- **WHEN** `mfsusie()` is called without specifying
  `cross_modality_prior`, or with
  `cross_modality_prior = cross_modality_prior_independent()`
- **THEN** the joint log-BF for variant j and effect l SHALL equal
  `sum_m lbf_{l,j,m}` and the fit output SHALL match a per-modality
  independent fit run separately and post-combined

#### Scenario: contract for future implementations

- **WHEN** a non-default `cross_modality_prior` is supplied
- **THEN** the prior object SHALL implement `apply(modality_lbfs,
  model_state) -> joint_lbf` returning a length-J numeric vector for
  effect l, and the IBSS SER step SHALL invoke this method per effect
  per iteration with no other code changes in the IBSS loop
