# Rename: genomic-specific terminology to functional-regression / Bayesian-VS general terms

## Why

mfsusieR was written with genomic-fine-mapping terminology baked into
both the public API and internal source: "SNP" for X columns,
"modality" for Y list elements, "trait" in occasional prose. The
package itself is a general method for Bayesian variable selection in
functional regression with mixed functional and scalar responses.
Genomic naming narrows the audience and adds friction for non-genetics
applications.

mvSuSiE uses "outcome" for Y-side indexing. fsusie paper uses
"variable" for X-side indexing. Adopting these aligns mfsusieR with
the SuSiE family while making the API readable for any
functional-regression user.

A second issue caught during this audit: the public-API value names
`per_scale_modality`, `per_modality`, `shared_per_modality` are
internally inconsistent (`shared_per_modality` reads as "shared OR
per?"). The semantic axis is just per-scale vs per-outcome variance.

A third issue: mfsusieR currently defaults to `per_scale_modality` for
both `prior_variance_scope` and `residual_variance_method`. fsusieR
and mvf.susie.alpha both default to single per-outcome (the IS prior
and homoskedastic residual, respectively). The fsusie paper Methods
section is explicit: "We used the IS prior unless otherwise mentioned"
and "we make the additional modeling assumption that the residual
variances do not depend on location". mfsusieR's default should match.

## What changes

Public API (hard rename, pre-1.0, no users affected):
- Argument `cross_modality_prior` -> `cross_outcome_prior`
- Argument `prior_variance_scope` value vocabulary
  `c("per_scale_modality", "per_modality")` -> `c("per_outcome", "per_scale")`,
  default flips to `"per_outcome"` to match fsusieR / mvf / paper
- Argument `residual_variance_method` -> `residual_variance_scope`,
  values `c("per_scale_modality", "shared_per_modality")` ->
  `c("per_outcome", "per_scale")`, default flips to `"per_outcome"`

Internal data and model fields:
- `data$predictor_weights` -> `data$xtx_diag` (susieR convention)
- `data$csd_X` -> `data$csd` (mvsusieR convention)
- `data$T_padded` -> `data$T_basis` (FDA convention; pairs with `T_m`)
- `model$cross_modality_combiner` -> `model$cross_outcome_combiner`
- `fit$mf_meta` -> `fit$dwt_meta`

Internal helpers and S3 generics:
- `combine_modality_lbfs.*` -> `combine_outcome_lbfs.*`
- `cross_modality_prior_independent()` -> `cross_outcome_prior_independent()`
- `mf_per_modality_bhat_shat` -> `mf_per_outcome_bhat_shat`
- `mf_invert_per_modality` -> `mf_invert_per_outcome`
- Internal mixsqp penalty `nullweight` -> `mixsqp_null_penalty`
  (disambiguates from public `null_prior_weight` mixture-component weight)

Prose (roxygen, comments, vignettes):
- "SNP", "SNPs" -> "variable", "variables"
- "modality", "modalities" -> "outcome", "outcomes"
- "trait", "traits" -> "outcome", "outcomes"

## Impact

- Affected specs: `mf-public-api`, `mf-prior` (delta).
- Affected code: every R/ source file, every test, every vignette.
- Affected tests: 3 test files set the legacy `shared_per_modality`
  / `per_scale_modality` value names explicitly; these update.
  Default-flip means C1/C2/C3 contracts now run with the correct
  per-outcome upstream-matching default automatically.
- No deprecation path: pre-1.0, no external users.
