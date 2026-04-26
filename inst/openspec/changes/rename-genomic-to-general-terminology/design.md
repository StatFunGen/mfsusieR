# Design

Hard rename, no `lifecycle::deprecate_warn()` path. Rationale:
mfsusieR is pre-1.0 with no released users. Adding a deprecation
shim doubles the API surface to maintain for zero benefit.

## Default flip

The default values for `prior_variance_scope` and the renamed
`residual_variance_scope` flip from `per_scale` to `per_outcome`.
Justification:

- fsusie paper Methods section: "We used the IS prior unless
  otherwise mentioned"; residual variance is homoskedastic per
  outcome.
- fsusieR `susiF()` `match.arg(prior)` defaults to first vector
  value `"mixture_normal"` = IS prior. Note: fsusieR roxygen has a
  doc bug claiming the default is `"mixture_normal_per_scale"`; the
  code is the source of truth.
- mvf.susie.alpha `multfsusie()` defaults to `prior =
  "mixture_normal"` = IS prior.

mfsusieR's previous `per_scale` default was a quiet divergence. The
rename undoes it.

## Independence of the two scope arguments

`prior_variance_scope` and `residual_variance_scope` are set
separately. No coupling. All four combinations are valid:

| prior        | residual     | Source                      |
|--------------|--------------|-----------------------------|
| per_outcome  | per_outcome  | fsusieR / mvf default       |
| per_outcome  | per_scale    | mfsusieR per-scale residual generalisation |
| per_scale    | per_outcome  | SPS prior (fsusie paper) + homoskedastic residual |
| per_scale    | per_scale    | mfsusieR full per-scale generalisation |

The new value vocabulary is parallel: same two strings on both
arguments.

## Naming choices

- `outcome` for Y-side index (was `modality`): matches mvSuSiE's
  canonical term in code; broader than "modality" (multi-omics
  flavor); generic for any functional-regression task.
- `variable` for X-side index (was `SNP`): matches general Bayesian
  variable selection literature; parallel with the `prior_weights`
  argument over the p variables.
- `xtx_diag` for `predictor_weights`: mirrors susieR's internal
  naming for `colSums(X^2)`. Disambiguates from `prior_weights`.
- `csd` for `csd_X`: drops the redundant `_X` suffix since the data
  class only has one column-SD field; matches mvsusieR.
- `T_basis` for `T_padded`: FDA-canonical; pairs with the kept
  `T_m` (original column count of Y[[m]]).
- `dwt_meta` for `mf_meta`: descriptive of the field's content
  (DWT pipeline metadata).
- `cross_outcome_*` for `cross_modality_*`: applies through the
  parameter, the model field, and the S3 generic; unified.
- `mixsqp_null_penalty` for the internal `nullweight` parameter:
  disambiguates from the public-API `null_prior_weight`. Two
  different concepts that previously shared a confusingly-similar
  name.

## What stays

- Index variables: `m` (outcome index), `j` (variable index), `l`
  (effect index), `s` (scale index), `t` (position index). These
  match the manuscript subscripts exactly.
- `mf_*` prefix on internal helpers: the package is mfsusieR.
- `M` (number of outcomes), `p` (number of variables), `L` (number
  of effects), `T_m` (positions per outcome), `T_basis[m]` (basis
  size for outcome m): scalar dimensions.
- `mf_individual` data class name: changing it would cascade
  through the susieR S3 dispatch registration in `.onLoad`. Not
  worth it.
