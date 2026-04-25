## ADDED Requirements

### Requirement: public mfsusie() entry function

The package SHALL export `mfsusie(X, Y, pos = NULL, L = 10, ...)` as the
sole supported entry for individual-data fits, with arguments named per
CLAUDE.md naming rules (snake_case, no abbreviations beyond
`pip`, `cs`, `lbf`, `elbo`, `kl`).

#### Scenario: signature matches design

- **WHEN** `args(mfsusie)` is inspected
- **THEN** the argument list SHALL include `X`, `Y`, `pos`, `L`,
  `prior_variance_scope`, `prior_variance_grid_multiplier`,
  `prior_variance_grid`, `null_prior_weight`, `cross_modality_prior`,
  `residual_variance_method`, `estimate_residual_variance`,
  `estimate_prior_mixture_weights`, `mixture_weight_method`,
  `max_padded_log2`, `max_iter`, `tol`, `coverage`, `min_abs_corr`,
  `filter_credible_sets`, `wavelet_filter_number`, `wavelet_family`,
  `standardize`, `intercept`, `verbose`, `track_fit`, `model_init`,
  with no abbreviations from `mvf.susie.alpha::multfsusie` retained
  (`nullweight`, `gridmult`, `max_SNP_EM`, `max_scale`, `cal_obj`,
  `cov_lev`, `min_purity`, `filter_cs`, `filter.number`, `family` are
  forbidden in the public signature)

#### Scenario: input validation

- **WHEN** `mfsusie()` is called with `Y` that is not a list of numeric
  matrices, with `length(Y) != length(pos)`, with `nrow(Y[[m]]) !=
  nrow(X)` for any m, or with `length(pos[[m]]) != ncol(Y[[m]])` for any
  m
- **THEN** the function SHALL stop with an informative error message
  identifying the offending argument

### Requirement: fit-object class hierarchy

`mfsusie()` SHALL return an object of class `c("mfsusie", "susie")`
whose fields conform to the data-class spec.

#### Scenario: class and fields

- **WHEN** the user inspects `class(fit)` and `names(fit)` for a
  successful run
- **THEN** `class(fit)` SHALL be `c("mfsusie", "susie")`, and
  `names(fit)` SHALL include at minimum `alpha`, `mu`, `mu2`, `KL`,
  `lbf`, `lbf_variable`, `pi`, `pi_V`, `V_grid`, `null_prior_weight`,
  `sigma2`, `elbo`, `niter`, `converged`, `fitted`, `intercept`, `pip`,
  `cs`, `wavelet_meta`, `call`

### Requirement: standard S3 methods on the fit

The package SHALL provide `coef.mfsusie`, `predict.mfsusie`,
`summary.mfsusie`, `print.mfsusie`. Inherited methods from `susie` (e.g.
`susieR::susie_get_pip`, `susieR::susie_get_cs`) SHALL work where their
semantics match the multi-modality functional setting.

#### Scenario: predict returns reconstructed curves

- **WHEN** `predict(fit, newx)` is called with a valid `newx`
- **THEN** the function SHALL return a list of M numeric matrices, each
  of dimension `nrow(newx) x T_m`, computed by inverse DWT of the
  posterior-mean wavelet coefficients

#### Scenario: coef returns wavelet-space coefficients

- **WHEN** `coef(fit, space = "measurement")` is called
- **THEN** the function SHALL return a list of M numeric matrices, each
  of dimension `J x T_m`, in the original measurement space; with
  `space = "wavelet"` the returned matrices SHALL be in the
  wavelet-coefficient basis with dimension `J x T_padded[m]`

### Requirement: defaults are documented and stable

Default argument values SHALL be documented in roxygen with a citation
to the source (CLAUDE.md naming rule, manuscript reference, or
`mvf.susie.alpha` legacy value). Defaults SHALL NOT change between
patch releases.

#### Scenario: roxygen documents defaults

- **WHEN** `?mfsusie` is rendered after `roxygen2::roxygenise()`
- **THEN** every argument's documentation block SHALL state the default
  value and a one-line rationale or citation
