## ADDED Requirements

### Requirement: public mfsusie() entry function

The package SHALL export `mfsusie(X, Y, pos = NULL, L = 10, ...)` as
the multi-modality entry for individual-data fits, with arguments
named per CLAUDE.md naming rules (snake_case, no abbreviations
beyond `pip`, `cs`, `lbf`, `elbo`, `kl`).

#### Scenario: signature matches design

- **WHEN** `args(mfsusie)` is inspected
- **THEN** the argument list SHALL include `X`, `Y`, `pos`, `L`,
  `L_greedy`, `prior_variance_scope`,
  `prior_variance_grid_multiplier`, `prior_variance_grid`,
  `null_prior_weight`, `cross_modality_prior`,
  `residual_variance_method`, `estimate_residual_variance`,
  `estimate_prior_variance`, `estimate_prior_method`,
  `estimate_prior_mixture_weights`, `mixture_weight_method`,
  `check_null_threshold`, `prior_tol`, `max_padded_log2`,
  `max_iter`, `tol`, `coverage`, `min_abs_corr`,
  `filter_credible_sets`, `wavelet_filter_number`,
  `wavelet_family`, `standardize`, `intercept`, `precompute_cache`,
  `save_residuals`, `n_thread`, `verbose`, `track_fit`, `model_init`

#### Scenario: forbidden names absent from mfsusie()

- **WHEN** `formalArgs(mfsusie)` is inspected
- **THEN** none of the following names from `mvf.susie.alpha::multfsusie`
  or `fsusieR::susiF` SHALL appear: `nullweight`, `gridmult`,
  `max_scale`, `max_SNP_EM`, `max_step_EM`, `max_step`, `cal_obj`,
  `cov_lev`, `min_purity`, `filter_cs`, `filter.number`, `family`,
  `init_pi0_w`, `tol_null_prior`, `lbf_min`, `posthoc`, `cor_small`,
  `thresh_lowcount`, `greedy`, `backfit`, `multfsusie.obj`,
  `quantile_trans`, `post_processing`. A test in
  `tests/testthat/test_public_api_naming.R` SHALL assert this
  programmatically.

#### Scenario: input validation

- **WHEN** `mfsusie()` is called with `Y` that is not a list of
  numeric matrices, with `length(Y) != length(pos)`, with
  `nrow(Y[[m]]) != nrow(X)` for any m, or with
  `length(pos[[m]]) != ncol(Y[[m]])` for any m
- **THEN** the function SHALL stop with an informative error message
  identifying the offending argument

### Requirement: public fsusie() thin wrapper

The package SHALL export `fsusie(Y, X, pos = NULL, L = 10, ...)` as
a thin wrapper that accepts a single phenotype as input and forwards
to `mfsusie()` with `M = 1`. Argument order matches `fsusieR::susiF`
`(Y, X, ...)` for drop-in migration.

#### Scenario: fsusie() canonicalizes single-phenotype input

- **WHEN** `fsusie(Y, X, pos)` is called with `Y` either a numeric
  matrix `n x T` or a numeric vector of length `n`, and `X` an
  `n x J` numeric matrix
- **THEN** the wrapper SHALL canonicalize the input as
  `Y_canonical = list(Y_or_matrix)`,
  `pos_canonical = list(pos %||% seq_len(T))`,
  forward to `mfsusie(X, Y_canonical, pos_canonical, L, ...)`,
  and return the resulting `c("mfsusie", "susie")` fit object

#### Scenario: fsusie() rejects multi-modality-only arguments

- **WHEN** `fsusie()` is called with an argument that is not
  meaningful for the `M = 1` case (e.g., `cross_modality_prior`
  set to a non-NULL value)
- **THEN** the wrapper SHALL stop with an informative error
  message naming the offending argument

#### Scenario: fsusie() drop-in compatibility

- **WHEN** `fsusie(Y_matrix, X, pos)` and
  `mfsusie(X, list(Y_matrix), list(pos))` are called with the
  same fixed seed
- **THEN** the two fits SHALL be element-wise identical at
  tolerance `1e-12` (the wrapper is a canonicalizer; no numerical
  changes)

### Requirement: forbidden names absent from fsusie()

The `fsusie()` thin wrapper SHALL NOT propagate any of the legacy
non-snake_case argument names from `fsusieR::susiF` or
`mvf.susie.alpha::multfsusie` into its public signature.

#### Scenario: forbidden names absent from fsusie() formalArgs

- **WHEN** `formalArgs(fsusie)` is inspected
- **THEN** none of the names listed in the `mfsusie()` forbidden
  list above (e.g., `nullweight`, `gridmult`, `max_scale`,
  `max_SNP_EM`, `cov_lev`, `min_purity`, `filter_cs`,
  `filter.number`, `family`, `init_pi0_w`, `tol_null_prior`,
  `lbf_min`, `posthoc`, `cor_small`, `thresh_lowcount`, `greedy`,
  `backfit`, `multfsusie.obj`, `quantile_trans`, `post_processing`)
  SHALL appear in the `fsusie()` signature

### Requirement: fit-object class hierarchy

`mfsusie()` and `fsusie()` SHALL each return an object of class
`c("mfsusie", "susie")` whose fields conform to the data-class spec.

#### Scenario: class and fields

- **WHEN** the user inspects `class(fit)` and `names(fit)` for a
  successful run
- **THEN** `class(fit)` SHALL be `c("mfsusie", "susie")`, and
  `names(fit)` SHALL include at minimum `alpha`, `mu`, `mu2`, `KL`,
  `lbf`, `lbf_variable`, `pi`, `pi_V`, `V_grid`,
  `null_prior_weight`, `sigma2`, `elbo`, `niter`, `converged`,
  `fitted`, `residuals` (NULL when `save_residuals = FALSE`),
  `intercept`, `pip`, `cs`, `csd_X`, `wavelet_meta`, `call`

### Requirement: standard S3 methods on the fit

The package SHALL provide `coef.mfsusie`, `predict.mfsusie`,
`fitted.mfsusie`, `summary.mfsusie`, `print.mfsusie`. Inherited
methods from `susie` (e.g., `susieR::susie_get_pip`,
`susieR::susie_get_cs`) SHALL work where their semantics match the
multi-modality functional setting.

#### Scenario: coef returns per-effect curves

- **WHEN** `coef(fit)` is called
- **THEN** the function SHALL return a list of length M, each entry
  an `L x T_m` numeric matrix of per-effect curves in the
  measurement space, computed by inverse DWT of the alpha-weighted
  posterior-mean wavelet coefficients divided by the cached
  `fit$csd_X`. The function SHALL NOT require X or Y as inputs.

#### Scenario: predict returns reconstructed curves on new data

- **WHEN** `predict(fit, newx)` is called with a valid `newx`
  matrix
- **THEN** the function SHALL return a list of M numeric matrices,
  each of dimension `nrow(newx) x T_m`, computed by inverse DWT of
  the posterior-mean wavelet coefficients summed across effects

#### Scenario: fitted returns sum-of-effects on training data

- **WHEN** `fitted(fit)` is called
- **THEN** the function SHALL return a list of M `n x T_m` numeric
  matrices read from `fit$fitted` (which was filled at IBSS
  finalize time); the function SHALL NOT require X as input

### Requirement: defaults are documented and stable

Default argument values SHALL be documented in roxygen with a
citation to the source (CLAUDE.md naming rule, manuscript reference,
`mvf.susie.alpha` legacy value, or `fsusieR::susiF` legacy value).
Defaults SHALL NOT change between patch releases.

#### Scenario: roxygen documents defaults

- **WHEN** `?mfsusie` and `?fsusie` are rendered after
  `roxygen2::roxygenise()`
- **THEN** every argument's documentation block SHALL state the
  default value and a one-line rationale or citation

#### Scenario: post-processing rationale documented

- **WHEN** `?mfsusie` is rendered
- **THEN** the function-level documentation SHALL include a
  short paragraph in statgen-writing-style prose explaining that
  smoothing and credible bands are NOT fit-time arguments but are
  separate post-fit transforms (see `?mf_post_smooth`,
  `?mf_credible_bands`), and explaining the residual-storage
  contract (`save_residuals = TRUE` by default)
