# mvsusieR (refactor-s3) S3 paradigm

This is paradigm reference #1 for mfsusieR: how to extend `susieR` to handle a
non-scalar response (multi-trait Y) by registering S3 methods on new data and
model classes, while leaving susieR's IBSS loop untouched.

Source root: `~/GIT/mvsusieR` on branch `refactor-s3`. All citations are
relative to that root.

## 1. Package layout

The R/ directory does not mirror susieR file names. New files, organized by
data path:

- `mvsusie.R` - public API (`mvsusie`, `mvsusie_ss`, `mvsusie_rss`) and
  internal `mvsusie_core`, `mvsusie_workhorse`.
- `mvsusie_constructors.R` - data and model constructors and the matching
  `initialize_susie_model.*` methods.
- `individual_data_methods.R` (~1290 lines) - S3 methods for the
  `mv_individual` data class.
- `sufficient_stats_methods.R` - overrides for the `mv_ss` data class.
- `mixture_prior.R` - `create_mixture_prior`, `create_mash_prior`, canonical
  covariance constructors.
- `model_methods.R` - methods on the `mvsusie` model class
  (`get_prior_variance_l.mvsusie`, etc.).
- `mvsusie_get_functions.R` - output formatting and LFSR.
- `predict.mvsusie.R` - `coef.mvsusie`, `predict.mvsusie`.
- `missing_y_*.R` - 3D missing-data path.
- `zzz.R` - `.onLoad` registers all S3 methods in `susieR`'s namespace and
  caches `mashr` and `susieR` internals.

`mvsusie()` is at `R/mvsusie.R:349`. `mvsusie_ss()` at `R/mvsusie.R:523`,
`mvsusie_rss()` at `R/mvsusie.R:431`. Internal orchestrator
`mvsusie_core()` at `R/mvsusie.R:674`. The actual loop runs inside
`susieR::susie_workhorse()` (called at `R/mvsusie.R:93`).

## 2. Class hierarchy (S3)

Three new classes:

- `mv_individual` (`R/mvsusie_constructors.R:151`) inherits from `individual`.
  Created by `create_mvsusie_data()`. Carries Y (n by R), residual covariance,
  per-pattern caches, missing-data structures.
- `mv_ss` (`R/mvsusie_constructors.R:331`) inherits from `mv_individual` and
  `ss`. Created by `create_mvsusie_ss_data()`. Adds `XtX`, `XtY`, `YtY`.
- `mvsusie` (`R/mvsusie_constructors.R:425`) inherits from `susie`. Constructed
  by `initialize_susie_model.mv_individual()` (`R/mvsusie_constructors.R:341`)
  or `.mv_ss()`. Stores `alpha` (L by J), `mu` (L by J by R for R > 1),
  `mu2_cache`, `V_structure` (list of K covariance matrices), `pi_V` (mixture
  weights), `pi_V_posterior` (per-effect, per-variable, per-component
  posteriors), `lbf_outcome`, `conditional_lfsr`, plus the standard susie
  fields.

A fourth, `mash_prior` (`R/mixture_prior.R:244`), is the prior object: a list
of K RxR covariance matrices, mixture weights, and a null weight.

Public S3 methods on the fit:

- `coef.mvsusie` (`R/predict.mvsusie.R:20`) - returns a (J+1) by R coefficient
  matrix.
- `predict.mvsusie` (`R/predict.mvsusie.R:68`) - applies coefficients to new X.

Internal S3 dispatch is done by registering methods on the data classes inside
`zzz.R:62-111`. The methods registered for `mv_individual` cover the full
extension surface: `ibss_initialize`, `SER_posterior_e_loglik`,
`calculate_posterior_moments`, `cleanup_model`, `compute_kl`,
`compute_residuals`, `compute_ser_statistics`, `em_update_prior_variance`,
`get_cs`, `get_fitted`, `get_intercept`, `get_objective`, `get_scale_factors`,
`get_var_y`, `get_variable_names`, `get_zscore`, `initialize_fitted`,
`initialize_susie_model`, `loglik`, `neg_loglik`, `trim_null_effects`,
`update_fitted_values`, `update_model_variance`, `update_variance_components`.
The `mv_ss` registrations override the residual-related and intercept-related
ones for the sufficient-statistics path. Methods on the `mvsusie` model
(`R/model_methods.R:2-33`, registered at `zzz.R:114-125`) cover
`get_prior_variance_l`, `set_prior_variance_l`, `get_alpha_l`,
`get_posterior_mean_l`, `get_posterior_mean_sum`, `get_posterior_moments_l`.

## 3. Delegation to susieR

The IBSS loop is reused, not reimplemented. `mvsusie_workhorse` builds the
`params` list and a multivariate `data` object, then calls
`susieR::susie_workhorse(data, params)` (`R/mvsusie.R:93`). Every per-iteration
quantity is reached by S3 dispatch on the data class.

`zzz.R:42-57` caches a small set of `mashr` and `susieR` internals so
mvsusieR can call them without `:::` at runtime: `mashr::calc_lik_rcpp`,
`mashr::calc_sermix_rcpp`, `mashr::compute_null_loglik_from_matrix`,
`mashr::compute_alt_loglik_from_matrix_and_pi`, `mashr::expand_cov`,
`susieR::get_var_y`, `susieR::initialize_susie_model`,
`susieR::initialize_fitted`, `susieR::SER_posterior_e_loglik`,
`susieR::compute_kl`, `susieR::warning_message`, `susieR::mem_used_gb`,
`susieR::check_semi_pd`. `susie_get_cs` is called explicitly at
`R/mvsusie.R:911` and `R/mvsusie.R:1053` to build credible sets.
`ashr::compute_lfsr` is imported in `R/mvsusie_get_functions.R:2`.

There is no major reimplemented routine. The multivariate variants of
residual variance estimation, prior variance estimation, and ELBO computation
are S3 methods (`estimate_residual_variance_mv*`,
`em_update_prior_variance.mv_individual`, `compute_multivariate_elbo*`); they
are dispatched into susieR's loop, not running outside it.

## 4. IBSS loop differences vs susieR

The shape of one outer iteration is unchanged: residuals, SER stats, posterior
moments, KL, fitted update, variance component update, ELBO, convergence check.
What changes is the math inside each generic.

Multivariate additions, all entered via S3 methods on `mv_individual` /
`mv_ss`:

- Residual covariance is an R by R matrix. `update_variance_components.mv_individual`
  estimates it from the expected sufficient statistics
  ((Y'Y - 2 B' X'Y + B' X'X B) / N), and stores the inverse and the spectral
  factors used in subsequent SER calls.
- The mixture prior plugs in inside `compute_ser_statistics.mv_individual` /
  `.mv_ss`. When `model$V_structure` is a list of K matrices and
  `model$pi_V` is set, the method calls `mashr::calc_sermix_rcpp`; otherwise it
  uses `mashr::calc_lik_rcpp` with a single covariance.
- Posterior mixture weights `pi_V_posterior` are L by J by K (or LxJx(K+1)
  with a null component). They are written by the SER routine each iteration.
- After the L sweep, mixture weights `pi_V` themselves are updated via
  mixsqp inside `update_model_variance.mv_individual` (the function
  `update_mixture_weights_and_prune` does the EM/mixsqp step and prunes
  vanishing components).
- For missing Y, a 3D path activates: per missing-data pattern, a separate
  precision matrix is used in a GLS-like SER. When residual variance estimation
  is on and missing data are present, the workhorse switches to a block
  coordinate ascent over (alpha, mu, V) versus sigma2 (`R/mvsusie.R:835-886`).
- A per-outcome filter computes `lbf_outcome` and `conditional_lfsr` so that
  outcomes whose conditional log-BF is below `min_outcome_lbf` are not
  borrowed across by a multi-trait effect.

## 5. Prior object interface

`mash_prior` is the canonical mvsusieR prior (`R/mixture_prior.R:244-320`).
Fields: `xUlist` (K covariance matrices, scaled or unscaled), `pi` (mixture
weights), `null_weight`, `usepointmass`, `R`. Constructed three ways:

- `create_mash_prior(Ulist, grid)` - takes mashr unscaled bases and a scaling
  grid.
- `create_mixture_prior(fitted_g | mixture_prior | R)` - higher-level wrapper
  (`R/mixture_prior.R:66`); accepts a fitted mash object, a user list, or just
  the number of outcomes (then generates canonical singletons +
  heterogeneity).
- Direct construction with a precomputed `xUlist`.

Generics on the prior: `scale_prior_variance.mash_prior`
(`R/mvsusie.R:728`) for runtime rescaling, plus the constructor and the
mashr-side `expand_cov`. The dispatch into the SER lives at the data-class
boundary, not the prior boundary: the methods on `mv_individual` /
`mv_ss` look at `model$V_structure` and decide which mashr Rcpp routine to
call. A third party could supply a different prior class, but they would
need their own SER-stats method that knows how to consume it; the current
design is not abstracted at the prior-class level.

## 6. Public API signature

`mvsusie()` (`R/mvsusie.R:349`):

```
mvsusie(X, Y, L = 10, prior_variance = 0.2,
        residual_variance = NULL, prior_weights = NULL,
        standardize = TRUE, intercept = TRUE,
        estimate_residual_variance = TRUE,
        estimate_prior_variance = TRUE,
        estimate_prior_method = "optim",
        estimate_prior_mixture_weights = TRUE,
        mixture_weight_method = "mixsqp",
        check_null_threshold = 0, prior_tol = 1e-9,
        model_init = NULL,
        missing_y_method = "approximate",
        coverage = 0.95, min_abs_corr = 0.5,
        compute_univariate_zscore = FALSE,
        precompute_cache = TRUE, n_thread = 1,
        max_iter = 100, tol = 1e-3, verbose = TRUE,
        track_fit = FALSE,
        min_outcome_lbf = 0)
```

Diff vs `susieR::susie()`: adds `estimate_prior_mixture_weights`,
`mixture_weight_method`, `missing_y_method`, `n_thread`, `precompute_cache`,
`min_outcome_lbf`, `compute_univariate_zscore`. `prior_variance` accepts a
scalar (R = 1), an R by R matrix, or a `mash_prior` object. Drops
`null_weight` (folded into the prior), `refine`, `unmappable_effects`,
`residual_variance_lowerbound`/`upperbound`, `n_purity`, `convergence_method`.
Default `tol` is 1e-3 instead of 1e-4.

## 7. Return object

`format_mvsusie_output()` (`R/mvsusie_get_functions.R:16-141`) builds the final
list. Differences from a plain susie fit:

- `mu` and `mu2_diag` are L by J by R arrays for R > 1 (matrices for R = 1).
- `sigma2` is an R by R matrix for R > 1 (scalar for R = 1).
- New: `V_structure`, `null_weight`, `prior_mixture_weights`,
  `posterior_mixture_weights` (L by J by K or by K+1),
  `lbf_outcome` (L by R), `conditional_lfsr` (L by J by R),
  `single_effect_lfsr` (L by R), `lfsr` (J by R, collapsed),
  `outcome_names`, `variable_names`, `X_column_scale_factors`.
- Standard fields preserved: `alpha`, `KL`, `lbf`, `lbf_variable`, `pi`,
  `Xr`, `elbo`, `niter`, `converged`, `fitted`, `intercept`, `sets`,
  `pip`, `trace`.

## 8. Evidence pointers

- `R/mvsusie.R:349`, `R/mvsusie.R:523`, `R/mvsusie.R:431` - public entries.
- `R/mvsusie.R:674`, `R/mvsusie.R:93` - core / workhorse.
- `R/mvsusie_constructors.R:21`, `:271`, `:341`, `:425` - data and model
  constructors.
- `R/mixture_prior.R:66`, `:244` - prior constructors.
- `R/mvsusie.R:728` - `scale_prior_variance.mash_prior`.
- `R/zzz.R:37-126` - method registration into susieR's namespace.
- `R/individual_data_methods.R:9` (`get_var_y`),
  `~:400-600` (SER stats and dispatch into mashr Rcpp),
  `~:650-770` (residual variance estimation),
  `~:950-1050` (mixture weight update),
  `~:1100-1200` (`em_update_prior_variance`).
- `R/sufficient_stats_methods.R:5-150` - SS-path overrides.
- `R/mvsusie_get_functions.R:16-141` - output formatting.
- `R/mvsusie.R:835-886` - block coordinate ascent for missing-Y +
  residual-variance estimation.
- `R/predict.mvsusie.R:20`, `:68` - coef and predict.

## Implications for mfsusieR

mvsusieR is the cleanest available demonstration that susieR's IBSS loop can
be reused unchanged when a new response structure is added through S3
dispatch on the data class plus a registered model class. mfsusieR can adopt
the same structure: a `multifunctional_individual` (and matching `*_ss`)
data class and an `mfsusie` model class inheriting from `susie`, with method
registration in `zzz.R`. The mixture-prior dispatch shows how to compose a
data class with a prior structure that sits in `model$V_structure` and is
read by `compute_ser_statistics.*`; this is the mechanism mfsusieR will use
to plug in scale-specific wavelet priors. The block coordinate ascent for
missing-Y + residual variance is a precedent for any iteration shape we
need beyond plain IBSS, including a possible per-scale residual variance
update that the FDR investigation may require.
