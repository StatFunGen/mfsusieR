## Context

mfsusieR is one S3 package on top of the `susieR` backbone that
replaces both `mvf.susie.alpha::multfsusie` (multi-modality functional
SuSiE) and `fsusieR::susiF` (single-modality functional SuSiE) by
making each a special case of one model. Architecture follows the
S3 conventions of `mvsusieR/refactor-s3` (paradigm reference #1) and
`susieAnn` (paradigm reference #2). The math anchor is the manuscript
at `/home/gw/GIT/MultifSuSiE_Manuscript`. Phase 1 paradigm notes
(`inst/notes/paradigms/`) catalogued the four packages and identified
six candidate sources for the FDR miscalibration Anjing observed (see
`mvf-original.md` section 8). This design fixes the architecture so
Phase 3 implementation has a stable target and Phase 4 reference
tests have clear contracts.

The model: M modalities, each with N samples and T_m positions.
Wavelet transform per modality maps Y_m to D_m of size N x T_m. Per
modality, fit a SuSiE-style sum-of-single-effects on D_m with a
per-(scale, modality) scale-mixture-of-normals prior. Bayes factors
combine across modalities by product (modality-independence
assumption, manuscript online_method line 41).

Three apple-to-apple equivalence contracts hold the port together:

| Contract | Inputs | Reference | Tolerance |
|---|---|---|---|
| C1: scalar SuSiE | `M = 1`, `T_1 = 1`, single-component prior | `susieR::susie` | `<= 1e-10` |
| C2: single-modality functional SuSiE | `M = 1`, `T_1 > 1` | `fsusieR::susiF` | `<= 1e-8` |
| C3: multi-modality functional SuSiE | `M >= 1`, legacy variance mode | `mvf.susie.alpha::multfsusie` | `<= 1e-8` |

When a port source has a known bug, mfsusieR fixes it, the
corresponding test asserts the deviation explicitly, and the
OpenSpec change that authorized the deviation is cited inline.

## Goals / Non-Goals

**Goals:**

- Replace `mvf.susie.alpha::multfsusie` AND `fsusieR::susiF` with one
  S3 implementation delegating to `susieR::susie_workhorse` for the
  IBSS loop.
- One public surface: `mfsusieR::mfsusie()` for multi-modality and
  `mfsusieR::fsusie()` for single-modality, the latter being a thin
  wrapper that calls the former with `M = 1`.
- Storage shape: per-(scale, modality) prior variance and mixture
  weights natively supported in `mf_individual` and the `mfsusie`
  model object.
- Single data class collapses degenerate cases: `M = 1`, `T_m = 1`,
  and ragged `T_m` all flow through the same code.
- Public API follows CLAUDE.md naming rules. No abbreviations beyond
  `pip`, `cs`, `lbf`, `elbo`, `kl`. Snake_case throughout.
- Pluggable seam for cross-modality covariance prior. Default fits
  modalities independently per the manuscript; a future change
  provides a mash-style implementation.
- PIPs are a deterministic function of the *final* alpha state. CS
  filtering happens before PIP computation, fixing the
  Phase-1-identified PIP-after-filter bug.
- Per-(scale, modality) residual variance is the v1 default. Legacy
  shared-per-modality mode preserved as `residual_variance_method =
  "shared_per_modality"` for Phase 4 contract C3 fidelity tests.
- Decoupled post-processing. Smoothing and credible-band computation
  are functions on the fit, not arguments to `mfsusie()`. The fit
  stores per-modality residuals by default so post-processors take
  only the fit.

**Non-Goals:**

- Sufficient-statistics path. Will be a follow-up change once the
  individual-data path stabilizes.
- Concrete cross-modality covariance prior. v1 ships only the seam;
  the mash-style implementation is a separate change.
- FDR investigation (Phase 5) and the fixes that follow it (Phase 6).
  This change makes those investigations possible without dictating
  their outcome.
- Wholesale Rcpp / RcppArmadillo migration. Phase 7 only, after the
  architecture is locked. Hot-path C++ via `cpp11` is permitted
  during Phase 3 for narrow utility kernels (per design.md D14)
  with a pure-R reference oracle and an apple-to-apple
  C++-vs-pure-R unit test at `<= 1e-12`.
- Deprecation of the read-only reference packages. `susieR`,
  `mvsusieR`, `susieAnn`, `mvf.susie.alpha`, and `fsusieR` remain
  read-only references for the lifetime of mfsusieR, per CLAUDE.md
  hard rule #1.
- `fsusieR::EBmvFR`. A different model (EB multivariate functional
  regression with no SuSiE structure). Out of scope for v1.

## Decisions

### D1. Single data class `mf_individual` covers all degenerate cases

`mvsusieR/refactor-s3` collapses `R = 1` vs `R > 1` inside one
`mv_individual` class via shape-aware branches. mfsusieR follows the
same pattern: `mf_individual` is the only data class in v1, and `M =
1`, `T_m = 1`, ragged T_m across modalities are all handled by shape
branches. `mf_ss` is deferred.

Alternatives considered:

- One class per modality count and one per per-position count:
  rejected, combinatorial explosion.
- Two classes (`mf_individual` for M > 1 functional, `mf_individual_uni`
  for univariate-only): rejected, duplicate method registrations and
  duplicate test surface.

### D2. Storage shape

`mf_individual` carries:

- `Y` as `list[M]` of `n x T_m` matrices (ragged T_m allowed).
- `pos` as `list[M]` of position vectors of length `T_m`.
- `D` as `list[M]` of `n x T_m` wavelet-coefficient matrices
  (post-DWT, post-padding, post-column-scaling). Cached at
  construction.
- `scale_index` as `list[M]` of integer vectors mapping each column of
  `D[[m]]` to a scale level (output of the ported
  `gen_wavelet_indx(log2(T_m))` helper in `R/utils_wavelet.R`).
- `T_padded` as `integer[M]` recording each modality's padded length
  (next power of 2 above original `T_m`).
- `csd_X` as `numeric[J]`. Per-SNP column scale of `X` cached at
  construction. Used by `coef.mfsusie` and the post-processors when
  reconstructing per-effect curves.
- `residuals` as `list[M]` of `n x T_m` matrices when `save_residuals
  = TRUE` (default). Filled at the end of the IBSS loop with `Y_m -
  X * sum_l (alpha_l . mu_l)` projected back to the measurement
  space. NULL when `save_residuals = FALSE`.
- `X`, `n`, `J` standard.

The `mfsusie` model object adds:

- `alpha` as `L x J` matrix (same as susieR).
- `mu`, `mu2` as `list[L]`, each entry a `list[M]`, each inner entry
  a `J x T_padded[m]` matrix. The nested list lets ragged T_m
  coexist; we do not flatten into a 3D array.
- `pi_V` as `list[M]` of `S_m x K` matrices (`pi_{k,s,m}` per
  manuscript). `S_m = log2(T_padded[m])`.
- `V_grid` as `list[M]` of length-K grids of prior variances
  `sigma_{k,s,m}^2`. Stored once per modality; per-scale variance is
  optional and stored as `list[M]` of `S_m x K` matrices when used.
- `sigma2` as either a scalar per modality (legacy mode) or a
  `list[M]` of length-`S_m` vectors (per-scale-per-modality default).
  A single field holds whichever shape is active; the dispatch on
  `residual_variance_method` is at update-component time.
- `lbf`, `lbf_variable`, `KL`, `pip` standard susieR shapes.

### D3. Class hierarchy

```
mf_individual            (data class; inherits from `individual`)
  -> dispatched into susieR::susie_workhorse via S3 methods

mfsusie                  (model class; inherits from `susie`)
  -> coef.mfsusie, predict.mfsusie, fitted.mfsusie,
     summary.mfsusie, print.mfsusie

mf_prior_scale_mixture   (default prior)
  -> per-(scale, modality) scale-mixture-of-normals
mf_prior_cross_modality  (plug-in seam, no default impl in v1)
  -> mash-style covariance across modalities, optional layer
```

### D4. S3 dispatch surface

mfsusieR registers methods on `mf_individual` for every generic that
mvsusieR overrides on `mv_individual`. We do not invent new generics.
The methods registered into `susieR`'s namespace via `zzz.R::.onLoad`:

- `initialize_susie_model.mf_individual` - construct the `mfsusie`
  model with the per-(scale, modality) shapes from D2.
- `compute_residuals.mf_individual` - residuals across (modality,
  scale, position).
- `compute_ser_statistics.mf_individual` - per modality, compute
  per-(scale, position) `Bhat`, `Shat`. Combines per-modality log-BFs
  by sum (manuscript algorithms.tex line 20: product on the BF scale,
  sum on the log-BF scale).
- `calculate_posterior_moments.mf_individual` - mixture posterior over
  K components, per (modality, scale, position). Manuscript
  eq:post_f_mix, eq:post_f2_mix.
- `compute_kl.mf_individual` - KL divergence for the mixture prior.
- `loglik.mf_individual`, `neg_loglik.mf_individual` - marginal
  log-likelihood used by the EM step on mixture weights.
- `update_fitted_values.mf_individual` - store the contribution of
  effect l back into the running fit.
- `update_variance_components.mf_individual` - residual variance
  update. Branch on `residual_variance_method`: `"per_scale_modality"`
  (default) or `"shared_per_modality"` (legacy fidelity mode).
- `update_model_variance.mf_individual` - mixture-weight EM step
  using mixsqp, per the S x M factorization in manuscript derivation
  line 216.
- `get_objective.mf_individual` - ELBO. Manuscript
  eq:elbo_frorm_mean_feild.
- `Eloglik.mf_individual` - expected log-likelihood, separate so we
  can unit-test it in isolation against the manuscript's eq:ERSS.
- `get_var_y.mf_individual`, `get_intercept.mf_individual`,
  `get_fitted.mf_individual`, `get_zscore.mf_individual`,
  `get_variable_names.mf_individual` - standard accessor overrides.
- `initialize_fitted.mf_individual`, `cleanup_model.mf_individual`,
  `trim_null_effects.mf_individual` - lifecycle methods.
- `get_cs.mf_individual` - returns CSs from the standard susieR
  routine (we delegate, no override in v1).

Methods on the `mfsusie` model class (per-effect getter API used by
some tooling):

- `get_alpha_l.mfsusie`, `get_posterior_mean_l.mfsusie`,
  `get_posterior_mean_sum.mfsusie`,
  `get_posterior_moments_l.mfsusie`, `get_prior_variance_l.mfsusie`,
  `set_prior_variance_l.mfsusie`.

The `mfsusie` model class also carries the standard S3 quartet:
`coef.mfsusie`, `predict.mfsusie`, `fitted.mfsusie`,
`summary.mfsusie`, `print.mfsusie`. `coef` and `predict` are pure
transforms of the fit (no X/Y needed, except `predict` which takes
`newx`). `fitted` reads `fit$fitted` (cached at IBSS finalize time).

### D5. Port map (port sources -> mfsusieR)

For every behavior in the two port sources, we record the action and
the paradigm chosen.

| Port source location | Behavior | Decision | Paradigm |
|---|---|---|---|
| `mvf.susie.alpha/R/multfsusie.R:489-589` | IBSS outer loop | Delegate to `susieR::susie_workhorse` | mvsusieR |
| `mvf.susie.alpha/R/multfsusie.R:519` | `cal_Bhat_Shat_multfsusie` | Reimplement as `compute_ser_statistics.mf_individual` | mvsusieR |
| `mvf.susie.alpha/R/EM.R:31-120` | `EM_pi_multsusie` | Reimplement as `update_model_variance.mf_individual` using mixsqp directly. Manuscript eq:S_weigthed_ash_pb. | mvsusieR |
| `mvf.susie.alpha/R/multfsusie.R:541` | `update_multfsusie` (posterior moment update) | Reimplement as `calculate_posterior_moments.mf_individual` | mvsusieR |
| `mvf.susie.alpha/R/operation_on_multfsusie_obj.R:1990-2003` | `update_cal_pip.multfsusie` | Reimplement as `mf_get_pip` (private). Called once at the end of `ibss_finalize.mf_individual`, after all CS operations. | bespoke (fixes PIP-after-filter bug) |
| `mvf.susie.alpha/R/operation_on_multfsusie_obj.R:1452-1454` | `check_cs`/`merge_effect`/`discard_cs` | Reimplement, with the strict ordering contract: filter -> recompute alpha-derived state -> compute PIPs -> finalize | bespoke |
| `mvf.susie.alpha/R/computational_routine.R:395-431` | `estimate_residual_variance` | Reimplement as `update_variance_components.mf_individual`, branched on `residual_variance_method`. Default per-(scale, modality). | mvsusieR with manuscript-aligned default |
| `fsusieR/R/utils.R` (`remap_data`, `colScale`) | Position remapping, column standardization | Port verbatim into `R/utils_wavelet.R`. Behaviour-preserving copy with snake_case naming and a code-quality audit (D13). | bespoke (port-with-audit) |
| `fsusieR/R/computational_functions.R` (`cal_Bhat_Shat`, full-Y case) | Per-position marginal OLS regression | NOT ported locally. Use the new `susieR::compute_marginal_bhat_shat(X, Y, ...)` helper added in the upstream coordination plan (Migration Plan). The same helper replaces mvsusieR's OLS path. mfsusieR's prior init and PR-group-5 per-effect SER both call it. | delegate to upstream susieR (after Migration Plan) |
| `fsusieR/R/wavelet_utils.R` (`gen_wavelet_indx`, DWT helpers) | Scale-index generation, ragged-T_m DWT plumbing | Port into `R/utils_wavelet.R`. `wavethresh::wd` and `wavethresh::wr` are math primitives and stay as library calls. | bespoke (port-with-audit) |
| `fsusieR/R/operation_on_prior.R` (`init_prior.default`) | Per-scale mixture-of-normals prior init | Port into `R/prior_scale_mixture.R`. | bespoke (port-with-audit) |
| `mvf.susie.alpha/R/multfsusie.R:308-316` and `fsusieR/R/susiF_workhorse.R` (DWT pipeline) | Forward DWT pipeline | `mf_dwt(Y_m, pos_m, max_padded_log2, filter_number, family)` in `R/dwt.R`. Calls the ported helpers + `wavethresh::wd`. Cached at data-class construction, not per iteration. | bespoke wrapper |
| `mvf.susie.alpha/R/operation_on_multfsusie_obj.R:2074-2099` and `fsusieR/R/operation_on_susiF_obj.R` (inverse DWT) | Inverse DWT | `mf_invert_dwt(W_m, ...)` in `R/dwt.R`. Called from `coef.mfsusie`, `predict.mfsusie`, and the post-processors. | bespoke wrapper |
| `fsusieR/R/operation_on_susiF_obj.R:1727-1809` (`update_cal_fit_func.susiF`, smash/TI/HMM branches) | Post-processing of effect curves | Port into `R/post_processing.R` as **decoupled** functions: `mf_post_smooth(fit, method = c("smash", "TI", "HMM"))` and `mf_credible_bands(fit, method = "TIWT", level = 0.95)`. Read residuals from `fit$residuals` when present; accept `(X, Y)` opt-in when `save_residuals = FALSE`. See D8c. | bespoke (port-with-audit) |
| `fsusieR/R/operation_on_prior.R` and `mvf.susie.alpha/R/operation_on_multfsusie_prior.R:17-138` | Prior object construction | Reimplement as `mf_prior_scale_mixture`. Per-(scale, modality) prior. Ported `init_prior` callable per modality. | mvsusieR (V_structure analogue) |
| `mvf.susie.alpha/R/multfsusie.R:184` `nullweight = 0.7` | Null component weight | Rename to `null_prior_weight`, default `2`. The original code uses `nullweight = 0.7` and then internally scales by `max(K)` at `EM.R:65`. With typical K = 3, the effective post-scaling value is ~2. Setting the public default directly to 2 makes the behavior interpretable; the internal multiplication-by-K is removed. | bespoke |
| `mvf.susie.alpha/R/multfsusie.R:172-208` and `fsusieR/R/susiF.R:276-309` | Public API | Reimplement as `mfsusie()` (multi-modality) and `fsusie()` (single-modality thin wrapper). Argument names per CLAUDE.md naming rules. | bespoke |
| Y_u (univariate) path throughout `mvf.susie.alpha` | Univariate trait special case | Treat as `T_m = 1` modality. Wavelet machinery short-circuits at construction. No separate code path. | mvsusieR (R = 1 collapse) |
| `mvf.susie.alpha/R/multfsusie.R:284, 287` | `max_scale` parameter | Keep as `max_padded_log2` argument with same meaning. | bespoke (renamed) |
| `mvf.susie.alpha/R/operation_on_multfsusie_obj.R:1716-1735` | Commented-out lfsr code | Drop. v1 does not compute LFSR. A follow-up change can add it. | drop |
| `mvf.susie.alpha/R/multfsusie.R:178-201` (backfit pass) | Re-visit each effect once per outer iter | Drop. susieR's `ibss_fit` already does this every iteration. | delegate to susieR |
| `mvf.susie.alpha/R/multfsusie.R:178-201` (greedy R-search) | Grow L from L_start until pruning | Drop. mvf.susie.alpha and fsusieR implement greedy as a per-iteration interleaved step (`+7`-effect grows + purity-based pruning inside the IBSS loop). mfsusieR uses the susieR-side wrapper instead: outer loop, linear `+L_greedy` steps, `min(lbf) < lbf_min` saturation criterion (slot-invariant, single-round verdict). `mfsusie()` passes `L_greedy` and `lbf_min` straight through to `susieR::susie_workhorse`. The two algorithms are not numerically equivalent; C2 / C3 fidelity tests run with greedy disabled on both sides (D11d). | delegate to susieR |
| `fsusieR/R/EBmvFR.R`, `EBmvFR_workhorse.R`, `operation_on_EBmvFR_obj.R` | EBmvFR algorithm | Out of scope. EBmvFR is a different model (EB multivariate functional regression with no SuSiE structure). Not ported. | drop (out of scope) |

### D6. Prior composition

Default prior is `mf_prior_scale_mixture`. Stored fields:

- `pi`: `list[M]` of `S_m x K` mixture-weight matrices `pi_{k,s,m}`.
- `V_grid`: `list[M]` of length-K grids `sigma_{k,s,m}^2`. Per-scale
  variance is the v1 default; collapsing to per-modality grids is a
  scalar multiplier and stored as a `list[M]` of `S_m`-vectors when
  the user requests `prior_variance_scope = "per_modality"`.
- `null_prior_weight`: scalar. Default per CLAUDE.md naming.
- `update_method`: "mixsqp" (default) or "em".

The cross-modality plug-in seam: at construction time, `mfsusie()`
accepts an optional `cross_modality_prior` argument. If non-NULL,
the SER-stats method calls
`combine_modality_lbfs(cross_modality_prior, modality_lbfs,
model_state)` to adjust per-modality log-BFs before they are
summed into a joint log-BF. The verb is `combine_modality_lbfs`
(an S3 generic) rather than `apply` to avoid collision with
`base::apply`. v1 provides one stub implementation:
`cross_modality_prior_independent()` whose method is a no-op (sums
log-BFs unchanged). A future change adds mash-style adjustments.
This is the susieAnn paradigm.

The per-modality lbf compute that feeds `combine_modality_lbfs` is
itself a per-modality `lapply` over `seq_len(M)` and is thus a
natural target for `future.apply::future_lapply` parallelisation
when M is large. Not in v1 default; flagged as a Phase 7
optimisation candidate.

### D7. Public API signatures

Two public entry points, both manuscript-aligned and CLAUDE.md
naming-rule compliant.

```r
mfsusie(
  X, Y, pos = NULL, L = 10,
  L_greedy = 3, lbf_min = 0.1,
  prior_variance_scope = c("per_scale_modality", "per_modality"),
  prior_variance_grid_multiplier = sqrt(2),
  prior_variance_grid = NULL,
  null_prior_weight = 2,
  cross_modality_prior = NULL,
  residual_variance_method = c("per_scale_modality", "shared_per_modality"),
  estimate_residual_variance = TRUE,
  estimate_prior_variance = TRUE,
  estimate_prior_method = c("optim", "EM", "uniroot"),
  estimate_prior_mixture_weights = TRUE,
  mixture_weight_method = c("mixsqp", "EM"),
  check_null_threshold = 0,
  prior_tol = 1e-9,
  max_padded_log2 = 10,
  max_iter = 100, tol = 1e-3,
  coverage = 0.95, min_abs_corr = 0.5,
  filter_credible_sets = TRUE,
  wavelet_filter_number = 10,
  wavelet_family = "DaubLeAsymm",
  standardize = TRUE, intercept = TRUE,
  precompute_cache = TRUE,
  save_residuals = TRUE,
  n_thread = 1,
  verbose = TRUE, track_fit = FALSE,
  model_init = NULL
)

fsusie(
  Y, X, pos = NULL, L = 10,
  ...
)
```

`fsusie()` accepts a single phenotype as input. `Y` can be a numeric
matrix (`n x T`) or a numeric vector of length `n`; the wrapper
canonicalizes it into the multi-modality shape `Y_canonical = list(Y_or_matrix)`,
sets `pos_canonical = list(pos %||% seq_len(T))`, and forwards to
`mfsusie(X, Y_canonical, pos_canonical, L = L, ...)`. Argument order
matches `fsusieR::susiF` `(Y, X, ...)` for drop-in migration. Every
`mfsusie()` argument that makes sense for the M = 1 case is forwarded
through `...`; arguments that imply multi-modality (e.g.,
`cross_modality_prior`) error if passed to `fsusie()`.

#### Naming convention alignment with susieR and mvsusieR

The signatures follow mvsusieR/refactor-s3 conventions, which
themselves follow susieR's master-branch conventions. Harmonized
choices:

- `tol = 1e-3` matches mvsusieR's default (susieR uses `1e-4` but
  mvsusieR loosened to `1e-3` for the iterative outer loop; we match
  mvsusieR because mfsusieR is a closer kin to mvsusieR's data-class
  pattern).
- `verbose = TRUE` matches mvsusieR's default (susieR uses `FALSE`).
- `mixture_weight_method = c("mixsqp", "EM")` capitalization matches
  mvsusieR exactly.
- `estimate_prior_method = c("optim", "EM", "uniroot")` matches
  mvsusieR (susieR has "optim", "EM", "simple"; the multivariate
  path needed `uniroot` instead of `simple`).
- `check_null_threshold`, `prior_tol`, `n_thread` are inherited from
  mvsusieR with identical names and defaults.
- `precompute_cache` matches mvsusieR.
- `save_residuals` is new to mfsusieR. It is `TRUE` by default; the
  per-modality residuals are needed by the decoupled post-processors
  (D8c). Memory cost is `O(N * sum(T_padded))`, comparable to the
  DWT cache. Opt-out via `save_residuals = FALSE` when the user is
  fitting many models on large data and accepts the verbose
  post-processor signature `mf_post_smooth(fit, X, Y, method)`.
- Argument ordering mirrors mvsusieR: data inputs first, model size,
  prior arguments grouped together, residual-variance arguments,
  estimation method choices, CS options, system knobs (verbose,
  n_thread, track_fit) at the end. `post_processing` is gone (D8c).

#### Renamed arguments

| Source argument | mfsusieR | Note |
|---|---|---|
| `multfsusie`, `susiF` | `mfsusie`, `fsusie` | function names |
| `nullweight` | `null_prior_weight` | default 2 (post-scaling-equivalent of original 0.7 * K = 3) |
| `gridmult` | `prior_variance_grid_multiplier` | |
| `max_scale` | `max_padded_log2` | clearer meaning |
| `cov_lev` | `coverage` | matches susieR |
| `min_purity` | `min_abs_corr` | matches susieR |
| `filter_cs` | `filter_credible_sets` | snake_case, no abbrev |
| `filter.number` | `wavelet_filter_number` | snake_case |
| `family` | `wavelet_family` | scope-prefixed |
| `cal_obj` | dropped | ELBO is computed implicitly when `track_fit = TRUE` or `convergence_method = "elbo"` |
| `max_SNP_EM` | dropped | rolled into `max_iter` and `mixture_weight_method` |
| `greedy`, `backfit` | dropped | susieR's IBSS does both natively |
| `max_step_EM`, `max_step` | dropped | redundant with `max_iter` |
| `init_pi0_w`, `tol_null_prior`, `lbf_min`, `e`, `posthoc`, `cor_small`, `thresh_lowcount` | dropped | over-parameterization; defaults are documented and there is no evidence users tune these |
| `L_start` | renamed to `L_greedy` (default 3) | passed through to `susieR::susie_workhorse`. Saturation threshold `lbf_min` (default 0.1) is also exposed; both go straight through. |
| `multfsusie.obj` | renamed to `model_init` | matches susieR convention |
| `post_processing` | **removed from `mfsusie()` and `fsusie()`** | D8c: post-processing is a separate function on the fit, not a fit-time argument |
| `quantile_trans` | dropped | over-parameterization in `mvf.susie.alpha`; not in `fsusieR::susiF` |

#### Forbidden arguments in mfsusieR

The following names from `mvf.susie.alpha::multfsusie` and
`fsusieR::susiF` MUST NOT appear in any mfsusieR public function
signature, per the naming rules in CLAUDE.md and the contract in
`specs/mf-public-api/spec.md`:

`multfsusie`, `susiF`, `nullweight`, `gridmult`, `max_scale`,
`max_SNP_EM`, `max_step_EM`, `max_step`, `cal_obj`, `cov_lev`,
`min_purity`, `filter_cs`, `filter.number`, `family` (without
`wavelet_` scope), `init_pi0_w`, `tol_null_prior`, `lbf_min`,
`posthoc`, `cor_small`, `thresh_lowcount`, `greedy`, `backfit`,
`multfsusie.obj`, `quantile_trans`, `post_processing`.

A test in `tests/testthat/test_public_api_naming.R` will assert that
none of these names appear in `formalArgs(mfsusie)` or
`formalArgs(fsusie)` so the rule is machine-checkable.

#### Behaviors dropped because susieR provides them natively

| Original location | Behavior | Why dropped |
|---|---|---|
| `multfsusie.R:178-201` (backfit pass) | Re-visit each effect once per outer iter | susieR's `ibss_fit` already does this every iteration |
| `cal_obj` argument controlling whether to compute ELBO | Toggle ELBO calculation | susieR computes ELBO when `convergence_method = "elbo"`; redundant flag |
| `max_step`, `max_step_EM` (separate iteration caps) | Multiple iteration limits | susieR has a single `max_iter`; one limit suffices |
| Bespoke convergence check in `test_stop_cond` | Convergence | susieR's `check_convergence` covers ELBO and PIP modes |
| Custom validate-prior step | Sanity check on V[l] | susieR's `validate_prior` already runs after each per-effect sweep |
| Hand-coded residual computation per modality | Per-effect partial residuals | dispatched through susieR's `compute_residuals` generic; mfsusieR only writes the method body |

The `mfsusie` fit object: class `c("mfsusie", "susie")`. Fields:
`alpha`, `mu`, `mu2`, `KL`, `lbf`, `lbf_variable`, `pi`, `pi_V`,
`V_grid`, `null_prior_weight`, `sigma2`, `elbo`, `niter`,
`converged`, `fitted` (list of M reconstructed curves), `residuals`
(list of M residual matrices when `save_residuals = TRUE`, else
NULL), `intercept` (length M), `pip` (length J), `cs` (list of CSs
with per-effect index), `model_init`, `csd_X`, `wavelet_meta`,
`call`.

### D8. Wavelet pathway organization

DWT is performed at data-class construction, not in the IBSS loop.
The `mf_dwt(Y_m, pos_m, max_padded_log2, filter_number, family)`
helper, defined in `R/dwt.R`, calls in order:

1. `remap_data(Y_m, pos_m, max_padded_log2)` (ported into
   `R/utils_wavelet.R`) to handle ragged or sparsely-sampled
   functional inputs.
2. `colScale(remapped_Y_m)` (ported into `R/utils_wavelet.R`) to
   compute per-column scale factors (cached on `mf_individual` for
   inverse DWT).
3. `wavethresh::wd(.)` per row (math primitive, library call), with
   the ported `gen_wavelet_indx(log2(T_padded[m]))` to produce
   `scale_index[[m]]`.

Inverse DWT is in `mf_invert_dwt(W_m, ...)`, called only from
`coef.mfsusie`, `predict.mfsusie`, and the post-processors (D8c).
For univariate traits (`T_m = 1`), both helpers short-circuit and
pass values through unchanged.

The ported helpers (`remap_data`, `colScale`, `gen_wavelet_indx`,
`init_prior.default` in D6) are behaviour-preserving copies of the
fsusieR originals subjected to the port-quality audit (D13). Style
and naming improvements are allowed; numerical drift is caught by
the C2 fidelity test as a safety net.

### D8b. Caching and `mf_individual` visibility

`mf_individual` is *internal* (not exported). Users do not construct
or hand around `mf_individual` objects. The constructor
`create_mf_individual` is in `R/data_class.R` but only called from
`mfsusie()`. This matches mvsusieR (`mv_individual` is internal,
`create_mvsusie_data` is internal, see
`mvsusieR/R/mvsusie_constructors.R`), susieR
(`individual_data_constructor` is internal), and susieAnn (per-block
susie fits are constructed inside the EM loop, never user-handed).

The user-facing knobs are two booleans:

- `precompute_cache = TRUE` (mvsusieR parity). When TRUE, the DWT is
  computed at the start of `mfsusie()` and cached on the internal
  `mf_individual` for the duration of the fit; the per-SNP column
  scale `csd_X` is cached; per-modality `t(X) %*% X / n` style
  sufficient statistics are cached when reused across iterations.
  Memory grows by `O(N * sum(T_padded[m]))` for the DWT cache and
  `O(J^2)` for any X covariance cache.

- `save_residuals = TRUE`. When TRUE, the per-modality residuals
  `Y_m - X * sum_l (alpha_l . mu_l)` are stored on the returned fit
  object (`fit$residuals`) so post-processors can read them without
  the user passing `(X, Y)` again. Memory grows by `O(N *
  sum(T_padded[m]))`. Opt-out (`save_residuals = FALSE`) is
  documented in `D8c` and changes only the post-processor calling
  convention, not numerics.

#### Cross-fit caching of DWT

Users who fit many models on the same Y (e.g., grid search over L,
prior, or coverage) may want to amortize the DWT cost. None of the
paradigm references expose a user-handed cached object; the right
idiom is to call `mfsusie()` repeatedly with `precompute_cache =
TRUE` and rely on internal caching. A future change can add a
function-level memoization keyed on `digest::digest(Y, pos)` if
benchmarks show this is a real bottleneck. Not in v1; flagged as a
candidate Phase 7 perf optimization.

### D8c. Decoupled post-processing API

Post-processing of effect curves (smashing, TI, HMM smoothing,
credible bands) is NOT a fit-time argument to `mfsusie()` or
`fsusie()`. It is a set of separate functions that take a fit object
and return a fit object with augmented fields.

Public signatures:

```r
mf_post_smooth(fit, method = c("smash", "TI", "HMM"), ..., X = NULL, Y = NULL)
mf_credible_bands(fit, method = "TIWT", level = 0.95, ..., X = NULL, Y = NULL)
```

The `method` argument is a single dispatch entry per
post-processor; the per-method implementations live in
`R/post_processing.R`. The single dispatch entry was chosen for
namespace cleanliness and to make adding a fourth smoother a
non-breaking change.

**Residual contract.** When `fit$residuals` is non-NULL (the default
case, `save_residuals = TRUE`), the post-processor reads residuals
from the fit and ignores `X` and `Y` if they are also passed
(returning a warning if they disagree with the fit's own X/Y
identity). When `fit$residuals` is NULL, the user MUST pass `X` and
`Y` explicitly; the post-processor recomputes the residuals on the
fly and proceeds.

**Why decoupled.** The roxygen on every post-processor explains, in
statgen-writing-style prose, that this design separates fitting (a
variational estimate of effect locations and posterior moments) from
summarizing (a smoothed point estimate of the effect curves). The
fit object is sufficient for most downstream analyses (`coef`,
`predict`, `pip`, `cs`); smoothing and credible bands are
optional refinements. Decoupling lets the user add or change the
smoother without rerunning the IBSS loop, lets new smoothers ship
without revising `mfsusie()`'s signature, and matches the
susieR/mvsusieR convention where `coef`, `predict`, `summary` are
post-fit accessors rather than fit-time arguments.

**Equivalence contract.** For each `method`, the apple-to-apple
contract is

```r
mfsusie(X, Y, pos) |> mf_post_smooth(method = m)
                                ≡   to <= 1e-8
fsusieR::susiF(Y, X, pos, post_processing = m)
```

with the standard caveat that intentional fixes of fsusieR bugs
break the equivalence and the test asserts the deviation explicitly.

### D8d. The `fsusie()` thin wrapper

`fsusie(Y, X, pos, ...)` is a public function in mfsusieR that
exists to give former fsusieR users a drop-in migration path. The
implementation is short:

```r
fsusie <- function(Y, X, pos = NULL, ...) {
  Y_canonical   <- list(if (is.matrix(Y)) Y else matrix(Y, ncol = length(Y)))
  pos_canonical <- list(pos %||% seq_len(ncol(Y_canonical[[1]])))
  forbidden <- c("cross_modality_prior")
  dots <- list(...)
  if (any(names(dots) %in% forbidden)) {
    stop("Arguments not meaningful for fsusie(): ",
         paste(intersect(names(dots), forbidden), collapse = ", "))
  }
  mfsusie(X = X, Y = Y_canonical, pos = pos_canonical, ...)
}
```

The wrapper preserves the `(Y, X, ...)` argument order from
`fsusieR::susiF` for drop-in compatibility; everything else is
forwarded. The numerical contract is the C2 equivalence at tolerance
`<= 1e-8` (D11b). When fsusieR has a bug we fix in mfsusie's M = 1
path, the C2 test asserts the deviation with a comment citing the
OpenSpec change that authorized it.

### D9. PIP / CS ordering

Strict ordering inside `ibss_finalize.mf_individual`:

1. Run final IBSS sweep, finalize alpha and posterior moments.
2. Compute credible sets from current alpha.
3. Apply purity filter (`min_abs_corr`) and any user-requested CS
   filtering.
4. If filtering removed effects, the current `alpha` is updated to
   drop those rows.
5. Compute `pip` from the *post-filter* alpha as
   `pip_j = 1 - prod_l (1 - alpha[l, j])`.
6. Attach `cs` and `pip` to the fit object.

This order is the contract enforced by spec
`mf-credible-sets/spec.md`.

### D10. Modularity: only what is unique to mfsusieR

The package SHALL contain only code that is *unique to* multi-modality
multi-functional SuSiE plus the routines we ported from fsusieR
because mfsusieR replaces fsusieR. Anything that already lives in
`susieR`, `mixsqp`, `ashr`, or `wavethresh` SHALL be delegated, not
reimplemented. The mvsusieR package is the layout model: ~25 S3
method implementations on `mv_individual` and a small number of new
files; everything else is delegated.

The mfsusieR layout SHALL match this shape:

| Concern | Where it lives | Notes |
|---|---|---|
| IBSS outer loop | `susieR::susie_workhorse` | not reimplemented |
| Per-effect single_effect_update orchestration | `susieR` | not reimplemented |
| Convergence check | `susieR::check_convergence` | not reimplemented |
| Prior-variance validation | `susieR::validate_prior` | not reimplemented |
| `validate_prior` calls and the per-effect sweep | `susieR::ibss_fit` | not reimplemented |
| Forward DWT (per-row, per-modality) | `wavethresh::wd` | math primitive, library call |
| Inverse DWT (per-row, per-modality) | `wavethresh::wr` | math primitive, library call |
| Wavelet scale indexing | `R/utils_wavelet.R` | ported from fsusieR (audited) |
| Position remapping for ragged/sparse functional input | `R/utils_wavelet.R` | ported from fsusieR (audited) |
| Per-column scale factors for inverse DWT | `R/utils_wavelet.R` | ported from fsusieR (audited) |
| Per-scale mixture-of-normals prior init | `R/prior_scale_mixture.R` | ported from fsusieR (audited) |
| Mixsqp convex solve | `mixsqp::mixsqp` | called directly |
| Ash-style mixture grid construction | `ashr` helpers | called directly |
| smash / TI / HMM smoothers | `R/post_processing.R` | ported from fsusieR (audited) |
| TIWT credible-band computation | `R/post_processing.R` | ported from fsusieR (audited) |
| Credible-set construction | `susieR::susie_get_cs` | called directly via the standard susie accessor |
| **Per-(scale, modality) residual variance update** | mfsusieR | unique |
| **Per-(scale, modality) prior variance and mixture weight storage** | mfsusieR | unique |
| **Modality-product joint Bayes factor** | mfsusieR | unique |
| **Cross-modality prior plug-in seam** | mfsusieR | unique |
| **`mf_individual` data-class with ragged T_m** | mfsusieR | unique |
| **DWT cache and inverse helpers across modalities** | mfsusieR | unique (thin wrappers around `wavethresh`) |
| **CS-then-PIP ordering invariant** | mfsusieR | unique (fixes the original code's bug) |
| **`mfsusie` model class and S3 method registrations** | mfsusieR | unique |
| **Decoupled post-processing API (`mf_post_smooth`, `mf_credible_bands`)** | mfsusieR | unique (re-architects fsusieR's fit-time post-processing into post-fit transforms) |
| **`fsusie()` thin wrapper** | mfsusieR | unique (single-modality entry point) |

Any future PR that adds code to mfsusieR SHALL be reviewed against
this table. If the proposed code duplicates something in this
column's "not reimplemented" rows, the reviewer SHALL reject in favor
of delegation. Routines marked "ported from fsusieR (audited)" SHALL
have a corresponding line-range citation in `inst/notes/refactor-exceptions.md`
header for any line of the fsusieR original that was intentionally
not carried over (D11e).

### D11. Three apple-to-apple equivalence contracts

mfsusieR is held to three numerical equivalence contracts. Each is
its own test file; each is a hard contract that the fit must satisfy
unless an OpenSpec change explicitly authorizes a deviation that
fixes a port-source bug.

#### D11a. C1: scalar SuSiE (susieR degeneracy)

mfsusie SHALL reduce mathematically to plain `susieR::susie` under
the following degenerate inputs:

- `M = 1` (single modality)
- `T_1 = 1` (univariate trait, DWT short-circuits)
- `prior_variance_grid` of length 1 (single mixture component)
- `null_prior_weight = 0` (no null component)
- `cross_modality_prior = NULL`
- `prior_variance_scope = "per_modality"`
- `residual_variance_method = "shared_per_modality"`
- `L_greedy = NULL`
- no post-processing applied

Under these inputs the model reduces to:

```
D_{1, i, 1, 1} = sum_{r=1..L} x_i' * (gamma_r * f_{r, 1, 1, 1}) + e_{i, 1, 1, 1}
e ~ N(0, sigma^2)
gamma_r ~ Multinom(1, pi)
f_{r, 1, 1, 1} ~ N(0, sigma_0^2)
```

which is the standard susieR model `y = X * b + e` with `b_r ~
gamma_r * N(0, sigma_0^2)`. The DWT collapses to identity (length-1
input, no transform applied), per-modality storage collapses to
scalar storage, and the per-(scale, modality) updates collapse to
the standard SuSiE single-effect regression updates.

C1 is a hard unit test: a fit with the above inputs SHALL match
`susieR::susie(X, y, L = ..., scaled_prior_variance = sigma_0^2 /
var(y), null_weight = 0, ...)` at tolerance `<= 1e-10` on every
numeric output (alpha, mu, mu2, lbf, KL, sigma2, elbo, pip) and
exactly on credible-set membership.

The test serves four purposes: (i) it pins the degeneracy contract;
(ii) it validates the mathematical reduction without relying on a
port source as the reference; (iii) it catches accidental
reintroductions of code that breaks the reduction (e.g., a stray K
multiplier); (iv) it gives mfsusieR a non-port-source ground truth
for the parts of the model that susieR already validates.

#### D11b. C2: single-modality functional SuSiE (fsusieR/susiF degeneracy)

mfsusie SHALL match `fsusieR::susiF` numerically under inputs that
are equivalent to the single-modality functional case:

- `M = 1`, `T_1 > 1` (single functional trait)
- The remaining mfsusie arguments map onto fsusiF defaults; explicit
  alignment of `null_prior_weight`, `prior_variance_grid_multiplier`,
  `wavelet_filter_number`, `wavelet_family`, `max_padded_log2` is
  recorded in the test file header.

C2 is a hard unit test for `L in c(1, 5, 10)` on a fixed `(X, y)`
fixture and seed: a fit from `mfsusie(X, list(y_matrix), pos =
list(seq_len(T_1)), ...)` (or equivalently `fsusie(y_matrix, X, pos,
...)`) SHALL match `fsusieR::susiF(y_matrix, X, pos, ...)` at
tolerance `<= 1e-8` on every numeric output (alpha, mu, mu2, lbf,
lbf_variable, KL, sigma2, elbo, niter, pip) and exactly on
credible-set membership.

When fsusieR has a bug that mfsusieR fixes (e.g., the
PIP-after-CS-filter ordering), the C2 test asserts the deviation
explicitly with a comment citing this OpenSpec change. The fixed
side is the mfsusieR side; the test does not "expect" fsusieR's
buggy output.

The test file is `tests/testthat/test_fsusier_degeneracy.R`. The
header lists the `fsusieR/R/<file>.R#L<lo>-L<hi>` ranges being
compared.

#### D11c. C3: multi-modality functional SuSiE (mvf.susie.alpha fidelity)

For every ported numerical routine (per-effect SER stats, posterior
moments, KL, residual variance update, ELBO, mixture-weight EM, full
end-to-end fit), mfsusie SHALL match `mvf.susie.alpha::multfsusie`
at tolerance `<= 1e-8` on a fixture and a fixed seed when run with
`residual_variance_method = "shared_per_modality"` (legacy mode that
preserves the original code's variance treatment).

C3 is a set of unit tests, one per routine, each in its own
`tests/testthat/test_mvf_alpha_<routine>.R` file. Each file's
header lists the `mvf.susie.alpha/R/<file>.R#L<lo>-L<hi>` ranges
being compared. As with C2, intentional fixes break the equivalence
and the test asserts the deviation with a citation.

### D11d. Test tolerance philosophy: binary apple/orange

Per the code-refactor skill principles
(`~/Documents/obsidian/AI/general/agents/AGENT-code-refactor.md`),
unit-test tolerances are binary, not graduated:

- **Apple-to-apple comparison.** Same algorithm and same code path
  (e.g., contract C2 vs `fsusieR::susiF` at the matching arguments,
  contract C3 vs `mvf.susie.alpha::multfsusie` at
  `residual_variance_method = "shared_per_modality"`, contract C1 at
  the susieR degenerate inputs). Tolerance MUST be `<= 1e-8` (`<=
  1e-10` for the deterministic-intermediates portion of C1 and any
  pure-R-vs-Rcpp test added in Phase 7). Anything looser is a bug to
  investigate, not a tolerance to relax.
- **Apple-to-orange comparison.** Genuinely different algorithms
  (e.g., mfsusie default per-(scale, modality) variance vs the legacy
  shared-per-modality variance, OR mfsusie vs an Rcpp re-implementation
  that diverges in iteration order). Do NOT compare numerically.
  Smoke test only: assert the call returns a fit, the fit has the
  documented shape, no NaNs, ELBO is monotone.
- **Tolerances `1e-2`, `5e-2`, `1e-6` for apple-to-apple are
  forbidden** in this repo. They hide bugs.

Concrete contract table:

| Test | Comparison type | Tolerance |
|---|---|---|
| C1 vs `susieR::susie` (D11a) | apple-to-apple | `<= 1e-10` |
| C2 vs `fsusieR::susiF` (D11b) | apple-to-apple | `<= 1e-8` |
| C3 vs `mvf.susie.alpha::multfsusie` (D11c) | apple-to-apple | `<= 1e-8` |
| Default per-(scale, modality) vs legacy variance | apple-to-orange | smoke test |
| Phase 3 cpp11 kernel vs pure-R oracle (D14) | apple-to-apple | `<= 1e-12` |
| Phase 7 Rcpp vs pure-R reference | apple-to-apple | `<= 1e-10` |
| Post-processor smash/TI/HMM vs `fsusieR::susiF(post_processing = ...)` | apple-to-apple | `<= 1e-8` |
| `mfsusie(L_greedy = K)` vs `fsusieR::susiF(greedy = TRUE)` | apple-to-orange | smoke test (different greedy algorithms) |
| `mfsusie(L_greedy = K)` vs `mvf.susie.alpha::multfsusie(greedy = TRUE)` | apple-to-orange | smoke test (different greedy algorithms) |

C2 and C3 apple-to-apple contracts SHALL be evaluated with the
reference's greedy disabled (`fsusieR::susiF(greedy = FALSE,
backfit = FALSE)` and `mvf.susie.alpha::multfsusie(greedy = FALSE,
backfit = FALSE)`) AND mfsusieR's `L_greedy = NULL`. Both sides
fit at the same fixed `L`, so the algorithmic difference (susieR's
outer-loop wrapper vs the legacy per-iteration interleaved greedy)
is taken out of the comparison. Greedy-L behaviour itself is
covered by smoke tests + property tests (does the wrapper recover
the simulated K?), not by numerical equivalence.

### D11e. Refactor-exceptions doctrine

Every line of the two port sources MUST be accounted for during the
port:

- `mvf.susie.alpha::multfsusie` and the supporting
  `operation_on_multfsusie_*.R`, `EM.R`, `ELBO_mutlfsusie.R`,
  `computational_routine.R`, `utils_wavelet_transform.R`, and the
  parts of `utils.R` and `utils_formatting.R` that the IBSS path
  touches.
- `fsusieR::susiF` and the supporting
  `operation_on_susiF_obj.R`, `wavelet_utils.R`, the parts of
  `operation_on_prior.R` that init the scale-mixture prior, the
  per-method smoothers in `operation_on_susiF_obj.R` (`TI_regression`,
  `smash_regression`, `HMM_regression`), and the parts of `utils.R`
  needed by the susiF path.
- `fsusieR::EBmvFR` and its supporting files (`EBmvFR.R`,
  `EBmvFR_workhorse.R`, `operation_on_EBmvFR_obj.R`) are entirely
  out of scope and SHALL be marked as such in
  `inst/notes/refactor-exceptions.md` once at the file level (no
  per-line walk required for OOS files).

If a line is intentionally NOT ported (because it implements a known
bug, is dead code, is replaced by a delegated upstream helper, or is
part of an out-of-scope file), the omission is documented in
`inst/notes/refactor-exceptions.md`. Each entry has the form:

```
- <port_source>/R/<file>.R:<lo>-<hi>
  Behavior: <one-line description of what the original lines do>
  Decision: omit | replaced-by-<upstream> | deferred-to-<phase> |
            out-of-scope-EBmvFR
  Reason: <one-paragraph justification, citing OpenSpec change or
          manuscript section>
```

The reviewer pass on each Phase 3 PR confirms that any omitted lines
are entered in this file. PRs that omit lines without an entry are
blocked.

### D12. Manuscript cross-references in code; original-code references in tests and dev notes only

Per the code-refactor skill principle: do NOT reference the original
implementations (`mvf.susie.alpha`, `fsusieR`, `susiF`,
`multfsusie`, "old code", historical structures) in main package
code comments or roxygen. Original-code references are
project-internal provenance metadata, not user-facing documentation.
They go in:

1. **`inst/notes/refactor-exceptions.md`** for omitted lines (D11e).
2. **Test file headers**: each
   `tests/testthat/test_<routine>.R` file that compares mfsusieR
   against a port source carries a header block listing the original
   file paths and line ranges. Format:
   ```
   # apple-to-apple comparison against:
   #   mvf.susie.alpha/R/multfsusie.R#L519-L541  (cal_Bhat_Shat path)
   #   mvf.susie.alpha/R/EM.R#L31-L120          (mixture-weight EM)
   ```
   For C2 tests:
   ```
   # apple-to-apple comparison against:
   #   fsusieR/R/susiF.R#L276-L310               (public entry)
   #   fsusieR/R/susiF_workhorse.R#L<lo>-L<hi>   (IBSS body)
   ```
3. **Dev notes** (e.g., `inst/notes/sessions/*.md`,
   `inst/notes/paradigms/mvf-original.md`) where free-form prose
   citation is appropriate.

In *main* package code (anything under `R/`):

- Manuscript citations are required and go in roxygen as
  `@references` blocks with the format
  `Manuscript: methods/<file>.tex eq:<label>` (one citation per
  line). The built-in `@references` tag is used because roxygen2
  rejects unknown tags; the `Manuscript:` prefix makes the
  citations greppable for the reviewer checklist scan.
- Where the manuscript contains a typo or ambiguity, follow the
  reference with a `@manuscript_note` block describing the
  discrepancy and the chosen interpretation.
- `@references_original` tags are FORBIDDEN in `R/` source. So are
  bare strings `mvf.susie.alpha`, `multfsusie`, `fsusieR`, `susiF`,
  "original implementation". The reviewer rejects any PR that
  introduces them. Provenance is available in the test file header
  and `refactor-exceptions.md`.

This separates user-facing documentation (manuscript math) from
implementation provenance (which port-source line was carried over).

### D13. Port-quality audit step

When the author pass copies a routine from a port source
(`mvf.susie.alpha/R/*.R` or `fsusieR/R/*.R`), the author SHALL run a
code-quality audit on the copied routine before advancing to the
reviewer pass. This is a substep of the Claude-native review loop;
the procedure is recorded in
`inst/notes/review-loop-methodology.md` (step 1, "Author pass /
Port-quality audit substep"). Two acceptable mechanisms:

1. Invoke the `simplify` skill on the changed files.
2. Spawn an `Agent` (Explore subagent) with a focused prompt of the
   form: "audit this ported routine for naming, dead code, redundant
   branches, idiomatic R, and consistency with the rest of the
   mfsusieR style; do not change behaviour."

Style improvements logged in the commit message; behaviour-preserving
only. Numerical changes are out of scope and are caught by the
apple-to-apple equivalence test for the relevant contract (C2 for
fsusieR-sourced routines, C3 for mvf.susie.alpha-sourced routines)
as a safety net. Any post-audit numerical drift fails the contract
test and reverts the audit edit.

This step is required of routines copied from port sources. It is
not required of routines that are bespoke to mfsusieR (e.g., the
per-(scale, modality) residual variance update, the cross-modality
prior plug-in seam).

### D14. C++ acceleration via cpp11 in Phase 3 hot paths

The Non-Goals section originally deferred all C++ to Phase 7. That
stance is loosened for narrow utility kernels on the per-effect
SER hot path. The rules:

- **Header**: `LinkingTo: cpp11` only. NOT Rcpp, NOT
  RcppArmadillo, NOT cpp11armadillo. Wholesale matrix-algebra
  acceleration (eigendecompositions, solves, Cholesky) is still
  Phase 7 territory and uses RcppArmadillo when it lands. Phase 3
  cpp11 kernels SHALL be element-wise / dense-array operations on
  `cpp11::doubles_matrix<>` and `cpp11::doubles` only.
- **Pure-R reference oracle**. Every cpp11 kernel SHALL ship with
  a pure-R counterpart in `R/reference_implementations.R` (mvsusieR
  pattern, `mvsusieR/R/reference_implementations.R`). The R version
  is named `<kernel>_R` (e.g., `mixture_log_bf_per_scale_R`); the
  cpp11 version is named `<kernel>` and is the production callable.
- **Apple-to-apple test at `<= 1e-12`**. Each cpp11 kernel has a
  `test_<kernel>.R` test file that asserts agreement with the `_R`
  oracle at tolerance `<= 1e-12` on randomized inputs, fixed seed.
  This catches drift introduced by FMA/SIMD reordering and forces
  the C++ formula to stay literal-equivalent to the R formula.
  This tolerance is tighter than the C2/C3 contract floor of `1e-8`
  to reduce the chance that C++/R drift contaminates a downstream
  contract test.
- **Where cpp11 lives**: `src/<kernel>.cpp`, one kernel per file
  unless they share state. R-side wrappers in `R/<topic>.R` (NOT
  inside the pure-R reference file). Generated R bindings via
  `cpp11::cpp_register()` in `R/cpp11.R` (auto-generated, not
  hand-edited).
- **Out-of-scope for D14**: anything that would call BLAS at the
  C++ level. `X %*% b` and `crossprod(X, R)` stay in R via base
  BLAS. cpp11 is for the per-(SNP, position, mixture-component)
  closed-form arithmetic that is not BLAS-shaped.

The kernels currently authorized under D14 (PR group 5b):

- `mixture_log_bf_per_scale` — per-(SNP, scale) log-Bayes factor
  for the mixture-of-normals prior. Pure-R oracle:
  `mixture_log_bf_per_scale_R`.
- `mixture_posterior_per_scale` — per-(SNP, position) posterior
  mean and second moment under the mixture-of-normals prior.
  Pure-R oracle: `mixture_posterior_per_scale_R`.

Adding a kernel under D14 in a future PR group SHALL be done by
amending this list and the cpp11 task block in tasks.md, NOT
silently. The reviewer pass rejects new `src/*.cpp` files whose
authorization is not recorded here.

**Caches deliberately NOT introduced (measured, not theoretical):**

Per the 5c efficiency review, two redundant computations were
flagged as "high severity" candidates for caching: (a) `Xb_l_m =
X %*% (alpha_l * mu_l[[m]])` is computed in both
`update_fitted_values.mf_individual` and (via `compute_kl`)
`SER_posterior_e_loglik.mf_individual` for the same `(l, iter, m)`;
(b) `colSums(R_m * R_m)` is computed in both `compute_kl` and
`SER_posterior_e_loglik`. We profiled both before deferring:

| Item                | Per-call cost | Per-IBSS-step impact |
|---------------------|--------------:|---------------------:|
| One full IBSS step  |        97 ms |                100% |
| One `X %*% b` call  |       2.4 ms |                ~2.5% |
| One `colSums(R^2)`  |     <0.1 ms |                <0.1% |

(Profile: J=1000, M=1, T=128, n=500.) The `Xb_l_m` cache saves
~2.5% per IBSS step, the `colSums(R^2)` cache <1%. Neither
justifies the schema change (`model$Xb_current[[m]]`) and the
deviation from susieR's pattern (susieR also recomputes both).
Reconsider if profvis in Phase 7 shows different distributions on
larger fixtures.

The motivation is that the per-effect SER step is called L * iter
times per fit, and each call walks J * M * S mixture cells. Pure
R is BLAS-bound only on `X %*% b` and `crossprod(X, R)`; the
mixture inner loop is element-wise and does not benefit from BLAS.
Profiling on the C2 fixture (J = 1000, M = 3, T = 128, K = 10,
L = 10, ~50 IBSS iters) shows the pure-R version at roughly 5-10s
per fit; the cpp11 version is expected at < 1s. Phase 7 may
revisit and unify this with RcppArmadillo when MASH-style
cross-modality priors land.

## Risks / Trade-offs

- *Per-(scale, modality) residual variance default departs from the
  original code.* Trade-off: aligns with manuscript and is the most
  plausible structural fix for the FDR miscalibration, but C3
  fidelity tests then need explicit `residual_variance_method =
  "shared_per_modality"` to compare numerics. Mitigation: the legacy
  mode is a documented argument in v1, not removed; the C3 test
  suite uses it; default-mode tests assert *calibration* properties
  rather than byte-for-byte match.
- *Single data class with ragged T_m may complicate vectorized inner
  loops.* Trade-off: cleaner R-level code, slower on uniformly-sized
  data than a 3D-array path. Mitigation: hot loops are flagged for
  Phase 7 Rcpp work; the data-class shape does not constrain
  implementation choices.
- *Cross-modality plug-in seam without a default implementation.*
  Risk of the seam being under-specified or wrong-shape because
  there is no concrete user. Mitigation:
  `cross_modality_prior_independent()` is shipped as the trivial
  implementation and is exercised in tests; the seam contract is in
  `mf-prior/spec.md`.
- *Porting fsusieR helpers + smoothers expands Phase 3 scope.* Trade-off:
  v1 release date pushed out vs a runtime dependency on fsusieR
  forever. Mitigation: each port lands in its own PR group with the
  D13 audit; the C2 contract serves as a numerical safety net so
  ports stay behaviour-preserving.
- *Decoupled post-processing requires the fit object to carry
  residuals by default.* Trade-off: larger fit objects (an extra
  `O(N * sum(T_padded))` of memory) vs simpler post-processor API.
  Mitigation: opt-out via `save_residuals = FALSE` keeps a smaller
  fit and surfaces the verbose `mf_post_smooth(fit, X, Y, method)`
  call signature. Default is to keep residuals because the typical
  genomics use case is small enough to not notice the cost.
- *EM on mixture weights via mixsqp may not converge for K large or
  pathological data.* Mitigation: fallback to `mixture_weight_method
  = "EM"`, documented in roxygen.
- *Naming break from `multfsusie` / `susiF` to `mfsusie` / `fsusie`
  will surprise existing users.* Mitigation: README documents the
  rename, points existing fsusieR users at `mfsusieR::fsusie()` as
  the migration path; we do not maintain a shim, since this is
  v1.0.0 and there are no users of mfsusieR.

## Migration Plan

This change applies cleanly: there is no prior `inst/openspec/specs/`
content to merge with. After Phase 3 implementation lands and Phase
4 tests pass, the change is archived. Rollback is `git revert` of the
applied commits.

### Upstream susieR `compute_marginal_bhat_shat` (planned)

A new helper `compute_marginal_bhat_shat(X, Y, predictor_weights =
NULL, sigma2 = NULL)` lands in `susieR/R/univariate_regression.R`.
Returns `list(Bhat = J x T, Shat = J x T)` for `Y` either a vector
or matrix, treating each `(X column, Y column)` pair as an OLS
regression with X assumed column-centred (no intercept; `n - 1`
denominator for `Shat`). Rfast-accelerated `colVars` internally,
matrixStats fallback.

Shared by:

- susieR's `compute_ser_statistics.individual` (T = 1 path).
  Refactor cosmetic; behaviour preserved bit-for-bit.
- mvsusieR's OLS path in `compute_ser_statistics.mv_individual`
  (when GLS `data$svs` is NULL). The GLS path keeps its own
  `compute_betahat`. Refactor cosmetic; bit-for-bit preserved.
- mfsusieR's prior init in `R/prior_scale_mixture.R` (full-Y per
  modality) AND PR-group-5 per-effect SER in
  `R/individual_data_methods.R::compute_ser_statistics.mf_individual`
  (residual per modality).

This consolidates the per-position marginal OLS regression that
fsusieR's `cal_Bhat_Shat`, mvsusieR's `compute_betahat` (OLS
branch), and susieR's existing scalar SER all reimplement
separately. Authorised by Gao on 2026-04-25 as a CLAUDE.md hard
rule #1 override (parallel to the L_greedy override). Tracked
outside this OpenSpec change but a hard prerequisite for PR
group 4 of `add-mfsusier-s3-architecture`.

### Upstream susieR `L_greedy` and `lbf_min` (landed)

`susieR::susie_workhorse` carries the greedy outer loop (linear
`+L_greedy` step, `min(lbf) < lbf_min` saturation criterion,
warm-start across rounds via the existing `params$model_init`
mechanism). The change was authorized by Gao on 2026-04-25 as an
override of CLAUDE.md hard rule #1 and landed on `~/GIT/susieR`
master. mfsusieR's `mfsusie()` puts `L_greedy` and `lbf_min` in
the params list; `susie_workhorse` handles the rest. No graceful-
degradation shim, no commit-hash pin (the susieR change is in our
local fork; if it doesn't merge upstream before mfsusieR ships,
mfsusieR's `Imports` line will pin the local version).

The mvf.susie.alpha and fsusieR per-iteration interleaved greedy
algorithm (with `+7`-effect grows and purity-based dummy
detection) is NOT ported. C2 / C3 fidelity tests run with greedy
disabled on both sides per D11d.

## Open Questions

- *Is `cleanup_model.mf_individual` strictly needed in v1?* mvsusieR
  has it; mfsusieR may not need anything beyond the default. Decide
  during Phase 3.
- *Which fsusieR smoother helpers carry over cleanly and which need
  refactoring during the port?* Resolve in PR group 6b
  (post-processing) by walking
  `fsusieR/R/operation_on_susiF_obj.R:1727-1809` line by line under
  the D13 audit.

Resolved questions (recorded for the design history):

- *Default value of `null_prior_weight`* -> 2 (post-scaling-equivalent
  of `mvf.susie.alpha`'s `nullweight = 0.7` * `max(K) = 3`).
- *Should `mfsusie()` accept a precomputed `mf_individual`?* -> No.
  All three paradigm references keep the data class internal; we
  match that. Cross-fit DWT caching is a candidate Phase 7 perf
  optimization.
- *Smoothers* -> Ported into `R/post_processing.R`, decoupled from
  `mfsusie()`'s argument list, accessed via `mf_post_smooth(fit,
  method)`.
- *Greedy R-search* -> Kept, per manuscript online_method step 2.
  susieR does not provide it; mfsusieR wraps `susie_workhorse` with
  the grow-until-prune logic after the upstream `L_greedy`
  generalization lands.
- *fsusieR as a dependency* -> No. mfsusieR replaces fsusieR. The
  routines we need are ported into `mfsusieR/R/` under the D13
  audit; fsusieR is a port source and a numerical reference, not a
  runtime dependency.
- *fsusieR::EBmvFR* -> Out of scope. Different model, no SuSiE
  structure.
- *Argument order of `fsusie()`* -> `(Y, X, pos, ...)` matching
  `fsusieR::susiF` for drop-in compatibility, even though
  `mfsusie()` takes `(X, Y, pos, ...)` in mvsusieR/susieR style.
