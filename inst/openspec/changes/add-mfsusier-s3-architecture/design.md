## Context

mfsusieR ports William Denault's `mvf.susie.alpha::multfsusie` onto the
`susieR` backbone, organized using the S3 conventions of
`mvsusieR/refactor-s3` (paradigm reference #1) and `susieAnn` (paradigm
reference #2). The math anchor is the manuscript at
`/home/gw/GIT/MultifSuSiE_Manuscript`. Phase 1 paradigm notes
(`inst/notes/paradigms/`) catalogued the four packages and identified six
candidate sources for the FDR miscalibration Anjing observed (see
`mvf-original.md` section 8). This design fixes the architecture so Phase 3
implementation has a stable target and Phase 4 reference tests have a clear
contract.

The model: M modalities, each with N samples and T_m positions. Wavelet
transform per modality maps Y_m to D_m of size N x T_m. Per modality, fit a
SuSiE-style sum-of-single-effects on D_m with a per-(scale, modality)
scale-mixture-of-normals prior. Bayes factors combine across modalities by
product (modality-independence assumption, manuscript online_method line
41).

## Goals / Non-Goals

**Goals:**

- Replace `mvf.susie.alpha::multfsusie` with a clean S3 implementation
  delegating to `susieR::susie_workhorse` for the IBSS loop.
- Storage shape: per-(scale, modality) prior variance and mixture weights
  natively supported in `mf_individual` and the `mfsusie` model object.
- Single data class collapses degenerate cases: `M = 1`, `T_m = 1`, and
  T_m varying across modalities all flow through the same code.
- Public API follows CLAUDE.md naming rules. No abbreviations beyond
  `pip`, `cs`, `lbf`, `elbo`, `kl`. Snake_case throughout.
- Pluggable seam for cross-modality covariance prior. Default fits
  modalities independently per the manuscript; a future change provides a
  mash-style implementation.
- PIPs are a deterministic function of the *final* alpha state. CS
  filtering happens before PIP computation, fixing the Phase-1-identified
  PIP-after-filter bug.
- Per-(scale, modality) residual variance is the v1 default. Legacy
  shared-per-modality mode preserved as `residual_variance_method =
  "shared_per_modality"` for Phase 4 byte-for-byte fidelity tests against
  `mvf.susie.alpha`.

**Non-Goals:**

- Sufficient-statistics path. Will be a follow-up change once the
  individual-data path stabilizes.
- Concrete cross-modality covariance prior. v1 ships only the seam; the
  mash-style implementation is a separate change.
- FDR investigation (Phase 5) and the fixes that follow it (Phase 6).
  This change makes those investigations possible without dictating their
  outcome.
- Rcpp ports. Phase 7 only, after the architecture is locked.
- Deprecation of mvf.susie.alpha. mvf.susie.alpha stays a read-only
  reference for the lifetime of mfsusieR, per CLAUDE.md hard rule #1.

## Decisions

### D1. Single data class `mf_individual` covers all degenerate cases

`mfsusieR/refactor-s3` has separate `mv_individual` and `mv_ss`; mvsusieR
collapses R=1 vs R>1 inside one class via shape-aware branches. mfsusieR
follows the same pattern: `mf_individual` is the only data class in v1,
and `M = 1`, `T_m = 1`, ragged T_m across modalities are all handled by
shape branches. `mf_ss` is deferred.

Alternatives considered:

- One class per modality count and one per per-position count: rejected,
  combinatorial explosion.
- Two classes (`mf_individual` for M > 1 functional, `mf_individual_uni`
  for univariate-only): rejected, duplicate method registrations and
  duplicate test surface.

### D2. Storage shape

`mf_individual` carries:

- `Y` as `list[M]` of `n x T_m` matrices (ragged T_m allowed).
- `pos` as `list[M]` of position vectors of length `T_m`.
- `D` as `list[M]` of `n x T_m` wavelet-coefficient matrices (post-DWT,
  post-padding, post-column-scaling). Cached at construction.
- `scale_index` as `list[M]` of integer vectors mapping each column of
  `D[[m]]` to a scale level (output of
  `fsusieR::gen_wavelet_indx(log2(T_m))`).
- `T_padded` as `integer[M]` recording each modality's padded length
  (next power of 2 above original T_m).
- `X`, `n`, `J` standard.

The `mfsusie` model object adds:

- `alpha` as `L x J` matrix (same as susieR).
- `mu`, `mu2` as `list[L]`, each entry a `list[M]`, each inner entry a
  `J x T_padded[m]` matrix. The nested list lets ragged T_m coexist; we
  do not flatten into a 3D array.
- `pi_V` as `list[M]` of `S_m x K` matrices (`pi_{k,s,m}` per
  manuscript). `S_m = log2(T_padded[m])`.
- `V_grid` as `list[M]` of length-K grids of prior variances
  `sigma_{k,s,m}^2`. Stored once per modality; per-scale variance is
  optional and stored as `list[M]` of `S_m x K` matrices when used.
- `sigma2` as either a scalar per modality (legacy mode) or a
  `list[M]` of length-`S_m` vectors (per-scale-per-modality default). A
  single field holds whichever shape is active; the dispatch on
  `residual_variance_method` is at update-component time.
- `lbf`, `lbf_variable`, `KL`, `pip` standard susieR shapes.

### D3. Class hierarchy

```
mf_individual            (data class; inherits from `individual`)
  -> dispatched into susieR::susie_workhorse via S3 methods

mfsusie                  (model class; inherits from `susie`)
  -> coef.mfsusie, predict.mfsusie, summary.mfsusie, print.mfsusie

mf_prior_scale_mixture   (default prior)
  -> per-(scale, modality) scale-mixture-of-normals
mf_prior_cross_modality  (plug-in seam, no default impl in v1)
  -> mash-style covariance across modalities, optional layer
```

### D4. S3 dispatch surface

mfsusieR registers methods on `mf_individual` for every generic that
mvsusieR overrides on `mv_individual`. We do not invent new generics.
The methods registered into `susieR`'s namespace via `zzz.R::.onLoad`:

- `initialize_susie_model.mf_individual` - construct the `mfsusie` model
  with the per-(scale, modality) shapes from D2.
- `compute_residuals.mf_individual` - residuals across (modality, scale,
  position).
- `compute_ser_statistics.mf_individual` - per modality, compute
  per-(scale, position) `Bhat`, `Shat`. Delegates to fsusieR's per-scale
  routines for the wavelet-coefficient regression. Combines per-modality
  log-BFs by sum (manuscript algorithms.tex line 20: product on the BF
  scale, sum on the log-BF scale).
- `calculate_posterior_moments.mf_individual` - mixture posterior over K
  components, per (modality, scale, position). Manuscript eq:post_f_mix,
  eq:post_f2_mix.
- `compute_kl.mf_individual` - KL divergence for the mixture prior.
- `loglik.mf_individual`, `neg_loglik.mf_individual` - marginal
  log-likelihood used by the EM step on mixture weights.
- `update_fitted_values.mf_individual` - store the contribution of effect
  l back into the running fit.
- `update_variance_components.mf_individual` - residual variance update.
  Branch on `residual_variance_method`: `"per_scale_modality"` (default)
  or `"shared_per_modality"` (legacy fidelity mode).
- `update_model_variance.mf_individual` - mixture-weight EM step using
  mixsqp, per the S x M factorization in manuscript derivation line 216.
- `get_objective.mf_individual` - ELBO. Manuscript eq:elbo_frorm_mean_feild.
- `Eloglik.mf_individual` - expected log-likelihood, separate so we can
  unit-test it in isolation against the manuscript's eq:ERSS.
- `get_var_y.mf_individual`, `get_intercept.mf_individual`,
  `get_fitted.mf_individual`, `get_zscore.mf_individual`,
  `get_variable_names.mf_individual` - standard accessor overrides.
- `initialize_fitted.mf_individual`, `cleanup_model.mf_individual`,
  `trim_null_effects.mf_individual` - lifecycle methods.
- `get_cs.mf_individual` - returns CSs from the standard susieR routine
  (we delegate, no override in v1).

Methods on the `mfsusie` model class (for the per-effect getter API used
by some tooling):

- `get_alpha_l.mfsusie`, `get_posterior_mean_l.mfsusie`,
  `get_posterior_mean_sum.mfsusie`,
  `get_posterior_moments_l.mfsusie`, `get_prior_variance_l.mfsusie`,
  `set_prior_variance_l.mfsusie`.

### D5. Delegation map (`mvf.susie.alpha::multfsusie` -> mfsusieR)

For every behavior in the original code, we record the action and the
paradigm chosen.

| Original location | Behavior | Decision | Paradigm |
|---|---|---|---|
| `multfsusie.R:489-589` | IBSS outer loop | Delegate to `susieR::susie_workhorse` | mvsusieR |
| `multfsusie.R:519` | `cal_Bhat_Shat_multfsusie` | Reimplement as `compute_ser_statistics.mf_individual`, delegate per-scale to fsusieR | mvsusieR |
| `EM.R:31-120` | `EM_pi_multsusie` | Reimplement as `update_model_variance.mf_individual` using mixsqp directly. Manuscript eq:S_weigthed_ash_pb. | mvsusieR |
| `multfsusie.R:541` | `update_multfsusie` (posterior moment update) | Reimplement as `calculate_posterior_moments.mf_individual` | mvsusieR |
| `operation_on_multfsusie_obj.R:1990-2003` | `update_cal_pip.multfsusie` | Reimplement as `mf_get_pip` (private). Called once at the end of `ibss_finalize.mf_individual`, after all CS operations. | bespoke (fixes Phase 1 PIP-after-filter bug) |
| `operation_on_multfsusie_obj.R:1452-1454` | `check_cs`/`merge_effect`/`discard_cs` | Reimplement, with the strict ordering contract: filter -> recompute alpha-derived state -> compute PIPs -> finalize | bespoke |
| `computational_routine.R:395-431` | `estimate_residual_variance` | Reimplement as `update_variance_components.mf_individual`, branched on `residual_variance_method`. Default per-(scale, modality). | mvsusieR with manuscript-aligned default |
| `multfsusie.R:308-316` | DWT pipeline | Port to `mf_dwt` private helper (R/dwt.R). Uses `fsusieR::remap_data`, `fsusieR::colScale`, `wavethresh::wd`, `fsusieR::gen_wavelet_indx`. Cached at data-class construction, not per iteration. | bespoke wrapper, fsusieR delegation |
| `operation_on_multfsusie_obj.R:2074-2099` | Inverse DWT | Port to `mf_invert_dwt` private helper. Same pieces, called from `predict.mfsusie` and from CS-band computation. | bespoke wrapper |
| `operation_on_multfsusie_obj.R:1456-1487, 2401-2477` | smash / TI / HMM post-processing | Keep. Delegate to fsusieR's smash/TI/HMM helpers where they exist; reimplement only where fsusieR does not provide an equivalent. Argument `post_processing = c("none", "smash", "TI", "HMM")` on `mfsusie()`. TIWT band computation per manuscript online_method step 6 is part of the same path. | mvsusieR (delegation), bespoke wrapper |
| `operation_on_multfsusie_prior.R:17-138` | Prior object construction | Reimplement as `mf_prior_scale_mixture`. Per-(scale, modality) prior. fsusieR's `init_prior.default` is delegated to per modality. | mvsusieR (V_structure analogue) |
| `multfsusie.R:184` `nullweight = 0.7` | Null component weight | Rename to `null_prior_weight`, default `2`. The original code uses `nullweight = 0.7` and then internally scales by `max(K)` at `EM.R:65` (`nullweight_scaled <- nullweight * max(K)`). With typical K = 3 mixture components, the effective post-scaling value is ~2. Setting the public default directly to 2 makes the behavior interpretable without requiring users to reason about K. The internal multiplication-by-K is removed. | bespoke |
| `multfsusie.R:172-208` | Public API | Reimplement as `mfsusie()`. Argument names per CLAUDE.md naming rules. | bespoke |
| Y_u (univariate) path throughout | Univariate trait special case | Treat as `T_m = 1` modality. Wavelet machinery short-circuits at construction (no DWT). No separate code path. | mvsusieR (R = 1 collapse) |
| `multfsusie.R:284, 287` | `max_scale` parameter and remapping | Keep as `max_padded_log2` argument with same meaning. | bespoke (renamed) |
| `R/operation_on_multfsusie_obj.R:1716-1735` | Commented-out lfsr code | Drop. v1 does not compute LFSR. A follow-up change can add it. | drop |
| `multfsusie.R:178-201` | Backfit (revisit each effect each iteration) | Drop. susieR's `ibss_fit` already revisits every effect each outer iteration via `single_effect_update` (see `susieR/R/iterative_bayesian_stepwise_selection.R:177-189`). The explicit re-fit pass in `mvf.susie.alpha` is redundant with susieR's existing IBSS. | delegate to susieR |
| `multfsusie.R:178-201` | Greedy R-search (grow L from L_start until pruning) | Generalize susieR. Per Gao's direction (overriding CLAUDE.md hard rule #1 for this case), the greedy-grow-L behavior is added to susieR core under the new argument `L_greedy = NULL`. When `L_greedy` is `NULL` (default), susieR's behavior is unchanged. When `L_greedy = K` for some integer K, susieR fits with L = K, checks for pruned effects, grows L by 1, and re-fits, repeating until one effect prunes or `L` (the upper bound) is reached. mfsusieR adds `L_greedy` as a passthrough argument with default `3` (per manuscript online_method step 2). The susieR generalization is its own work item, tracked outside this OpenSpec change but a hard prerequisite for Phase 3. | delegate to susieR (after generalization) |

### D6. Prior composition

Default prior is `mf_prior_scale_mixture`. Stored fields:

- `pi`: `list[M]` of `S_m x K` mixture-weight matrices `pi_{k,s,m}`.
- `V_grid`: `list[M]` of length-K grids `sigma_{k,s,m}^2`. Per-scale
  variance is the v1 default; collapsing to per-modality grids is a
  scalar multiplier and stored as a `list[M]` of `S_m`-vectors when the
  user requests `prior_variance_scope = "per_modality"`.
- `null_prior_weight`: scalar. Default per CLAUDE.md naming.
- `update_method`: "mixsqp" (default) or "em".

The cross-modality plug-in seam: at construction time, `mfsusie()` accepts
an optional `cross_modality_prior` argument. If non-NULL, the SER-stats
method calls `cross_modality_prior$apply(modality_lbfs, model_state)` to
adjust per-modality log-BFs before they are summed into a joint log-BF.
v1 provides one stub implementation: `cross_modality_prior_independent()`
which is a no-op (sums log-BFs unchanged). A future change adds
mash-style adjustments. This is the susieAnn paradigm.

### D7. Public API signature

Manuscript-aligned, CLAUDE.md naming-rule compliant. Defaults are
suggestions; final values land with this proposal.

```r
mfsusie(
  X, Y, pos = NULL, L = 10,
  L_greedy = 3,
  prior_variance_scope = c("per_scale_modality", "per_modality"),
  prior_variance_grid_multiplier = sqrt(2),
  prior_variance_grid = NULL,
  null_prior_weight = 2,
  cross_modality_prior = NULL,
  residual_variance_method = c("per_scale_modality", "shared_per_modality"),
  estimate_residual_variance = TRUE,
  estimate_prior_mixture_weights = TRUE,
  mixture_weight_method = c("mixsqp", "em"),
  post_processing = c("none", "smash", "TI", "HMM"),
  max_padded_log2 = 10,
  max_iter = 100, tol = 1e-4,
  coverage = 0.95, min_abs_corr = 0.5,
  filter_credible_sets = TRUE,
  wavelet_filter_number = 10,
  wavelet_family = "DaubLeAsymm",
  standardize = TRUE, intercept = TRUE,
  precompute_cache = TRUE,
  verbose = FALSE, track_fit = FALSE,
  model_init = NULL
)
```

#### Renamed arguments

| `mvf.susie.alpha` | `mfsusieR` | Note |
|---|---|---|
| `multfsusie` | `mfsusie` | function name |
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
| `greedy`, `backfit` | dropped | susieR's IBSS does both natively (see D5) |
| `max_step_EM`, `max_step` | dropped | redundant with `max_iter` |
| `init_pi0_w`, `tol_null_prior`, `lbf_min`, `e`, `posthoc`, `cor_small`, `thresh_lowcount` | dropped | over-parameterization; defaults are documented and there is no evidence users tune these |
| `L_start` | renamed to `L_greedy` (default 3) | passed through to susieR after the susieR generalization lands; `NULL` means non-greedy (use L as fixed upper bound) |
| `multfsusie.obj` | renamed to `model_init` | matches susieR convention |

#### Forbidden arguments in mfsusieR

The following names from `mvf.susie.alpha::multfsusie` MUST NOT appear
in any mfsusieR public function signature or argument list, per the
naming rules in CLAUDE.md and the contract in
`specs/mf-public-api/spec.md`:

`multfsusie`, `nullweight`, `gridmult`, `max_scale`, `max_SNP_EM`,
`max_step_EM`, `max_step`, `cal_obj`, `cov_lev`, `min_purity`,
`filter_cs`, `filter.number`, `family` (without `wavelet_` scope),
`init_pi0_w`, `tol_null_prior`, `lbf_min`, `posthoc`, `cor_small`,
`thresh_lowcount`, `greedy`, `backfit`, `multfsusie.obj`.

A test in `tests/testthat/test-public-api-naming.R` will assert that
none of these names appear in `formalArgs(mfsusie)` so the rule is
machine-checkable.

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
`V_grid`, `null_prior_weight`, `sigma2`, `elbo`, `niter`, `converged`,
`fitted` (list of M reconstructed curves), `intercept` (length M), `pip`
(length J), `cs` (list of CSs with per-effect index), `model_init`,
`call`, `wavelet_meta` (for inverse-DWT helpers).

### D8. Wavelet pathway organization

DWT is performed at data-class construction, not in the IBSS loop. The
`mf_dwt(Y_m, pos_m, max_padded_log2, filter_number, family)` helper
wraps `fsusieR::remap_data`, `fsusieR::colScale`, `wavethresh::wd`, and
records `T_padded[m]`, `scale_index[[m]]`, and the column scales for
later inversion. Inverse DWT is in `mf_invert_dwt`, called only from
`predict.mfsusie` and `coef.mfsusie` and from the optional TIWT band
helper. For univariate traits (`T_m = 1`), the DWT short-circuits and
both helpers pass through.

### D8b. Caching and `mf_individual` visibility

`mf_individual` is *internal* (not exported). Users do not construct or
hand around `mf_individual` objects. The constructor
`create_mf_individual` is in `R/data_class.R` but only called from
`mfsusie()`. This matches mvsusieR (`mv_individual` is internal,
`create_mvsusie_data` is internal, see `mvsusieR/R/mvsusie_constructors.R`),
susieR (`individual_data_constructor` is internal), and susieAnn
(per-block susie fits are constructed inside the EM loop, never
user-handed). The pattern across all three paradigm references is
identical, so we follow it.

The user-facing knob is the boolean `precompute_cache = TRUE`. When
`TRUE`:

- DWT is computed at the start of `mfsusie()` and cached on the
  internal `mf_individual` for the duration of the fit.
- Per-modality `t(X) %*% X / n` style sufficient statistics are cached
  if they are reused across iterations (analogous to mvsusieR's
  eigendecomposition cache).
- Memory grows by `O(N * sum(T_padded[m]))` for the DWT cache and
  `O(J^2)` for any X covariance cache.

When `FALSE`, DWT and any sufficient statistics are recomputed on
demand (slower, lower memory). The default is `TRUE` for parity with
mvsusieR.

#### Cross-fit caching of DWT

Users who fit many models on the same Y (e.g., grid search over L,
prior, or coverage) may want to amortize the DWT cost. None of the
paradigm references expose a user-handed cached object; instead, the
right idiom is to call `mfsusie()` repeatedly with `precompute_cache =
TRUE` and rely on internal caching. A future change can add a
function-level memoization keyed on `digest::digest(Y, pos)` if
benchmarks show this is a real bottleneck. Not in v1; flagged as a
candidate Phase 7 perf optimization.

### D9. PIP / CS ordering

Strict ordering inside `ibss_finalize.mf_individual`:

1. Run final IBSS sweep, finalize alpha and posterior moments.
2. Compute credible sets from current alpha.
3. Apply purity filter (`min_abs_corr`) and any user-requested CS
   filtering.
4. If filtering removed effects, the current `alpha` is updated to drop
   those rows.
5. Compute `pip` from the *post-filter* alpha as
   `pip_j = 1 - prod_l (1 - alpha[l, j])`.
6. Attach `cs` and `pip` to the fit object.

This order is the contract enforced by spec
`mf-credible-sets/spec.md`.

### D10. Modularity: only what is unique to mfsusieR

The package SHALL contain only code that is *unique to* multi-modality
multi-functional SuSiE. Anything that already lives in `susieR`,
`fsusieR`, `mixsqp`, `ashr`, or `wavethresh` SHALL be delegated, not
reimplemented. mvsusieR is the model: it adds ~25 S3 method
implementations on `mv_individual` and a small number of new files
(`mvsusie_constructors.R`, `mixture_prior.R`, `predict.mvsusie.R`,
`zzz.R`); everything else is delegated to susieR or mashr.

The mfsusieR layout SHALL match the same shape:

| Concern | Where it lives | Notes |
|---|---|---|
| IBSS outer loop | `susieR::susie_workhorse` | not reimplemented |
| Per-effect single_effect_update orchestration | `susieR` | not reimplemented |
| Convergence check | `susieR::check_convergence` | not reimplemented |
| Prior-variance validation | `susieR::validate_prior` | not reimplemented |
| `validate_prior` calls and the per-effect sweep | `susieR::ibss_fit` | not reimplemented |
| Forward DWT (per-modality) | `wavethresh::wd` via `fsusieR` helpers | wrapped, not reimplemented |
| Inverse DWT | `wavethresh::wr` via `fsusieR` helpers | wrapped, not reimplemented |
| Wavelet scale indexing | `fsusieR::gen_wavelet_indx` | called directly |
| Per-scale mixture-of-normals prior | `fsusieR::init_prior.default` | called directly |
| Mixsqp convex solve | `mixsqp::mixsqp` | called directly |
| Ash-style mixture grid construction | `ashr` helpers | called directly |
| smash / TI / HMM smoothers | `fsusieR` if exposed there | wrapped, not reimplemented |
| Credible-set construction | `susieR::susie_get_cs` | called directly via the standard susie accessor |
| **Per-(scale, modality) residual variance update** | mfsusieR | unique |
| **Per-(scale, modality) prior variance and mixture weight storage** | mfsusieR | unique |
| **Modality-product joint Bayes factor** | mfsusieR | unique |
| **Cross-modality prior plug-in seam** | mfsusieR | unique |
| **`mf_individual` data-class with ragged T_m** | mfsusieR | unique |
| **DWT cache and inverse helpers across modalities** | mfsusieR | unique (thin wrappers) |
| **CS-then-PIP ordering invariant** | mfsusieR | unique (fixes the original code's bug) |
| **`mfsusie` model class and S3 method registrations** | mfsusieR | unique |

Any future PR that adds code to mfsusieR SHALL be reviewed against this
table. If the proposed code duplicates something in this column's "not
reimplemented" rows, the reviewer SHALL reject in favor of delegation.

### D11. Mathematical degenerate case: mfsusie reduces to susie

mfsusie SHALL reduce mathematically to plain `susieR::susie` under the
following degenerate inputs:

- `M = 1` (single modality)
- `T_1 = 1` (univariate trait, DWT short-circuits)
- `prior_variance_grid` of length 1 (single mixture component, no
  scale-mixture)
- `null_prior_weight = 0` (no null component in the mixture)
- `cross_modality_prior = NULL` (no cross-modality prior)
- `prior_variance_scope = "per_modality"` (no per-scale variance)
- `residual_variance_method = "shared_per_modality"` (single global
  residual variance)
- `post_processing = "none"` (no smoothing)
- `L_greedy = NULL` (fixed L)

Under these inputs, the model reduces to:

```
D_{1, i, 1, 1} = sum_{r=1..L} x_i' * (gamma_r * f_{r, 1, 1, 1}) + e_{i, 1, 1, 1}
e ~ N(0, sigma^2)
gamma_r ~ Multinom(1, pi)
f_{r, 1, 1, 1} ~ N(0, sigma_0^2)
```

which is the standard susieR model `y = X * b + e` with
`b_r ~ gamma_r * N(0, sigma_0^2)`. The DWT collapses to identity
(length-1 input, no transform applied), per-modality storage collapses
to scalar storage, and the per-(scale, modality) updates collapse to
the standard SuSiE single-effect regression updates.

This degenerate case is a hard unit test: a fit with the above inputs
SHALL match `susieR::susie(X, y, L = ..., scaled_prior_variance =
sigma_0^2 / var(y), null_weight = 0, ...)` at tolerance 1e-10 on every
numeric output (alpha, mu, mu2, lbf, KL, sigma2, elbo, pip, cs).

The test serves four purposes: (i) it pins the degeneracy contract;
(ii) it validates the mathematical reduction without relying on
`mvf.susie.alpha` as the reference; (iii) it catches accidental
reintroductions of code that breaks the reduction (e.g., a stray K
multiplier); (iv) it gives mfsusieR a non-`mvf.susie.alpha` ground truth
for the parts of the model that susieR already validates.

The next-most-degenerate case worth pinning as a separate test:

- `M = 1`, `T_1 > 1` (single functional modality) reducing to fSuSiE.
  This requires fsusieR as a reference and is recorded as a follow-up
  test; it is not the susieR degeneracy.

### D12. Manuscript cross-references in code

Every non-trivial formula in the codebase carries a roxygen tag
`@manuscript_ref methods/<file>.tex eq:<label>`. Where the manuscript
contains a typo or ambiguity, the roxygen tag is followed by a
`@manuscript_note` block describing the discrepancy and the chosen
interpretation. For ported routines, an additional
`@references_original mvf.susie.alpha/R/<file>.R#L<lo>-L<hi>` is added
per CLAUDE.md.

## Risks / Trade-offs

- *Per-(scale, modality) residual variance default departs from the
  original code.* Trade-off: aligns with manuscript and is the most
  plausible structural fix for the FDR miscalibration, but Phase 4
  fidelity tests then need explicit `residual_variance_method =
  "shared_per_modality"` to compare numerics. Mitigation: the legacy
  mode is a documented argument in v1, not removed; the test suite uses
  it for all reference comparisons; default-mode tests assert
  *calibration* properties rather than byte-for-byte match.
- *Single data class with ragged T_m may complicate vectorized inner
  loops.* Trade-off: cleaner R-level code, slower on uniformly-sized
  data than a 3D-array path. Mitigation: hot loops are flagged for
  Phase 7 Rcpp work and can keep a per-modality vectorization
  internally; the data-class shape does not constrain implementation
  choices.
- *Cross-modality plug-in seam without a default implementation.* Risk
  of the seam being under-specified or wrong-shape because there is no
  concrete user. Mitigation: `cross_modality_prior_independent()` is
  shipped as the trivial implementation and is exercised in tests; the
  seam contract is in `mf-prior/spec.md` and any future implementation
  must satisfy it.
- *Smoothers (smash / TI / HMM) are kept and delegated to fsusieR
  where possible.* Trade-off: mfsusieR depends on fsusieR for the
  smoother backends, which means the v1 release tracks fsusieR
  versions for those code paths. Mitigation: a small fallback
  reimplementation lands for any smoother fsusieR does not currently
  expose; CI pins a fsusieR commit hash.
- *EM on mixture weights via mixsqp may not converge for K large or
  pathological data.* Mitigation: fallback to `update_method = "em"`,
  documented in roxygen.
- *Naming break from `multfsusie` to `mfsusie` will surprise existing
  users.* Mitigation: README documents the rename and points to
  `mvf.susie.alpha` for the legacy package; we do not maintain a
  shim, since this is v1.0.0 and there are no users of mfsusieR.

## Migration Plan

This change applies cleanly: there is no prior `inst/openspec/specs/`
content to merge with. After Phase 3 implementation lands and Phase 4
tests pass, the change is archived. Rollback is `git revert` of the
applied commits.

### External coordination: susieR `L_greedy` generalization

This proposal assumes a generalization of `susieR::susie_workhorse` to
accept an `L_greedy` argument (default `NULL`, no behavior change). The
generalization is OUT OF SCOPE for this OpenSpec change but is a hard
prerequisite for the public mfsusieR API to ship as designed.

Coordination plan:

1. Branch `feature/L_greedy` on `~/GIT/susieR` (created 2026-04-25,
   off `master`). Develop the `L_greedy` addition there: add the
   argument to `susie_workhorse`, `susie()`, and any other public
   entry that takes `L`. Behavior unchanged when `L_greedy = NULL`.
   Open a PR to `stephenslab/susieR` once the implementation is
   reviewed in-house.
2. Pin the susieR commit hash that contains the generalization in
   mfsusieR's `DESCRIPTION` `Imports: susieR (>= <hash-or-version>)`
   line. Until the upstream PR merges, mfsusieR pins to the
   feature-branch hash; once merged, the pin moves to the upstream
   master release that contains it.
3. If the upstream PR review takes longer than Phase 3, mfsusieR
   ships v1.0.0 with the pre-generalization graceful-degradation
   warning. Once upstream lands, a follow-up OpenSpec change
   (`enable-mfsusier-greedy-l-passthrough`) flips the warning to a
   real passthrough.

Until the susieR change lands, mfsusieR's `L_greedy` argument is
accepted but ignored with a one-time `lifecycle::deprecate_warn()`-style
warning. The `mf-public-api/spec.md` records this transitional state.

This is an explicit override of CLAUDE.md hard rule #1 ("susieR is
read-only reference") authorized by Gao on 2026-04-25. The override is
narrow: it covers only the addition of `L_greedy` and any internal
plumbing required to support it. All other susieR code remains
read-only.

## Open Questions

- *Is `cleanup_model.mf_individual` strictly needed in v1?* mvsusieR
  has it; mfsusieR may not need anything beyond the default. Decide
  during Phase 3.
- *Which fsusieR smoother helpers exist and which need
  reimplementation?* Resolve in PR group 7 (post-processing) by
  reading fsusieR's exported function list.

Resolved questions (recorded for the design history):

- *Default value of `null_prior_weight`* -> 2 (post-scaling-equivalent
  of `mvf.susie.alpha`'s `nullweight = 0.7` * `max(K) = 3`).
- *Should `mfsusie()` accept a precomputed `mf_individual`?* -> No.
  All three paradigm references keep the data class internal; we match
  that. Cross-fit DWT caching is a candidate Phase 7 perf optimization.
- *Smoothers* -> Kept, delegated to fsusieR where possible.
- *Greedy R-search* -> Kept, per manuscript online_method step 2.
  susieR does not provide it; mfsusieR wraps `susie_workhorse` with the
  grow-until-prune logic.
