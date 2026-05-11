# Design: unit-level APA step-SuSiE

## Locked decisions

- **Naming**: use `unit` / `analysis unit` everywhere in public API,
  documentation, and tests. Do not name the model after a particular
  unit-construction strategy.
- **Core transform**: APA uses a step-function / 0th-order
  trend-filtering basis. It is Haar-like, but not the existing DWT
  module. The step basis is a first-class transform module, analogous
  to DWT in software structure.
- **Default orientation**: use upstream step functions
  `S[j, t] = 1(x_j <= c_t)` so positive coefficients are directly
  interpretable as upstream/PAS contribution.
- **Intercept handling**: match susieR's intercept strategy. Do not
  add an intercept predictor to IBSS. Subtract offsets, then use
  weighted centering of `Y` and each step predictor so the fit is
  equivalent to a unit-specific baseline model.
- **Default effect prior**: simple one-sided Gaussian slab with a
  scalar `tau2`. If `tau2 = NULL`, initialize it once from marginal
  positive step estimates before IBSS. Keep it fixed during the first
  implementation's IBSS loop. No mixture, no grid, no scale mixture
  in this change.
- **Annotation prior**: annotations supply nonnegative prior weights
  over candidate PAS positions. They do not impose effect direction
  and do not hard-filter candidates. Strong annotation is represented
  through strong prior odds, controlled by an explicit strength
  parameter and a positive prior floor.
- **Sharing**: PAS locations are shared across analysis units through
  the common `alpha[l, t]`; effect magnitudes remain unit-specific as
  `beta_il`.
- **Unit weights**: default unit weights are all one. User-supplied
  nonzero unit weights are normalized to have mean one before use;
  zero weights are allowed and contribute no shared location evidence.
- **Candidate pre-scan**: candidate PAS can be supplied directly or
  produced by an optional high-recall `apa_prescan()` function. The
  pre-scan is an input-variable screening step, not final APA
  inference and not a strong coverage-derived prior.
- **Initialization**: the default warm start is a matrix-free
  marginal positive initializer. Optional `L0Learn` initialization is
  off by default, only runs when explicitly requested, must not change
  the posterior target, and must have a safe fallback.
- **Efficiency**: the default fit must use cumulative-sum step
  operators. Dense explicit `J x T` matrices are allowed only for
  tests, debugging, and bounded-size optional initialization.
- **PIP and credible sets**: use susieR posterior semantics for PIP
  and CS construction. Skip dense correlation-purity filtering by
  default; optional purity filtering can use a bounded-size
  step-basis correlation path.

## Review of existing code paths

### susieR trend filtering

The relevant susieR pattern is `susie_trendfilter(order = 0)`.
It constructs a special design object with matrix metadata and relies
on cumulative-sum helpers for `X %*% b`, `t(X) %*% y`, and column
norms. It also uses an optional MAD residual-variance initialization
to reduce poor local optima. The APA implementation should reuse the
operator idea, not the exact object: APA candidates may be a subset of
positions, candidate coordinates carry biological priors and outputs,
and the model has multiple analysis units with shared locations.

#### Direct susieR code map

Implementation agents should inspect and mirror the following
source-level patterns.

1. **Trend-filtering wrapper**:
   `susieR/R/susie_trendfilter.R::susie_trendfilter()` creates an
   empty sparse matrix, tags it with `attr(X, "matrix.type") =
   "tfmatrix"` and `attr(X, "order")`, and lets downstream generic
   multiplication detect the special design. APA should use the same
   separation of concerns: `apa_step_basis()` carries positions,
   candidate indices, orientation, and names; the IBSS methods call
   APA-specific operators rather than materializing `S`.
2. **Fast step operators**:
   `susieR/R/susie_trendfilter_utils.R::compute_tf_Xb()` and
   `compute_tf_Xty()` are cumulative-sum implementations of
   multiplication by the 0th-order trend-filtering design and its
   transpose. APA should implement the analogous
   `apa_step_Xb()` and `apa_step_Xtr()` using prefix or suffix sums,
   with an explicit orientation convention. Candidate subsets are
   handled by indexing the cumulative sums at candidate positions.
3. **Column statistics without dense matrices**:
   `compute_tf_d()`, `compute_tf_cm()`, `compute_tf_csd()`, and
   `susieR/R/susie_utils.R::compute_colstats()` show how centering,
   scaling, and column norms are cached as attributes. APA should
   expose the same invariants through `apa_step_colnorm()` and, when
   needed, weighted centered denominators for each analysis unit.
   The main fit should never need `crossprod(S)`.
4. **Generic multiplication interface**:
   `susieR/R/sparse_multiplication.R::compute_Xb()`,
   `compute_Xty()`, and `compute_MXt()` are the templates for
   dispatching dense, sparse, and special matrices through one
   interface. APA can either call its own `apa_step_*` helpers inside
   `compute_residuals.mf_apa_individual()` or use the same attribute
   dispatch style, but all tests must compare the fast path to
   `apa_step_explicit()`.
5. **Residual-variance initialization**:
   `susie_trendfilter()` optionally runs a first SuSiE fit with
   `residual_variance = 0.5 * median(abs(diff(y))/0.6745)^2`, then
   passes that fit as `model_init`. For APA this is optional and
   should be treated only as a warm start. If used, apply it per
   analysis unit to the working residual `Y_i - offset_i`, skip it when the user
   supplies `model_init`, and fall back to variance-based
   initialization when the MAD estimate is zero or non-finite.
6. **Coefficient-based initialization**:
   `susieR/R/susie_get_functions.R::susie_init_coef()` builds a
   one-hot `alpha`, `mu`, and `mu2` object from selected coefficient
   indices. APA's marginal and optional L0 initializers should mirror
   that shape after adapting it to mfsusieR's list-of-list posterior moments:
   `alpha[l, t_l] = 1`, `mu[[l]][[i]][t_l, 1] = beta_init[i,l]`,
   and `mu2[[l]][[i]][t_l, 1] = beta_init[i,l]^2`.
7. **Model-init plumbing**:
   `susieR/R/iterative_bayesian_stepwise_selection.R::ibss_initialize()`,
   `susieR/R/susie_utils.R::validate_init()`,
   `adjust_L()`, and `prune_single_effects()` define the rules for
   validating warm starts, expanding or preserving `L`, and resetting
   iteration-specific fields. mfsusieR already has the analogous
   `expand_model_init_to_L()` path in `R/ibss_methods.R`; APA should
   reuse that path and add only shape checks specific to candidate
   locations and per-unit positive effects.
8. **Optional L0Learn initialization vignette**:
   `susieR/vignettes/l0_initialization.Rmd` uses
   `L0Learn::L0Learn.cvfit(X, y, penalty = "L0")`, chooses
   `which.min(fit$cvMeans[[1]])`, extracts nonzero coefficients
   after removing the intercept, and converts them with
   `susie_init_coef()`. APA should follow this sequence only when
   `init_method = "l0learn"` is explicitly requested, only on the
   explicit sparse step matrix under the configured size threshold,
   and then combine selected candidates across analysis units.
9. **Posterior normalization**:
   `susieR/R/single_effect_regression.R::apply_ser_lbf()` and
   `susieR/R/susie_utils.R::lbf_stabilization()` /
   `compute_posterior_weights()` show the stable update
   `alpha_l = softmax(lbf + log(pi))`. APA's positive-prior
   log-BFs should feed into the same stable normalization logic after
   cross-unit logBF combination.

### mfsusieR

mfsusieR already has the right software pattern. The current
`mfsusie()` wrapper builds a transformed data object, initializes
prior objects, and then calls susieR's IBSS workhorse through S3
methods. The DWT module is one transform; APA should be another. The
core implementation should therefore be a new APA data subclass plus
small method overrides for places where the inherited methods assume
dense `X`, wavelet coefficient arrays, or symmetric normal effects.
The existing `prior_weights` and `combine_outcome_lbfs()` mechanisms
are direct extension points for candidate priors and cross-unit
sharing.

## Minimal implementation plan

The shortest stable implementation is the core APA phenotype model:
implement matrix-free step
operators, pre-scan, fixed or empirically initialized scalar `tau2`,
half-normal SER updates, weighted cross-unit logBF combination, and
phenotype extraction. This scope should produce correct APA phenotypes
using fixed candidate priors.

## Statistical model

All functions operate on one gene at a time. The gene subscript is
dropped in the code and documentation for lighter notation.

Let `i = 1, ..., n` index analysis units, `j = 1, ..., J` index
ordered observed positions/bins, and `t = 1, ..., T` index candidate
PAS / breakpoint positions. Let `x_j` be the coordinate of position
`j`, and `c_t` be the coordinate of candidate PAS `t`.

For each analysis unit `i`, after offsets and bias corrections are
applied outside or inside the wrapper, the conceptual Gaussian model is

```text
Y_i(j) = offset_i(j) + baseline_i
         + sum_{l=1}^L beta_il S[j, Z_l] + e_i(j),
e_i(j) ~ approximately Gaussian,
beta_il >= 0.
```

`Y_i(j)` is the working coverage outcome, not the raw count. It is the
coverage after library-size normalization, bias correction, and any
variance-stabilizing or residualizing transformation.
The APA fit starts after this preprocessing. The wrapper may carry
offsets, precision weights, and user-supplied corrected outcomes, but
this change does not implement a unique library-size, GC-bias, or
5'/3' bias correction pipeline.

The default upstream step basis is

```text
S[j, t] = 1(x_j <= c_t).
```

The baseline is handled exactly in the spirit of susieR's
individual-level constructor: it is not included as a variable in
IBSS. Define `ytilde_i(j) = Y_i(j) - offset_i(j)`, using zero offset
when none is supplied. Let `w_ij` be row precision weights, defaulting
to one. The APA data constructor computes

```text
ybar_i   = sum_j w_ij ytilde_i(j) / sum_j w_ij
Sbar_it  = sum_j w_ij S[j,t]      / sum_j w_ij
Yc_i(j)  = ytilde_i(j) - ybar_i
Sc_it(j) = S[j,t]      - Sbar_it.
```

The IBSS loop operates on the centered model

```text
Yc_i(j) = sum_{l=1}^L beta_il Sc_i,Z_l(j) + e_i(j).
```

This is algebraically equivalent to fitting a unit-specific baseline
outside the sparse effects. The stored baseline is metadata for
reconstruction and diagnostics, not a single-effect candidate.

The single-effect location variable is shared:

```text
Z_l in {1, ..., T}
Pr(Z_l = t) = pi_t.
```

The annotation prior supplies

```text
omega_t = h_eta(F_t),
pi_t    = omega_t / sum_s omega_s.
```

The unit-specific effect prior is

```text
beta_il | Z_l = t ~ N_+(0, tau2)
```

where `N_+` is a zero-mean normal slab truncated to positive values.
The first version uses one positive scalar `tau2` for each fitted
gene. If the user does not supply `tau2`, initialize it once from
marginal positive step estimates before the IBSS loop. Per-unit or
per-group scale mixtures are intentionally out of scope.

Strong annotation should be encoded as prior odds, not as forced usage.
For example, a score vector `q_t` can be mapped to

```text
pi_t = (1 - T * prior_floor) * softmax(prior_strength * q)_t
       + prior_floor.
```

Large `prior_strength` creates strong prior odds. The positive
`prior_floor` keeps all candidates available to the data. A user may
choose an extremely small floor for an annotation-dominated analysis,
but the default must remain soft rather than hard. This prior is over
locations conditional on an effect being active; it must not force
every analysis unit to use an annotated PAS.

## Relation to Haar and trend filtering

The APA basis is exactly the 0th-order trend-filtering changepoint
basis, up to orientation and candidate restriction. In susieR,
`susie_trendfilter(order = 0)` creates a special trend-filtering
matrix and uses cumulative-sum algebra so operations such as `X' y`
are linear in the number of positions.

Haar wavelets also represent piecewise-constant signals, and the
default mfsusieR wavelet filter can be Haar. However, a candidate PAS
at an arbitrary transcript coordinate is not generally one coefficient
in a fixed dyadic Haar DWT. For APA phenotyping, the coefficient index
must be the candidate PAS itself because priors, PIPs, credible sets,
and usage outputs are candidate-level. Therefore the implementation
should reuse the mfsusieR transform-module pattern, but should not
force the APA model into the existing DWT module.

## Data layout

The APA wrapper maps the problem into mfsusieR's multi-response shape:

- observations `N = J`: ordered positions/bins;
- predictors `p = T`: candidate PAS step functions;
- outcomes `M = n`: analysis units;
- each outcome matrix `Y[[i]]` is `J x 1`.

The APA data class must be the S3 subclass:

```r
c("mf_apa_individual", "mf_individual", "individual")
```

It should carry enough information for existing generics without
requiring dense `data$X` multiplication. It must include:

- `basis`: step-basis metadata and candidate indices;
- `D`: list of unit response matrices, each `J x 1`;
- `Y_center`: unit-specific weighted mean of `Y - offset`;
- `step_center`: per-unit, per-candidate weighted step means
  `Sbar_it`;
- `centered_d`: per-unit, per-candidate centered weighted
  denominators `d_it`;
- `residuals`: list of unit residual matrices;
- `xtx_diag_list`: per-unit candidate column norms;
- `na_idx`: per-unit complete-position indices;
- `unit_weights`: normalized numeric length `n`;
- `row_weights`: row precision weights for each unit, if supplied;
- `sigma2_init`: residual-variance initialization per unit;
- `apa_diagnostics`: centering, low-information, scaling, and
  fallback diagnostics;
- `prior`: location prior object;
- any standard fields expected by inherited mfsusieR methods.

## Step-basis module

File: `R/apa_step_basis.R`.

Required constructor:

```r
apa_step_basis <- function(pos,
                           candidates,
                           orientation = c("upstream", "downstream"),
                           candidate_names = NULL)
```

Required operations:

```r
apa_step_Xb <- function(basis, b)
apa_step_Xtr <- function(basis, r, weights = NULL)
apa_step_colnorm <- function(basis, weights = NULL)
apa_step_center_cache <- function(basis, weights = NULL)
apa_step_explicit <- function(basis, sparse = TRUE)
```

`apa_step_explicit()` exists for tests, debugging, and bounded-size
optional initialization only. The fitting path must call the operator
functions.

For upstream orientation, implementation sketch:

```text
apa_step_Xtr:
  z_j = r_j or weights_j * r_j
  prefix = cumsum(z_j over ordered positions)
  return prefix[candidate_position_index]

apa_step_Xb:
  work = numeric(J)
  add coefficient b_t at the position index corresponding to c_t
  return reverse cumulative sum over ordered positions

apa_step_colnorm:
  prefix_w = cumsum(weights), or seq_len(J) when weights is NULL
  return prefix_w[candidate_position_index]
```

Downstream orientation uses suffix sums or the equivalent sign /
complement transformation. Duplicate candidates are disallowed unless
the constructor explicitly collapses them with a documented rule.

`apa_step_center_cache()` precomputes all basis quantities that do not
change across IBSS iterations:

```text
Q_i       = sum_j w_ij
Sbar_it  = sum_j w_ij S[j,t] / Q_i
d_it     = sum_j w_ij (S[j,t] - Sbar_it)^2.
```

For upstream steps, these are obtained from prefix sums at candidate
position indices:

```text
P_it     = sum_{j:x_j <= c_t} w_ij
Sbar_it  = P_it / Q_i
d_it     = P_it - P_it^2 / Q_i.
```

Candidates with `d_it` below `min_denominator` or with too little
effective weight on either side should be flagged as low-information
for that unit. The fit may keep the candidate for other units, but
its SER statistics for the low-information unit must be made
uninformative.

Candidate mapping rule: if `candidates` are integer indices, use them
directly after bounds checking. If they are coordinates, map each
candidate to the last ordered position satisfying `x_j <= c_t` for
upstream orientation, and to the first valid downstream-side position
for downstream orientation. Candidates with no observed position on
one side are flagged or removed according to the same
low-information rule. The constructor must store both original
coordinates and mapped position indices.

## Candidate PAS pre-scan module

File: `R/apa_prescan.R`.

Required public helper:

```r
apa_prescan <- function(Y,
                        pos,
                        offset = NULL,
                        weights = NULL,
                        grid = NULL,
                        unit_weights = NULL,
                        tau2_scan = NULL,
                        min_side_weight = NULL,
                        min_denominator = NULL,
                        region_half_width = 100,
                        merge_distance = 25,
                        max_peaks = NULL,
                        annotations = NULL,
                        tail_sites = NULL,
                        motif_sites = NULL,
                        keep_local_background = TRUE,
                        ...)
```

The helper returns a candidate set for the main APA wrapper. It is
optional: users may still supply `candidates` directly. The pre-scan
must not fit the final APA model, return final usage phenotypes, or
replace the Bayesian posterior.

Notation:

- `d_r`, `r = 1, ..., R`, is the scan grid.
- By default `grid = NULL` means dense scan: `R = J` and `d_r = x_r`.
- `G[j, r] = 1(x_j <= d_r)` uses the same upstream orientation as
  the default APA step basis.
- `m_ir` is the unit-level scan score, a positive-prior log Bayes
  factor.
- `shat2_ir` is the working sampling variance of the marginal
  `bhat_ir^scan`; it is not a score.

Offset and intercept rules:

1. `Y` is already the working coverage outcome.
2. If `offset` is supplied, subtract it before scanning. It must be a
   fixed position-specific nuisance correction with the same shape as
   `Y`, or a vector recycled across units by a documented rule.
3. If `offset` is `NULL`, use zero offset.
4. Always fit a weighted intercept in each marginal scan. The
   intercept is implemented by weighted centering and absorbs average
   coverage level, so no user-supplied average offset is needed.
5. Do not estimate a flexible smooth offset inside `apa_prescan()`;
   smooth GC, mappability, or broad 5'/3' bias correction belongs in
   preprocessing or in an externally supplied offset.

For one unit and one scan position, the regression is

```text
ytilde_ij = gamma_ir + beta_ir^scan G[j, r] + e_ij,
beta_ir^scan >= 0.
```

with `ytilde = Y - offset`. Weighted centering gives

```text
N_ir =
  sum_j w_ij (G[j,r] - Gbar_ir) (ytilde_ij - ybar_i)

D_ir =
  sum_j w_ij (G[j,r] - Gbar_ir)^2

bhat_ir^scan = N_ir / D_ir

shat2_ir = sigma2_i / D_ir.
```

`N_ir` is the centered weighted covariance between the step predictor
and working coverage. `D_ir` is the centered weighted information in
the step predictor. The denominator must be guarded by
`min_denominator`; scan positions with too little effective weight
upstream or downstream must be flagged and skipped using
`min_side_weight`.

The unit score is

```text
m_ir = log BF^+(bhat_ir^scan, shat2_ir, tau2_scan).
```

If `tau2_scan` is `NULL`, estimate it once from the scan marginal
statistics using the same positive-excess rule as the main
`tau2` initializer, with its own scan diagnostic label. This keeps
the scan scale on the working-outcome scale without introducing a
separate tuning grid.

The pooled score curve is

```text
M_r = sum_i unit_weights_norm[i] * m_ir.
```

Candidate regions are seeded from local maxima in the pooled curve.
By default, retain finite local maxima with `M_r > 0`. If `max_peaks`
is finite, keep the highest-scoring retained maxima up to that cap.
If a user wants a stratified sensitivity scan, the same function can
be called separately within those strata, but this is not part of the
core API. Local maxima are expanded to

```text
[d_peak - region_half_width, d_peak + region_half_width].
```

The default `region_half_width = 100` is a high-recall default for
region construction, not a single-nucleotide precision claim. It is
wider than common poly(A) signal and cleavage microheterogeneity
windows. Nearby positions within `merge_distance`, default exactly
25 nt, are collapsed using a documented representative rule to avoid
treating cleavage microheterogeneity as separate APA events. Users may
disable this merging by setting `merge_distance = NULL` or
`merge_distance <= 0`.

The returned object must contain:

- `candidates`: numeric candidate coordinates for `apa_susie()`;
- `candidate_metadata`: source labels such as `scan_region`,
  `scan_background`, `annotation`, `tail`, or `motif`;
- `regions`: scan-seeded candidate regions with peak position, source
  curve, peak score, expanded interval, and retained candidates;
- `scores`: optional scan score summaries or a lazily returned score
  object;
- `diagnostics`: skipped positions, denominator failures,
  low-side-weight failures, caps applied, and final candidate counts.

Within each expanded scan region, keep fine-scale positions and
flanking low-score positions when `keep_local_background = TRUE`.
These are not guaranteed null effects, but they provide nearby
competing variables and reduce conditioning on only pre-selected
maxima. Coverage-derived scan scores should normally be used for
candidate inclusion. Annotation, motif, tail-read support, and other
external evidence are the preferred inputs for strong final
`prior_weights`, to reduce double use of the same coverage data.

Implementation should use cumulative sums, not `R` separate dense
regressions. For dense `R = J`, compute for each unit
`P_ir = sum_{j:x_j <= d_r} w_ij`,
`H_ir = sum_{j:x_j <= d_r} w_ij * ytilde_ij`,
`T_i = sum_j w_ij * ytilde_ij`, and `Q_i = sum_j w_ij`; then
`N_ir = H_ir - P_ir * T_i / Q_i` and
`D_ir = P_ir - P_ir^2 / Q_i`. These are the same WLS quantities
above, written in prefix-sum form. This gives `O(n * J)` cost per gene
for a dense scan, or `O(n * (J + R))` for a custom grid. The scan can
be parallelized over genes and analysis units.

## Annotation prior module

File: `R/apa_prior_annotation.R`.

Required helpers:

```r
apa_prior_uniform <- function(candidates)
apa_prior_from_scores <- function(scores,
                                  prior_strength = 1,
                                  prior_floor = NULL)
apa_prior_from_annotations <- function(features,
                                       coefficients = NULL,
                                       intercept = 0,
                                       prior_strength = 1,
                                       prior_floor = NULL)
```

The output is always a numeric vector of length `T`, strictly
positive and summing to one. All helpers must validate finite values
and fall back to uniform only when requested explicitly.

Annotation features are evidence about candidate locations. They must
not encode effect signs. Hard filtering, if explicitly requested by a
user before calling the wrapper, is outside this module. The first
implementation should be deliberately simple: support direct
`prior_weights`, direct numeric `scores`, or a small fixed linear
score from supplied feature columns and supplied coefficients. It must
not train a genome-wide prediction model or classifier.

Implementation notes:

- `prior_strength` multiplies the annotation score before softmax.
- `prior_floor` defaults to a small positive value such as
  `.Machine$double.eps` or a documented larger value for numerical
  stability. The implementation must reject `prior_floor < 0` and
  `T * prior_floor >= 1`.
- Direct `prior_weights` supplied to the wrapper are normalized after
  adding or enforcing the same positive floor.
- For diagnostics, the prior object should store prior odds summaries
  such as `max(pi) / min(pi)` so users can see how strong annotation
  was in a fit.

## Positive prior module

File: `R/apa_prior_positive.R`.

Required helpers:

```r
apa_halfnormal_lbf <- function(bhat, shat2, tau2)
apa_halfnormal_moments <- function(bhat, shat2, tau2)
apa_halfnormal_ser <- function(bhat, shat2, tau2)
```

Numerical requirements:

- guard `shat2` and `tau2` at positive machine-safe lower bounds;
- compute `log(2 * Phi(x))` stably for large negative `x`;
- posterior second moments must be non-negative and finite;
- return zero-BF / zero-moment behavior for no-information candidates.

The downstream orientation path may sign-flip marginal estimates into
the upstream-positive convention, use the same positive-prior
functions, and convert back only if needed. The default APA phenotype
uses upstream positive coefficients.

## Initialization module

File: `R/apa_init.R`.

Required helpers:

```r
apa_marginal_init <- function(Y,
                              basis,
                              L,
                              unit_weights = NULL,
                              weights = NULL,
                              offset = NULL,
                              tau2 = NULL,
                              max_nonzero = L,
                              ...)

apa_estimate_tau2 <- function(bhat,
                              shat2,
                              unit_weights = NULL,
                              y_scale = NULL,
                              tau2_floor = NULL,
                              quantile = 0.5)

apa_l0learn_init <- function(Y,
                             basis,
                             L,
                             unit_weights = NULL,
                             positive = TRUE,
                             lambda_choice = c("cv_min", "cv_1se"),
                             max_nonzero = L,
                             max_explicit_entries = 5e7,
                             ...)

apa_init_from_step_coef <- function(candidate_index,
                                    beta_init,
                                    p,
                                    L)
```

Default marginal behavior:

1. Form `ytilde_i = Y_i - offset_i`, with zero offset when none is
   supplied.
2. Center `ytilde_i` by its weighted mean and use the centered step
   cache. The default initializer must follow the same no-explicit
   intercept convention as the main fit.
3. For each analysis unit, compute marginal weighted step estimates
   and sampling variances for all candidates using cumulative-sum
   operators and cached centered denominators, not a dense step
   matrix:

   ```text
   bhat_it = sum_j w_ij Sc_it(j) Yc_i(j) / d_it
   shat2_it = sigma2_i / d_it.
   ```

   Low-information `d_it` entries are treated as uninformative.
4. Compute a positive-prior marginal score for each unit-candidate
   pair, such as `logBF_it^+ = apa_halfnormal_lbf(bhat_it,
   shat2_it, tau2)`.
5. Pool candidate support across units as
   `score_t = sum_i unit_weights_norm[i] * logBF_it^+`, with unit
   weights defaulting to one and normalized before use.
6. Keep at most `min(L, max_nonzero, T)` candidates by pooled score.
7. Create a SuSiE initialization with one effect per retained
   candidate: `alpha[l, t_l] = 1`, with unit-specific starting
   effects from positive marginal estimates.
8. For a retained candidate with non-positive marginal evidence in
   unit `i`, initialize that unit's coefficient at zero.
9. Return an object accepted by the APA wrapper's `model_init`
   plumbing, plus metadata recording selected candidates, scores, and
   the initializer used.

Default `tau2` initialization:

1. If the user supplies a positive scalar `tau2`, use it unchanged
   and record `tau2_source = "user"`.
2. Otherwise, use the same marginal `bhat_it` and `shat2_it` computed
   for `apa_marginal_init()`.
3. For each unit-candidate pair, compute the positive excess signal
   estimate
   `e_it = max(bhat_it, 0)^2 - shat2_it`.
4. Keep finite `e_it > 0` values and weight them by
   normalized unit weights, defaulting to one.
5. Set `tau2` to a robust weighted quantile of the retained positive
   excess values, defaulting to the weighted median. Enforce
   `tau2 >= tau2_floor`.
6. If no finite positive excess values are available, set
   `tau2 = tau2_floor` and record
   `tau2_source = "floor_no_positive_excess"`.
7. The default `y_scale` is the median across analysis units of the
   finite weighted variance of centered `Y - offset`. The default
   `tau2_floor` should be
   `max(.Machine$double.eps, 1e-8 * y_scale)`, with a fallback
   `y_scale = 1` if no finite positive scale is available.
8. Store `tau2`, `tau2_source`, `tau2_floor`, the quantile used, and
   the number of positive excess estimates in fit diagnostics.

This rule is deliberately a one-time initialization, not empirical
Bayes inside IBSS. It gives the positive-prior SER update a stable
scale without running per-effect or in-loop scale optimization.

Residual variance initialization:

1. If the user supplies `residual_variance`, validate it as positive
   and use it as the initial `sigma2_i`.
2. Otherwise, for each unit use a weighted analogue of susieR
   trend-filtering's MAD-difference initializer on the centered
   working coverage:

   ```text
   sigma2_i^MAD = 0.5 * (median(|diff(Yc_i)|) / 0.6745)^2.
   ```

3. If the MAD value is zero or non-finite, fall back to the weighted
   variance of `Yc_i`.
4. If both are unavailable, use a small positive floor and record a
   diagnostic.
5. Residual variance may be updated by
   `update_model_variance.mf_apa_individual()` when
   `estimate_residual_variance = TRUE`; this is separate from the
   fixed `tau2` slab scale.

Optional L0Learn behavior:

1. If `L0Learn` is not installed, return `NULL` plus a diagnostic.
2. If `J * T > max_explicit_entries`, return `NULL` plus a diagnostic.
3. Build `S <- apa_step_explicit(basis, sparse = TRUE)`.
4. For each analysis unit `i`, form `ytilde_i = Y_i - offset_i`,
   with zero offset when none is supplied. If row weights are present,
   use the weighted least-squares equivalent
   `S_w = sqrt(w_i) * S` and `y_w = sqrt(w_i) * ytilde_i`;
   otherwise use `S` and `ytilde_i`.
5. Fit `L0Learn::L0Learn.cvfit(S_w, y_w, penalty = "L0", ...)`.
6. Choose the lambda index using `lambda_choice`. For `cv_min`, use
   `which.min(fit$cvMeans[[1]])`, matching the susieR vignette.
7. Extract coefficients with `coef(fit$fit, lambda = selected_lambda)`,
   remove the intercept, and map nonzero entries to candidate indices.
8. Under `positive = TRUE`, discard coefficients `<= 0`. Do not
   truncate negative coefficients to positive values in the default
   implementation.
9. Score each candidate across units as
   `score_t = sum_i unit_weights_norm[i] * max(beta_hat_it, 0)^2`,
   with unit weights defaulting to one and normalized before use.
10. Keep at most `min(L, max_nonzero, T)` candidates by score.
11. Create a SuSiE initialization with one effect per retained
    candidate: `alpha[l, t_l] = 1`. Store unit-specific posterior
    mean starts from the L0 coefficients when available.
12. For a retained candidate missing in unit `i`, initialize
    `beta_it` by the positive marginal WLS estimate at that candidate,
    or zero if the estimate is not positive.
13. Return an object accepted by the APA wrapper's `model_init`
    plumbing, plus metadata recording selected candidates, scores,
    lambda rule, and fallback diagnostics.

This module is allowed to materialize a sparse explicit step matrix
only for optional `init_method = "l0learn"` and only under the safety
threshold. The default marginal initializer and main fitting loop must
remain matrix-free.

## Cross-unit module

File: `R/apa_cross_unit.R`.

Constructor:

```r
cross_unit_prior_weighted <- function(unit_weights = NULL)
```

S3 method:

```r
combine_outcome_lbfs.mf_prior_cross_unit_weighted <-
  function(prior, outcome_lbfs, model_state)
```

The method first validates and normalizes weights. If `unit_weights`
is `NULL`, use all ones. If supplied, weights must be finite,
non-negative, length `n`, and not all zero. The positive weights are
rescaled so their mean is one:

```text
unit_weights_norm[i] = unit_weights[i] /
                       mean(unit_weights[k] : unit_weights[k] > 0).
```

Zero weights remain zero. This convention preserves the all-one
default and prevents arbitrary global rescaling of evidence.

The method returns

```r
Reduce("+", Map(function(w, lbf) w * lbf, unit_weights_norm, outcome_lbfs)).
```

## APA data subclass and wrapper

File: `R/apa_mfsusie.R`.

Candidate public wrapper:

```r
apa_susie <- function(Y,
                      pos,
                      candidates,
                      L = 10,
                      prior_weights = NULL,
                      annotations = NULL,
                      unit_weights = NULL,
                      tau2 = NULL,
                      orientation = c("upstream", "downstream"),
                      prior_strength = 1,
                      prior_floor = NULL,
                      init_method = c("marginal", "none", "l0learn"),
                      l0_control = list(),
                      weights = NULL,
                      offset = NULL,
                      residual_variance = NULL,
                      estimate_residual_variance = TRUE,
                      ...)
```

Input contract:

- `Y` is either a numeric matrix `J x n` or a list of `n` numeric
  vectors/matrices each length `J`;
- `Y` is the working coverage outcome, not raw counts; preprocessing
  and bias correction are performed upstream or supplied through
  `offset` and `weights`;
- `pos` is length `J` and ordered in transcript direction;
- `candidates` are positions or integer indices mappable to `pos`;
  they may be user-supplied directly or taken from
  `apa_prescan(...)$candidates`;
- `prior_weights` overrides annotation-derived weights when supplied;
- `annotations` is optional and must have one row per candidate;
- `unit_weights` is optional and length `n`; if supplied, it is
  normalized to positive mean one before being used in shared
  log-BF combination or initialization scores;
- `tau2` is scalar and positive; if `NULL`, initialize from marginal
  step estimates;
- `residual_variance` is optional and initializes per-unit residual
  variance; if `NULL`, initialize from centered working coverage using
  MAD-difference with weighted-variance fallback;
- `estimate_prior_variance` is not a first-version public APA
  argument; the slab scale is controlled by `tau2`;
- `prior_strength` and `prior_floor` control how strongly annotation
  affects candidate prior odds while keeping candidates available;
- `init_method = "marginal"` is the default matrix-free warm start;
- `init_method = "l0learn"` requests the optional L0 warm start.

The wrapper builds the APA data object, sets `prior_weights`, installs
the weighted cross-unit combiner, chooses the positive prior path,
creates `model_init` via `apa_marginal_init()` unless
`init_method = "none"` or a user-supplied `model_init` is provided,
optionally tries `apa_l0learn_init()` only when requested, and calls
`susieR::susie_workhorse()`.

The wrapper must not pass an explicit intercept variable into IBSS.
It stores weighted response centers and step centers on the APA data
object and implements all residual, fitted-value, and SER-statistic
methods using centered quantities.

## S3 methods to implement for `mf_apa_individual`

The subclass inherits all unchanged behavior from `mf_individual`.
Override a method only where dense multiplication, wavelet semantics,
or the Gaussian/mixed prior is wrong.

Required overrides:

- `compute_residuals.mf_apa_individual`
- `update_fitted_values.mf_apa_individual`
- `loglik.mf_apa_individual`
- `calculate_posterior_moments.mf_apa_individual`
- `SER_posterior_e_loglik.mf_apa_individual`
- `pre_loglik_prior_hook.mf_apa_individual` as a no-op for fixed
  `tau2`
- `post_loglik_prior_hook.mf_apa_individual` as a no-op for fixed
  `tau2`
- `update_variance_components.mf_apa_individual`
- `Eloglik.mf_apa_individual`
- `update_model_variance.mf_apa_individual` for residual-variance
  updates and derived-cache refresh, mirroring
  `update_model_variance.mf_individual` but using APA residual
  structures.

Also implement `compute_ser_statistics.mf_apa_individual()` explicitly,
even if it is thin. It should use the cached `apa_step_Xtr()` residual
inner products, `apa_step_colnorm()` denominators, unit-specific
weights, and the current scalar residual variance. This avoids an
accidental fall-through to dense or wavelet-specific assumptions.

Register the new S3 methods using the same mechanism as the
corresponding `mf_individual` methods. Methods for susieR generics
that mfsusieR already registers in `.onLoad` should be added there.
Methods such as `pre_loglik_prior_hook` and `post_loglik_prior_hook`,
which are currently registered through roxygen-generated NAMESPACE
`S3method(...)` entries for `mf_individual`, should be exported and
registered the same way for `mf_apa_individual`.

## Positive-prior SER update

Within `loglik.mf_apa_individual`, per-unit logBFs are computed from
`ser_stats$betahat[[i]]`, `ser_stats$shat2[[i]]`, and `tau2` via
`apa_halfnormal_lbf()`. These per-unit logBF vectors are combined by
`combine_outcome_lbfs()`. The resulting joint logBF is combined with
`model$pi` exactly as in `loglik.mf_individual` to update `alpha`.

For effect `l`, compute the centered partial residual

```text
Rc_il(j) = Yc_i(j) - sum_{l' != l} sum_t alpha_l't mu_il't Sc_it(j).
```

The SER sufficient statistics are

```text
x_it      = sum_j w_ij Sc_it(j) Rc_il(j)
d_it      = sum_j w_ij Sc_it(j)^2
betahat_it = x_it / d_it
shat2_it   = sigma2_i / d_it.
```

`d_it`, `Sbar_it`, low-information flags, and candidate mapping are
cached in the data object. At each IBSS update only residual-dependent
inner products `x_it` are recomputed using prefix or suffix sums.
Low-information unit-candidate pairs receive uninformative logBFs and
zero posterior moments for that unit.

Within `calculate_posterior_moments.mf_apa_individual`, per-unit
posterior means and second moments are computed by
`apa_halfnormal_moments()`.

The first implementation keeps `tau2` fixed during IBSS. The wrapper
must set the APA prior path so `pre_loglik_prior_hook.mf_apa_individual`
and `post_loglik_prior_hook.mf_apa_individual` return the current
`tau2` unchanged. Residual variance `sigma2` may still be updated by
`update_model_variance.mf_apa_individual()` when
`estimate_residual_variance = TRUE`; that hook is not the slab
variance update.

## Phenotype module

File: `R/apa_phenotype.R`.

Required helper:

```r
apa_phenotype <- function(fit,
                          type = c("usage", "balance", "expected_length", "all"),
                          cutpoints = NULL,
                          eps = 1e-8,
                          return_uncertainty = TRUE,
                          export_for_mixed_model = FALSE)
```

Expected returned object:

- `candidate_pip`: candidate-level PIP;
- `credible_sets`: PAS / breakpoint credible sets;
- `effect_mean`: unit-by-candidate posterior mean contribution;
- `usage`: unit-by-candidate normalized usage;
- `balance`: unit-by-cutpoint log proximal/distal balance, if requested;
- `expected_length`: unit-level expected 3'UTR coordinate, if requested;
- `usage_se` or `usage_var`: posterior uncertainty for usage when available;
- `logit_usage` and `logit_usage_se`: optional logit-scale summaries
  for downstream weighted Gaussian mixed models;
- `effective_n`, `effective_success`: optional pseudo-count summaries
  when supplied by preprocessing or requested through an approximation;
- `diagnostics`: denominator flags, low-information units, and basis metadata.

Candidate PIP should be computed by susieR's standard formula:

```text
PIP_t = 1 - prod_l (1 - alpha_lt).
```

The default credible sets are alpha-based CSs from
`susieR::susie_get_cs(fit, coverage = coverage)` without passing a
dense `X` matrix. This skips purity filtering by default, which is
appropriate for nearby breakpoint candidates. If a caller explicitly
requests purity and the candidate count is below a configured safety
threshold, the implementation may compute a step-basis `Xcorr` and
pass it to `susie_get_cs()`.

For candidate `t` in analysis unit `i`,

```text
A_it = sum_l alpha_lt * mu_lit.
```

where `mu_lit` is the positive-prior posterior mean conditional on
effect `l` selecting candidate `t`. Usage is

```text
U_it = A_it / sum_s A_is.
```

Expected length is

```text
EL_i = sum_t U_it c_t.
```

Balance for cutpoint `r`:

```text
B_i(r) = log((sum_{t: c_t <= r} U_it + eps) /
             (sum_{t: c_t >  r} U_it + eps)).
```

Cutpoints are reporting boundaries, not fitted breakpoints. If
`cutpoints = NULL`, the helper must not require a cutpoint and must
return no balance matrix unless `type = "balance"` was explicitly
requested. If cutpoints are provided, each cutpoint must lie on the
same transcript-oriented coordinate scale as `c_t`; otherwise the
helper must error. The result must keep one column per requested
cutpoint and record the cutpoint values in metadata.

If effective numerator-denominator summaries are supplied or can be
constructed upstream, `apa_phenotype()` may carry them through for
downstream mixed generalized linear models. Otherwise it returns
posterior means, standard errors, and optional precision weights for
weighted Gaussian or logit-scale mixed models.

For the first implementation, uncertainty fields that do not have a
clear implemented approximation may be returned as `NA` with a
diagnostic rather than using an undocumented delta-method shortcut.
The required phenotype is the posterior plug-in usage mean plus the
available posterior summaries from the fitted SuSiE object.

The package must not implement ISSAC-style association testing in this
change. It only exports quantities that such a downstream model could
consume, for example a posterior mean `U_it` with variance `V_it`, a
logit-scale pair `z_it = logit(U_it)` and `se_z_it`, or pseudo-count
fields `effective_success` and `effective_n`.

## Testing strategy

Unit tests must prove correctness module by module before testing the
full wrapper:

1. Step operators equal explicit step matrices for upstream and
   downstream orientations, weighted and unweighted.
2. Half-normal logBF and posterior moments match numerical
   integration on small examples.
3. Weighted cross-unit combiner equals manual weighted sums and
   handles zero weights.
4. Annotation priors are strictly positive, normalized, and do not
   drop candidates.
5. Annotation strength changes prior odds monotonically while
   preserving a positive prior floor.
6. Marginal initialization recovers valid candidate indices on a
   simulated breakpoint example, respects `L`, and does not use a
   dense step matrix.
7. Optional L0 initialization recovers valid candidate indices when
   explicitly requested, respects `L`, discards negative coefficients
   under the positive default, and falls back cleanly to marginal
   initialization when `L0Learn` is unavailable or the explicit-matrix
   threshold is exceeded.
8. `apa_susie()` recovers simulated shared breakpoints with
   unit-specific positive effects.
9. `apa_phenotype()` returns normalized usage and expected length
   with correct dimensions.
10. `apa_phenotype()` handles zero, one, and multiple cutpoints,
   validates cutpoint coordinates, and returns one balance column per
   cutpoint.
11. `apa_phenotype()` exports mixed-model compatible uncertainty fields
   when requested, without running association tests.
12. The default fitting path does not call `apa_step_explicit()`; tests
   may monkey-patch or count calls.

## Simulation validation plan

The implementation is considered statistically credible only after
the following simulations pass. These are not all CRAN tests; some can
be skipped or placed in a vignette/benchmark script, but their
expected outputs must be documented.

1. **No-signal null gene**: simulate corrected coverage with no
   breakpoint. Expected output: low maximum PIP, no confident credible
   set, usage missing or flagged low-information rather than a precise
   false APA phenotype.
2. **Single shared PAS**: simulate one breakpoint shared across units
   with unit-specific positive effects. Expected output: high PIP in a
   small neighborhood of the true candidate and unit-level usage
   increasing with the true effect.
3. **Multiple PAS in one region**: simulate two or more candidate
   breakpoints separated beyond the microheterogeneity merge distance.
   Expected output: the pre-scan retains the region and the final
   model can allocate posterior mass to multiple candidates instead of
   forcing one peak.
4. **Unit sharing benefit**: simulate weak per-unit evidence but a
   common breakpoint across units. Expected output: pooled/shared
   model has higher true-location PIP than fitting each unit
   independently, without forcing equal effect sizes.
5. **Annotation prior stress test**: simulate a strong but imperfect
   annotation. Expected output: higher prior odds and improved weak
   signal recovery when annotation is correct; recovery by data when
   annotation is wrong and the prior floor is not zero.
6. **Bias/offset test**: add a broad positional trend that is supplied
   as an offset. Expected output: false scan peaks are reduced after
   offset subtraction, while true sharp breakpoints remain detectable.
7. **Overlapping-unit weights**: duplicate or partially duplicate
   analysis units. Expected output: downweighting duplicates reduces
   overconfidence in PIPs relative to treating them as independent.
8. **Pre-scan recall**: simulate true PAS not exactly at the local
   maximum. Expected output: expanded candidate regions contain the
   true candidate and local background candidates.

## Code-generation supervision

Coding agents implementing this change must keep the implementation
modular:

- add APA code only in the APA files listed in this spec unless a
  small registration or documentation edit is required elsewhere;
- do not rewrite the existing DWT path;
- do not add `L0Learn` to `Imports`; it should be in `Suggests` and
  guarded by `requireNamespace`;
- keep the core fit and default marginal initializer working without
  `L0Learn`;
- expose diagnostics for every fallback: missing optional package,
  dense-matrix threshold exceeded, failed L0 fit, invalid prior
  weights, and skipped pre-scan positions.

Expected code-review checks:

- every matrix-free operator has an explicit-matrix parity test;
- every probability vector is checked for finite, positive, normalized
  entries;
- every returned fit has finite `alpha`, `mu`, `mu2`, `lbf`, `pip`,
  and phenotype summaries unless explicitly flagged missing;
- the default fitting path can be instrumented to prove
  `apa_step_explicit()` was not called;
- simulation scripts set seeds and record expected qualitative
  outputs, not only absence of errors.

## Performance expectations

For one gene, per-effect step operations must scale as

```text
O(n * (J + T))
```

not

```text
O(n * J * T).
```

The benchmark test should compare operator output to explicit matrix
output on a small example and include a non-CRAN performance smoke on a
larger sparse example. Exact wall-clock thresholds may be loose because CI
is noisy, but the implementation must not materialize dense design
matrices in the default fit path.

Do not impose fixed hard limits on `n`, `J`, or `T` in the default
matrix-free fitting path. Instead, record empirical scaling
diagnostics on the fit object: dimensions, initializer used, IBSS
iteration count, whether an explicit matrix was created, and fallback
reasons. Hard safety thresholds are appropriate only for optional
paths that explicitly build a step matrix, such as
`init_method = "l0learn"`.

## Documentation

Add roxygen documentation for all exported helpers and one vignette or
article stub explaining:

- the unit-level APA model;
- the 0th-order trend-filtering / Haar-like interpretation;
- why the implementation uses a dedicated step transform rather than
  the existing DWT path;
- annotation priors;
- strong annotation via prior odds, `prior_strength`, and `prior_floor`;
- positive effect priors;
- marginal initialization and optional L0Learn initialization;
- unit weights;
- phenotype outputs, including the interpretation of cutpoints;
- uncertainty exports for downstream mixed models, while making clear
  that QTL association is out of scope.

The docs must emphasize that the model is unit-agnostic and that unit
construction is the user's responsibility.
