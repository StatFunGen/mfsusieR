# apa-step-susie Specification

## ADDED Requirements

### Requirement: Unit-agnostic APA model

The package SHALL provide an APA fitting module whose public API and
documentation use analysis units as the response index. The module
SHALL NOT assume how analysis units were constructed.

#### Scenario: Matrix input with analysis units as columns

- **GIVEN** a numeric matrix `Y` with `J` rows and `n` columns
- **AND** ordered positions `pos`
- **AND** candidate PAS positions `candidates`
- **WHEN** `apa_susie(Y, pos, candidates, ...)` is called
- **THEN** the implementation SHALL treat each column as one analysis
  unit
- **AND** map positions to rows and candidates to predictors.

#### Scenario: List input with analysis units as entries

- **GIVEN** a list `Y` whose entries are numeric vectors of length `J`
- **AND** ordered positions `pos`
- **AND** candidate PAS positions `candidates`
- **WHEN** `apa_susie(Y, pos, candidates, ...)` is called
- **THEN** each list entry SHALL be treated as one analysis unit.

### Requirement: Working coverage outcome

The package documentation SHALL define `Y_i(j)` as the working coverage
outcome for analysis unit `i` and ordered position `j`.

#### Scenario: Documentation distinguishes raw and working coverage

- **GIVEN** the APA wrapper documentation
- **THEN** it SHALL state that `Y_i(j)` is not raw coverage
- **AND** it SHALL describe `Y_i(j)` as coverage after library-size
  normalization, bias correction, and optional variance-stabilizing or
  residualizing transformation.

#### Scenario: Bias correction is upstream of the APA fit

- **GIVEN** a call to `apa_susie()`
- **WHEN** `Y`, `offset`, or `weights` are supplied
- **THEN** the wrapper SHALL treat them as the working outcome and
  optional working-model components
- **AND** SHALL NOT require a built-in library-size, GC-bias, or
  5'/3' bias correction pipeline.

#### Scenario: Row precision weights

- **GIVEN** row weights supplied as `NULL`, a length-`J` vector, a
  `J x n` matrix, or a list of `n` length-`J` vectors
- **WHEN** `apa_susie()` validates inputs
- **THEN** it SHALL convert them to one finite non-negative
  length-`J` weight vector per analysis unit
- **AND** reject missing, negative, non-finite, or zero-total weights
  for any unit
- **AND** document that row weights are precision weights within a
  unit and are distinct from cross-unit evidence weights.

#### Scenario: No explicit intercept in IBSS

- **GIVEN** working coverage `Y`, optional `offset`, and optional row
  precision weights
- **WHEN** `apa_susie()` constructs the APA data object
- **THEN** it SHALL subtract the offset, using zero offset when absent
- **AND** compute unit-specific weighted response centers
- **AND** compute weighted candidate step centers
- **AND** run IBSS on centered working coverage and centered step
  predictors
- **AND** SHALL NOT add an explicit intercept column or baseline
  candidate to the IBSS variable set.

#### Scenario: Fitted coverage reconstruction

- **GIVEN** a fitted APA model with centered fitted step contribution
- **WHEN** fitted working coverage is requested
- **THEN** it SHALL reconstruct
  `offset_i(j) + ybar_i + sum_l sum_t alpha_lt * mu_lit * Sc_it(j)`
- **AND** SHALL keep the baseline outside the IBSS variable set.

### Requirement: Step-basis transform module

The package SHALL implement a step-basis module that represents
candidate PAS locations as step functions over ordered positions.

#### Scenario: Upstream step orientation

- **GIVEN** ordered position coordinates `x_j`
- **AND** candidate PAS coordinates `c_t`
- **WHEN** `orientation = "upstream"`
- **THEN** the basis SHALL represent `S[j, t] = 1(x_j <= c_t)`.

#### Scenario: Downstream step orientation

- **GIVEN** ordered position coordinates `x_j`
- **AND** candidate PAS coordinates `c_t`
- **WHEN** `orientation = "downstream"`
- **THEN** the public basis SHALL represent `S[j, t] = 1(x_j > c_t)`
- **AND** any internal transformed implementation SHALL pass explicit
  matrix parity tests for multiplication, transpose multiplication,
  and column norms under this public convention.

#### Scenario: Matrix-free multiplication

- **GIVEN** a step-basis object and coefficient vector `b`
- **WHEN** `apa_step_Xb(basis, b)` is called
- **THEN** the result SHALL equal explicit multiplication
  `S %*% b`
- **AND** the default fitting path SHALL NOT materialize dense `S`.

#### Scenario: Matrix-free transpose multiplication

- **GIVEN** a step-basis object and residual vector `r`
- **WHEN** `apa_step_Xtr(basis, r)` is called
- **THEN** the result SHALL equal explicit multiplication
  `crossprod(S, r)`.

#### Scenario: Weighted transpose multiplication

- **GIVEN** a step-basis object, residual vector `r`, and weights `w`
- **WHEN** `apa_step_Xtr(basis, r, weights = w)` is called
- **THEN** the result SHALL equal explicit multiplication
  `crossprod(S, w * r)`.

#### Scenario: Column norms

- **GIVEN** a step-basis object
- **WHEN** `apa_step_colnorm(basis)` is called
- **THEN** the result SHALL equal `colSums(S * S)` for the explicit
  basis.

#### Scenario: Candidate coordinate mapping

- **GIVEN** candidate coordinates not exactly equal to observed
  position coordinates
- **WHEN** `apa_step_basis(pos, candidates, orientation = "upstream")`
  is constructed
- **THEN** each candidate SHALL be mapped to the last observed
  position with `x_j <= c_t`
- **AND** the original candidate coordinate and mapped position index
  SHALL both be stored.

#### Scenario: Centered step cache

- **GIVEN** row weights for an analysis unit
- **WHEN** `apa_step_center_cache()` is called
- **THEN** it SHALL compute weighted step means `Sbar_it`
- **AND** centered denominators
  `d_it = sum_j w_ij * (S[j,t] - Sbar_it)^2`
- **AND** flag candidates with near-zero denominators or insufficient
  effective weight on either side
- **AND** SHALL compute these quantities without materializing a dense
  `J x T` step matrix.

#### Scenario: Trend-filtering helper parity

- **GIVEN** a full-candidate unweighted step basis
- **WHEN** `apa_step_Xb()`, `apa_step_Xtr()`, and
  `apa_step_colnorm()` are compared with an explicit 0th-order
  changepoint design
- **THEN** the fast helpers SHALL match the explicit design up to the
  documented upstream or downstream orientation convention
- **AND** tests SHALL NOT require materializing the explicit design
  outside bounded-size parity checks.

### Requirement: Trend-filtering and Haar documentation

The package documentation SHALL describe the APA step basis as a
0th-order trend-filtering / changepoint basis and SHALL explain its
relationship to, but distinction from, the existing Haar/DWT path.

#### Scenario: Dedicated transform, not DWT overload

- **GIVEN** user-facing documentation for `apa_susie()`
- **THEN** it SHALL state that the step basis is Haar-like because it
  represents piecewise-constant changes
- **AND** it SHALL state that arbitrary candidate PAS locations are
  not generally single coefficients in a fixed dyadic Haar DWT
- **AND** it SHALL state that APA uses a dedicated step-transform
  module rather than the existing DWT module as the default fit path.

### Requirement: Candidate PAS pre-scan

The package SHALL provide an optional high-recall candidate pre-scan
function that can construct candidate PAS inputs for the APA fitting
module. The pre-scan SHALL be documented as input-variable screening,
not final APA inference.

#### Scenario: Dense default scan

- **GIVEN** a working coverage matrix `Y` with `J` ordered positions
- **AND** no custom scan grid
- **WHEN** `apa_prescan(Y, pos, ...)` is called
- **THEN** the function SHALL scan all observed positions or bins as
  possible breakpoints
- **AND** treat `R = J` and `d_r = x_r` up to boundary positions that
  fail minimum information checks.

#### Scenario: Optional offset and fitted intercept

- **GIVEN** working coverage `Y`
- **AND** `offset = NULL`
- **WHEN** `apa_prescan()` computes marginal scan estimates
- **THEN** it SHALL use zero position-specific offset
- **AND** SHALL fit a weighted intercept in each unit-by-position
  marginal regression by weighted centering.

#### Scenario: Supplied offset

- **GIVEN** working coverage `Y`
- **AND** a finite supplied `offset` with compatible shape
- **WHEN** `apa_prescan()` computes marginal scan estimates
- **THEN** it SHALL subtract the supplied offset before fitting the
  marginal step regression
- **AND** SHALL still fit a weighted intercept by weighted centering
- **AND** SHALL document that the offset represents external
  position-specific nuisance correction, not an average coverage
  level.

#### Scenario: No internal smooth offset fitting

- **GIVEN** a call to `apa_prescan()`
- **THEN** the function SHALL NOT estimate a flexible smooth GC,
  mappability, or 5'/3' bias offset internally
- **AND** SHALL leave those corrections to preprocessing or to a
  supplied fixed offset.

#### Scenario: Scan score definitions

- **GIVEN** a scan position `d_r` and analysis unit `i`
- **WHEN** `apa_prescan()` computes marginal statistics
- **THEN** it SHALL compute `bhat_ir^scan` from the weighted centered
  regression on `1(x_j <= d_r)`
- **AND** compute `shat2_ir` as the working sampling variance of
  `bhat_ir^scan`
- **AND** compute the unit scan score as
  `m_ir = log BF^+(bhat_ir^scan, shat2_ir, tau2_scan)`
- **AND** SHALL NOT describe `shat2_ir` as a score.

#### Scenario: Default scan tau2

- **GIVEN** `tau2_scan = NULL`
- **WHEN** `apa_prescan()` computes scan scores
- **THEN** it SHALL estimate one scan slab scale from scan marginal
  statistics using the documented positive-excess rule
- **AND** record the scan slab scale and fallback status in
  diagnostics.

#### Scenario: Low-information scan positions

- **GIVEN** a scan position with too little effective weight upstream
  or downstream, or a near-zero centered denominator
- **WHEN** `apa_prescan()` scans that position
- **THEN** it SHALL skip or flag that position
- **AND** record the reason in diagnostics.

#### Scenario: Pooled score curve

- **GIVEN** unit-level scan scores
- **WHEN** pooled scores are requested
- **THEN** the function SHALL compute weighted sums across analysis
  units using normalized unit weights.

#### Scenario: Default peak retention

- **GIVEN** a pooled score curve with local maxima
- **WHEN** `apa_prescan()` selects peaks
- **THEN** it SHALL retain finite local maxima with positive pooled
  score by default
- **AND** if `max_peaks` is finite, retain at most `max_peaks`
  highest-scoring retained maxima.

#### Scenario: Candidate regions and local background

- **GIVEN** local maxima in the pooled score curve
- **WHEN** `apa_prescan()` constructs candidate regions
- **THEN** it SHALL expand each retained peak by
  `region_half_width`, defaulting to exactly 100 nt
- **AND** retain fine-scale and flanking local-background candidate
  positions inside expanded regions when requested
- **AND** SHALL NOT keep only the peak center.

#### Scenario: Microheterogeneity merge

- **GIVEN** nearby candidate positions within `merge_distance`
- **WHEN** `apa_prescan()` finalizes candidates
- **THEN** it SHALL collapse them using a documented representative
  rule when `merge_distance` is positive
- **AND** the default merge distance SHALL be exactly 25 nt
- **AND** distinct peaks outside the merge distance SHALL remain
  eligible as separate candidate PAS
- **AND** `merge_distance = NULL` or `merge_distance <= 0` SHALL
  disable this merging step.

#### Scenario: Candidate source union

- **GIVEN** scan-derived regions, annotations, tail-read sites, and
  motif sites
- **WHEN** `apa_prescan()` finalizes the candidate set
- **THEN** it SHALL return the union of these sources
- **AND** return metadata indicating which source or sources supported
  each candidate.

#### Scenario: Pre-scan efficiency

- **GIVEN** `n` analysis units, `J` positions, and dense scan
  `R = J`
- **WHEN** `apa_prescan()` computes marginal scan statistics
- **THEN** it SHALL use cumulative sums or equivalent linear-time
  operations
- **AND** SHALL NOT fit `R` separate dense regressions per unit.

### Requirement: Annotation-informed location prior

The package SHALL support candidate-level prior weights derived from
annotations or user-supplied scores.

#### Scenario: User-supplied scores

- **GIVEN** finite candidate scores
- **WHEN** `apa_prior_from_scores(scores)` is called
- **THEN** the result SHALL be strictly positive
- **AND** sum to one
- **AND** have the same length as `scores`.

#### Scenario: Direct prior weights

- **GIVEN** user-supplied finite `prior_weights`
- **WHEN** the APA wrapper normalizes location priors
- **THEN** the result SHALL be strictly positive and sum to one
- **AND** SHALL respect the same positive prior-floor rules used by
  annotation-derived priors.

#### Scenario: Annotation feature matrix

- **GIVEN** a numeric feature matrix with one row per candidate
- **WHEN** `apa_prior_from_annotations(features, coefficients)` is
  called
- **THEN** the result SHALL be strictly positive and sum to one.

#### Scenario: Annotation coefficients omitted

- **GIVEN** `coefficients = NULL`
- **AND** an annotation feature table without a numeric column named
  `score`
- **WHEN** `apa_prior_from_annotations(features, coefficients = NULL)`
  is called
- **THEN** the helper SHALL throw an informative error
- **AND** SHALL NOT invent feature weights.

#### Scenario: Annotation fixed linear score

- **GIVEN** annotation features and supplied coefficients
- **WHEN** `apa_prior_from_annotations(features, coefficients)` is
  called
- **THEN** coefficients SHALL either be named to match feature columns
  or have length equal to the number of feature columns
- **AND** the score SHALL be computed as
  `intercept + features %*% coefficients`.

#### Scenario: No trained annotation prediction model

- **GIVEN** annotation features for candidates
- **WHEN** first-version APA prior helpers are used
- **THEN** they SHALL use direct prior weights, direct numeric scores,
  or a supplied fixed linear score
- **AND** SHALL NOT train a genome-wide prediction model or
  classifier.

#### Scenario: Strong annotation through prior odds

- **GIVEN** finite candidate scores with at least two distinct values
- **WHEN** `apa_prior_from_scores(scores, prior_strength = a)` and
  `apa_prior_from_scores(scores, prior_strength = b)` are called with
  `b > a > 0`
- **THEN** the prior odds comparing a higher-scored candidate to a
  lower-scored candidate SHALL be no smaller under strength `b` than
  under strength `a`.

#### Scenario: Positive prior floor

- **GIVEN** finite annotation scores
- **WHEN** a positive `prior_floor` is supplied
- **THEN** every returned candidate prior probability SHALL be at
  least `prior_floor` up to numerical tolerance
- **AND** invalid floors, including negative values or
  `T * prior_floor >= 1`, SHALL produce an informative error.

#### Scenario: Annotation is not hard filtering

- **GIVEN** finite annotation features for all candidates
- **WHEN** an annotation prior is computed
- **THEN** every candidate SHALL retain positive prior probability
  unless the user supplied invalid input and requested failure.

#### Scenario: Annotation does not encode sign

- **GIVEN** annotation features
- **WHEN** the APA wrapper computes prior weights
- **THEN** no annotation-derived value SHALL be interpreted as a
  positive or negative effect direction.

#### Scenario: Annotation does not force unit usage

- **GIVEN** an annotation-informed location prior with strong prior
  odds for one candidate
- **WHEN** the model is fit
- **THEN** the annotation SHALL affect the shared location prior
  conditional on an effect being active
- **AND** SHALL NOT require every analysis unit to have nonzero usage
  of that candidate.

### Requirement: One-sided effect prior

The package SHALL implement a one-sided Gaussian prior for step-basis
effect magnitudes. Effect magnitudes SHALL be documented as `beta`,
while `alpha` SHALL refer to posterior location probabilities.

#### Scenario: Positive-prior log Bayes factor

- **GIVEN** finite marginal estimates `bhat`, sampling variances `shat2`,
  and positive scalar `tau2`
- **WHEN** `apa_halfnormal_lbf(bhat, shat2, tau2)` is called
- **THEN** it SHALL return the positive-prior log Bayes factor
  `log BF^+ = log BF^N + log(2 Phi(mu / sqrt(v)))`
- **AND** all returned values SHALL be finite for finite inputs.

#### Scenario: Positive-prior moments

- **GIVEN** finite marginal estimates `bhat`, sampling variances `shat2`,
  and positive scalar `tau2`
- **WHEN** `apa_halfnormal_moments(bhat, shat2, tau2)` is called
- **THEN** it SHALL return posterior means and second moments for
  the normal posterior truncated to `(0, Inf)`
- **AND** posterior means SHALL be non-negative
- **AND** second moments SHALL be no smaller than squared means up
  to numerical tolerance.

### Requirement: Shared PAS locations with unit-specific effects

The package SHALL share candidate PAS locations across analysis
units while allowing effect sizes to differ by unit.

#### Scenario: Shared alpha

- **GIVEN** a fitted APA model with `L` effects and `T` candidates
- **THEN** `alpha` SHALL have shape `L x T`
- **AND** the same `alpha[l, t]` SHALL be used for all analysis units.

#### Scenario: Unit-specific beta moments

- **GIVEN** a fitted APA model with `n` analysis units
- **THEN** posterior effect moments SHALL retain a unit dimension
- **AND** unit-specific `beta` effects SHALL be recoverable for
  phenotype construction.

### Requirement: Weighted cross-unit evidence

The package SHALL provide a cross-unit logBF combiner that supports
effective unit weights.

#### Scenario: Default weights

- **GIVEN** no user-supplied unit weights
- **WHEN** per-unit logBFs are combined
- **THEN** the combiner SHALL use weight one for every unit.

#### Scenario: Weighted combination

- **GIVEN** unit weights `w_i`
- **AND** per-unit logBF vectors `lbf_i`
- **WHEN** the weighted combiner is called
- **THEN** it SHALL normalize positive weights to mean one, leaving
  zeros as zeros
- **AND** return `sum_i w_i_norm * lbf_i`.

#### Scenario: Global weight scale invariance

- **GIVEN** two positive unit-weight vectors that differ only by a
  positive constant multiplier
- **WHEN** per-unit logBFs are combined
- **THEN** the combined logBF SHALL be identical up to numerical
  tolerance.

#### Scenario: Invalid weights

- **GIVEN** negative, missing, or non-finite unit weights
- **WHEN** the weighted combiner is constructed or called
- **THEN** the implementation SHALL throw an informative error.

### Requirement: Initialization

The package SHALL provide a matrix-free marginal warm start for the
APA step-SuSiE fit. Optional L0Learn-based initialization SHALL be
available only when explicitly requested. No initializer SHALL change
the posterior model being optimized.

#### Scenario: Default marginal initialization

- **GIVEN** `init_method = "marginal"` or the default wrapper call
- **WHEN** `apa_susie()` initializes the fit
- **THEN** it SHALL compute marginal positive step estimates and
  sampling variances using centered working coverage, cached centered
  step denominators, and matrix-free step operators
- **AND** it SHALL pool candidate support across analysis units using
  normalized unit weights
- **AND** it SHALL convert at most `L` retained candidates into a
  SuSiE-compatible `model_init`
- **AND** it SHALL NOT call `apa_step_explicit()` or construct a
  dense or sparse explicit step matrix.

#### Scenario: Successful L0Learn initialization

- **GIVEN** `init_method = "l0learn"`
- **AND** `L0Learn` is installed
- **AND** the explicit sparse step matrix is below the configured
  safety threshold
- **AND** row weights are absent or constant within every analysis
  unit
- **WHEN** `apa_susie()` initializes the fit
- **THEN** it SHALL run `L0Learn::L0Learn.cvfit(S, ytilde_i,
  penalty = "L0", ...)` for each analysis unit
- **AND** rely on L0Learn's intercept and remove that intercept after
  fitting
- **AND** select candidate PAS with nonzero coefficients
- **AND** convert at most `L` retained candidates into a
  SuSiE-compatible `model_init`.

#### Scenario: L0Learn non-constant row-weight fallback

- **GIVEN** `init_method = "l0learn"`
- **AND** row weights are non-constant within at least one analysis
  unit
- **WHEN** `apa_susie()` initializes the fit
- **THEN** it SHALL skip L0Learn initialization
- **AND** fall back to marginal initialization with a diagnostic.

#### Scenario: SuSiE-shaped warm start

- **GIVEN** selected candidates `t_l` and unit-specific starting
  effects `beta_init[i,l]`
- **WHEN** the initializer builds `model_init`
- **THEN** `alpha` SHALL be an `L x T` matrix with one active
  candidate per initialized effect
- **AND** `mu` and `mu2` SHALL use the mfsusieR list-of-list shape
  with `mu[[l]][[i]][t_l, 1] = beta_init[i,l]` and
  `mu2[[l]][[i]][t_l, 1] = beta_init[i,l]^2`
- **AND** all finite-shape, dimension, and probability checks SHALL
  pass before the object is passed to the IBSS workhorse.

#### Scenario: Positive initialization under upstream model

- **GIVEN** the default upstream-positive parameterization
- **AND** an initializer returns a non-positive coefficient for a candidate
- **WHEN** the initializer builds retained candidates
- **THEN** the initializer SHALL discard the non-positive coefficient
- **AND** SHALL NOT initialize a positive effect with a negative
  `beta`.

#### Scenario: Fixed tau2 in first implementation

- **GIVEN** `tau2 = NULL`
- **AND** marginal estimates `bhat_it` and sampling variances
  `shat2_it` are available from the matrix-free marginal initializer
- **WHEN** `apa_susie()` prepares the model
- **THEN** it SHALL compute
  `e_it = max(bhat_it, 0)^2 - shat2_it`
- **AND** keep finite positive excess estimates
- **AND** initialize one positive scalar `tau2` before the IBSS loop
  as a robust weighted quantile of those excess estimates, defaulting
  to the weighted median with normalized unit weights
- **AND** enforce a positive `tau2_floor` on the working-outcome scale
- **AND** define the default working-outcome scale as the median
  across units of finite weighted variances of centered `Y - offset`
- **AND** record `tau2`, `tau2_source`, `tau2_floor`, the quantile
  used, and the retained-excess count in diagnostics
- **AND** the APA prior-variance hooks SHALL leave `tau2` unchanged
  during IBSS
- **AND** `update_model_variance.mf_apa_individual()` SHALL remain
  responsible for residual-variance updates and derived-cache refresh
  when `estimate_residual_variance = TRUE`
- **AND** the implementation SHALL NOT perform grid, mixture, or
  in-loop `tau2` estimation in this change.

#### Scenario: User-supplied tau2

- **GIVEN** the user supplies a positive scalar `tau2`
- **WHEN** `apa_susie()` prepares the model
- **THEN** it SHALL use that value unchanged as the positive slab
  scale
- **AND** record `tau2_source = "user"`
- **AND** the APA prior-variance hooks SHALL leave `tau2` unchanged
  during IBSS.

#### Scenario: Tau2 floor fallback

- **GIVEN** `tau2 = NULL`
- **AND** no finite positive excess estimates are available
- **WHEN** `apa_susie()` prepares the model
- **THEN** it SHALL set `tau2 = tau2_floor`
- **AND** record `tau2_source = "floor_no_positive_excess"`
- **AND** report the fallback in diagnostics.

#### Scenario: Residual variance initialization

- **GIVEN** `residual_variance = NULL`
- **WHEN** `apa_susie()` prepares the model
- **THEN** it SHALL initialize each unit's residual variance from a
  MAD-difference estimate on centered working coverage when finite
  and positive
- **AND** fall back to weighted variance of centered working coverage
  when the MAD estimate is zero or non-finite
- **AND** use a positive floor with diagnostics when neither estimate
  is available.

#### Scenario: Residual variance update

- **GIVEN** `estimate_residual_variance = TRUE`
- **WHEN** `update_model_variance.mf_apa_individual()` is called
- **THEN** it SHALL update each unit's residual variance as
  `max(sigma2_floor_i, ERSS_i / sum_j w_ij)`
- **AND** `ERSS_i` SHALL be the weighted expected squared residual
  including posterior second-moment terms
- **AND** it SHALL NOT use only the squared posterior mean residual
- **AND** it SHALL leave the fixed slab scale `tau2` unchanged.

#### Scenario: Candidate validation and L cap

- **GIVEN** candidate validation leaves `T` informative candidates
- **WHEN** `apa_susie()` prepares the model
- **THEN** it SHALL set the effective number of effects to
  `min(L, T)` and record any adjustment in diagnostics.

#### Scenario: No informative candidates

- **GIVEN** no candidate maps to an informative step after validation
- **WHEN** `apa_susie()` prepares the model
- **THEN** it SHALL stop before IBSS with an informative error.

#### Scenario: Candidate union across units

- **GIVEN** different analysis units select different candidates in
  their marginal or optional L0 initializers
- **WHEN** the initializer combines them
- **THEN** it SHALL rank the union of candidates by a documented
  positive-evidence score, such as
  `sum_i unit_weights_norm[i] * logBF_it^+` for marginal initialization
  or `sum_i unit_weights_norm[i] * max(beta_hat_it, 0)^2` for L0
  initialization
- **AND** retain at most `min(L, max_nonzero, T)` candidates.

#### Scenario: L0Learn fallback

- **GIVEN** `init_method = "l0learn"`
- **AND** `L0Learn` is unavailable, the explicit matrix threshold is
  exceeded, the fit fails, or no valid positive candidate is selected
- **WHEN** `apa_susie()` initializes the fit
- **THEN** it SHALL fall back to marginal initialization
- **AND** record an informative diagnostic on the fit object.

### Requirement: APA data subclass

The package SHALL implement an APA-specific mfsusieR data subclass
that reuses inherited mfsusieR behavior where valid and overrides only
APA-specific computations.

#### Scenario: No dense design in default path

- **GIVEN** a normal call to `apa_susie()` with default marginal
  initialization
- **WHEN** the model is fit
- **THEN** dense `J x T` step design construction SHALL NOT be used
  by the fitting loop.

#### Scenario: S3 methods registered

- **GIVEN** the package has loaded
- **WHEN** an APA data object is passed to susieR generics used by the
  IBSS loop
- **THEN** APA-specific methods SHALL be dispatched for every method
  whose default implementation assumes dense regression design,
  wavelet coefficient arrays, or symmetric normal effects
- **AND** APA-specific prior hooks SHALL be registered through the
  same roxygen/NAMESPACE S3 mechanism used by the corresponding
  `mf_individual` hooks.

#### Scenario: Centered SER statistics

- **GIVEN** centered partial residuals for one effect
- **WHEN** `compute_ser_statistics.mf_apa_individual()` is called
- **THEN** it SHALL compute
  `x_it = sum_j w_ij * Sc_it(j) * Rc_il(j)`
- **AND** compute `bhat_it = x_it / d_it`
- **AND** compute `shat2_it = sigma2_i / d_it`
- **AND** use cached `d_it` from the centered step cache
- **AND** make low-information unit-candidate pairs uninformative.

### Requirement: APA phenotype extraction

The package SHALL provide phenotype extraction that converts posterior
step-basis effects into unit-level APA summaries.

#### Scenario: Candidate PIP

- **GIVEN** a fitted APA model with posterior `alpha`
- **WHEN** candidate PIPs are requested
- **THEN** the implementation SHALL use susieR posterior semantics,
  `PIP_t = 1 - prod_l (1 - alpha_lt)`.

#### Scenario: Default credible sets

- **GIVEN** a fitted APA model
- **WHEN** credible sets are requested with default options
- **THEN** the implementation SHALL construct alpha-based credible
  sets using susieR's `susie_get_cs()` semantics
- **AND** SHALL NOT materialize dense `S` or require correlation
  purity filtering by default.

#### Scenario: Credible-set diffuseness diagnostics

- **GIVEN** default credible sets are constructed
- **WHEN** credible-set diagnostics are returned
- **THEN** each credible set SHALL include size, `cs_fraction = size/T`,
  claimed coverage, and `cs_is_diffuse`
- **AND** the default diffuse rule SHALL be `cs_fraction > 0.5`
- **AND** diffuse credible sets SHALL be reported but SHALL NOT be
  treated as confident breakpoint calls.

#### Scenario: Optional bounded purity filtering

- **GIVEN** a user explicitly requests credible-set purity filtering
- **WHEN** candidate count is below the configured explicit-correlation
  safety threshold
- **THEN** the implementation MAY compute a step-basis `Xcorr` and
  pass it to `susie_get_cs()`
- **AND** SHALL otherwise skip purity filtering with a diagnostic.

#### Scenario: Candidate usage

- **GIVEN** a fitted APA model
- **WHEN** `apa_phenotype(fit, type = "usage")` is called
- **THEN** it SHALL return a unit-by-candidate usage matrix
- **AND** each informative unit's usage SHALL sum to one within
  numerical tolerance
- **AND** usage SHALL be normalized only over modeled step-derived
  candidate contributions
- **AND** the centered baseline/intercept component SHALL NOT be
  included as a distal PAS or in the usage denominator.

#### Scenario: Usage reportability

- **GIVEN** a fitted APA model
- **WHEN** `apa_phenotype()` returns usage
- **THEN** it SHALL include a unit-level `usage_reportable` flag
- **AND** the default rule SHALL require usage denominator at least
  `eps`, `max(candidate_pip) >= 0.5`, and at least one non-diffuse
  credible set
- **AND** non-reportable usage SHALL be returned as missing or flagged
  according to the documented output contract.

#### Scenario: Expected length

- **GIVEN** candidate coordinates
- **AND** candidate usage estimates
- **WHEN** expected length is requested
- **THEN** the helper SHALL return `sum_t usage_it * coordinate_t`
  for each analysis unit with reportable usage
- **AND** expected length SHALL be missing with diagnostics for
  non-reportable units.

#### Scenario: Cutpoints are optional reporting boundaries

- **GIVEN** a fitted APA model
- **AND** `cutpoints = NULL`
- **WHEN** usage, PIPs, credible sets, or expected length are
  requested
- **THEN** `apa_phenotype()` SHALL NOT require a proximal-distal
  cutpoint
- **AND** SHALL NOT return a balance matrix unless balance was
  explicitly requested.

#### Scenario: Proximal-distal balance

- **GIVEN** one or more cutpoints
- **WHEN** balance is requested
- **THEN** the helper SHALL return
  `log((proximal_usage + eps) / (distal_usage + eps))` for each
  analysis unit and cutpoint
- **AND** SHALL keep one output column per requested cutpoint
- **AND** SHALL record cutpoint values in output metadata.

#### Scenario: Invalid cutpoint coordinates

- **GIVEN** a cutpoint outside the coordinate scale used by candidate
  PAS positions
- **WHEN** balance is requested
- **THEN** the helper SHALL throw an informative error.

#### Scenario: Low-information unit

- **GIVEN** a unit with near-zero total posterior effect contribution
- **WHEN** usage is requested
- **THEN** usage for that unit SHALL be returned as missing
- **AND** a diagnostic flag SHALL be recorded.

#### Scenario: Mixed-model compatible output

- **GIVEN** effective numerator-denominator summaries or precision
  weights supplied by preprocessing
- **WHEN** `apa_phenotype(export_for_mixed_model = TRUE)` returns
  output
- **THEN** it SHALL preserve these summaries or weights where
  available
- **AND** SHALL return posterior uncertainty summaries such as
  `usage_var` or `usage_se`, and logit-scale summaries where
  numerically defined
- **AND** SHALL NOT run QTL, ISSAC, or other association testing.

#### Scenario: Unimplemented uncertainty summaries

- **GIVEN** an uncertainty field is requested but no documented
  approximation is implemented
- **WHEN** `apa_phenotype()` returns output
- **THEN** the field SHALL be returned as missing
- **AND** a diagnostic SHALL explain that the approximation is not
  implemented.

### Requirement: Performance scaling

The APA step-SuSiE fitting path SHALL scale linearly in the number of
positions plus the number of candidates per analysis unit for step
operations.

#### Scenario: Per-effect operator scaling

- **GIVEN** `n` analysis units, `J` positions, and `T` candidates
- **WHEN** one SER sweep requires step-basis multiplication and
  transpose multiplication
- **THEN** the intended cost SHALL be `O(n * (J + T))`
- **AND** SHALL NOT require `O(n * J * T)` dense multiplication in
  the default path.

#### Scenario: Scaling diagnostics without fixed size limits

- **GIVEN** any fitted APA model
- **WHEN** the fit object is returned
- **THEN** it SHALL record fit dimensions such as `n`, `J`, `T`, `L`,
  the number of IBSS iterations, and whether any explicit step matrix
  was created
- **AND** the default matrix-free fitting path SHALL NOT impose fixed
  hard limits on the number of analysis units, positions, or
  candidates
- **AND** fixed safety thresholds SHALL apply only to optional paths
  that materialize explicit step matrices, such as `init_method =
  "l0learn"`.

### Requirement: Validation and simulation evidence

The package SHALL include deterministic tests and documented
simulation checks that demonstrate correctness of the APA model and
its candidate-screening path.

#### Scenario: Operator parity tests

- **GIVEN** a small step-basis example
- **WHEN** matrix-free operators and explicit step matrices are both
  evaluated
- **THEN** `apa_step_Xb`, `apa_step_Xtr`, and `apa_step_colnorm`
  SHALL match explicit multiplication within numerical tolerance.

#### Scenario: Centered intercept-equivalence tests

- **GIVEN** a small weighted example
- **WHEN** centered APA marginal statistics are compared with an
  explicit weighted least-squares regression containing an intercept
- **THEN** `bhat` and `shat2` SHALL match within numerical tolerance.

#### Scenario: Row-weight validation tests

- **GIVEN** row precision weights supplied as `NULL`, a length-`J`
  vector, a `J x n` matrix, or a list of `n` length-`J` vectors
- **WHEN** the APA wrapper validates input
- **THEN** all valid shapes SHALL be converted to the internal
  per-unit representation
- **AND** missing, negative, non-finite, or zero-total weights for any
  analysis unit SHALL raise informative errors.

#### Scenario: Half-normal numerical tests

- **GIVEN** small finite values of `bhat`, `shat2`, and `tau2`
- **WHEN** half-normal logBFs and moments are computed
- **THEN** they SHALL match numerical integration within documented
  tolerance.

#### Scenario: Residual-variance update tests

- **GIVEN** a small weighted APA fit with fixed posterior
  responsibilities and moments
- **WHEN** residual variance is updated
- **THEN** the implementation SHALL match an explicit weighted ERSS
  calculation that includes posterior second-moment terms
- **AND** SHALL apply the documented per-unit variance floor.

#### Scenario: Candidate validation tests

- **GIVEN** more single effects than informative candidate PAS
- **WHEN** the wrapper validates candidates
- **THEN** it SHALL reduce `L` to the number of informative candidates
  and record this diagnostic
- **AND** if no informative candidate remains, it SHALL stop before
  IBSS with an informative error.

#### Scenario: Fitted reconstruction tests

- **GIVEN** a small centered APA fit with known offset, stored
  baseline, and posterior mean step effects
- **WHEN** fitted working coverage is reconstructed
- **THEN** the result SHALL equal `offset_i + ybar_i +
  sum_l sum_t alpha_lt * mu_lit * Sc_it`
- **AND** the baseline SHALL remain outside the IBSS variable set.

#### Scenario: No-signal simulation

- **GIVEN** simulated working coverage with no true breakpoint
- **WHEN** `apa_susie()` is fit
- **THEN** the model SHALL not return a high-confidence false PAS
  under the documented default thresholds
- **AND** low-information usage SHALL be flagged rather than reported
  as an over-precise phenotype.

#### Scenario: Phenotype reportability tests

- **GIVEN** fitted models with low maximum PIP, diffuse credible sets,
  or too little modeled step-derived abundance
- **WHEN** `apa_phenotype()` reports usage
- **THEN** `usage_reportable` SHALL be `FALSE`
- **AND** expected length SHALL be missing for non-reportable units
- **AND** the centered baseline/intercept SHALL NOT be included in the
  usage denominator.

#### Scenario: Shared-breakpoint simulation

- **GIVEN** simulated working coverage with one or more true shared
  PAS locations and unit-specific positive effects
- **WHEN** `apa_susie()` is fit
- **THEN** candidate PIPs SHALL concentrate near the true candidate
  locations
- **AND** estimated unit-level effects SHALL vary across analysis
  units.

#### Scenario: Pre-scan recall simulation

- **GIVEN** simulated multi-PAS coverage and a dense pre-scan
- **WHEN** candidate regions are constructed
- **THEN** the final candidate set SHALL include the true PAS
  locations or their candidate-resolution neighborhoods
- **AND** include local background candidates inside scan regions.

#### Scenario: Generated-code safeguards

- **GIVEN** the default APA fitting path
- **WHEN** tests instrument `apa_step_explicit()`
- **THEN** the default fit SHALL fail the test if it materializes the
  explicit `J x T` design
- **AND** optional dense construction SHALL be allowed only in tests,
  debugging, or bounded-size L0 initialization.

#### Scenario: Optional L0Learn fallback tests

- **GIVEN** `init_method = "l0learn"` and at least one analysis unit
  has non-constant row precision weights
- **WHEN** initialization is requested
- **THEN** the implementation SHALL skip L0Learn
- **AND** SHALL fall back to matrix-free marginal initialization with
  a diagnostic.

#### Scenario: PIP and credible-set parity tests

- **GIVEN** a fitted APA model with valid `alpha`
- **WHEN** APA PIPs and default credible sets are computed
- **THEN** PIPs SHALL match `susieR::susie_get_pip()`
- **AND** default credible sets SHALL match `susieR::susie_get_cs()`
  called without `X` or `Xcorr`.
