# Add a unit-level APA step-SuSiE module

## Why

mfsusieR already has the right architecture for modular transformed-space
regression: a data/basis construction layer, S3-dispatched SER updates, a
cross-outcome log-BF combiner, and post-fit summarizers. APA phenotype
construction needs the same modular structure, but the default transform is
not the existing DWT transform. It is a step-function / 0th-order
trend-filtering basis over 3'UTR positions, with candidate PAS locations as
variables. This basis is Haar-like because it represents piecewise-constant
coverage changes, but arbitrary PAS candidates are not generally single
coefficients in a fixed dyadic Haar DWT. The APA step basis should therefore
be its own transform module, analogous to the DWT module in architecture.

The goal is to add a small, self-contained APA module without turning
the existing general `mfsusie()` API into an APA-specific interface.
The model must be unit-agnostic: each response is an analysis unit, and
the package must not assume how that unit was produced. A unit may be a
bulk sample, an aggregated sample, a pseudobulk, a metaprofile, or any
other user-supplied coverage profile. Unit construction is outside the
core model.

## What changes

Add an additive APA module under separate files:

- `R/apa_prescan.R`
- `R/apa_step_basis.R`
- `R/apa_prior_annotation.R`
- `R/apa_prior_positive.R`
- `R/apa_cross_unit.R`
- `R/apa_l0_init.R`
- `R/apa_phenotype.R`
- `R/apa_mfsusie.R`

The public-facing entry point is `apa_susie()`. It constructs a
specialized data class and then calls the existing susieR/mfsusieR
IBSS machinery.

The fitted model is a unit-level multi-response step-SuSiE model. For
one gene, with the gene subscript dropped:

```text
Y_i(j) = offset_i(j) + sum_l beta_il S[j, Z_l] + e_i(j),
         beta_il >= 0.
```

where `i` indexes analysis units, `j` indexes observed ordered
positions/bins, `t` indexes candidate PAS locations, `Z_l` is the PAS
location for single-effect component `l`, and `beta_il` is the
unit-specific effect magnitude. The working outcome `Y_i(j)` is
corrected coverage after library-size normalization, bias correction,
and any variance-stabilizing or residualizing transformation. The APA
module starts from this working outcome; it may accept offsets,
weights, and user-supplied bias-corrected coverage, but it does not
define a unique preprocessing pipeline.

The first implementation deliberately avoids grid priors, mixture
priors, and learned scale mixtures. It uses a simple one-sided
Gaussian slab in the step-basis coefficient space:

```text
beta_il | Z_l = t ~ N_+(0, tau2)
```

or a user-supplied scalar `tau2`. Later empirical-Bayes scale updates
can be added as a separate change.

## Codebase review conclusions

The implementation should borrow ideas from two existing code paths
without copying their full complexity:

- **susieR trend filtering**: use the same principle as
  `susie_trendfilter(order = 0)`: represent the design by attributes
  and cumulative-sum algebra rather than by a dense matrix. The APA
  version differs because candidates may be a subset of positions,
  coefficients are constrained positive by orientation, and there are
  multiple analysis units with shared locations. The concrete code
  patterns to mirror are `compute_tf_Xb()`, `compute_tf_Xty()`,
  `compute_tf_d()`, `compute_colstats()`, `compute_Xb()`,
  `compute_Xty()`, the MAD residual-variance warm-start in
  `susie_trendfilter()`, and the `model_init` conventions used by
  `susie_init_coef()` and `ibss_initialize()`.
- **mfsusieR**: reuse the S3-dispatched IBSS backbone, the
  transform-module architecture, `prior_weights`, warm-start plumbing,
  and the cross-outcome log-BF combiner pattern. Do not route APA
  through the DWT transform.

## Model modules

### Candidate PAS pre-scan module

`R/apa_prescan.R` defines an optional high-recall screening function
that constructs candidate PAS inputs before the APA step-SuSiE fit:

```text
apa_prescan(Y, pos, offset = NULL, weights = NULL,
            grid = NULL, unit_weights = NULL,
            region_half_width = 100, merge_distance = 25,
            annotations = NULL, tail_sites = NULL,
            motif_sites = NULL, ...)
```

The pre-scan is not final inference and must not be treated as a hard
APA call. By default it scans every observed position or bin, so
`R = J` and `d_r = x_r`. A coarser custom grid is allowed as a speed
option.

For each scan position it fits a one-step weighted regression with an
intercept:

```text
Y_i(j) - offset_i(j) = gamma_ir + beta_ir^scan 1(x_j <= d_r) + e_i(j),
beta_ir^scan >= 0.
```

If no position-specific offset is supplied, the offset is zero. The
intercept is always fit by weighted centering and absorbs the
unit-specific average coverage level; users do not need to provide an
average offset. A supplied offset should represent an external
position-specific nuisance correction, such as GC, mappability, or
broad 5'/3' bias. The pre-scan must not estimate a flexible smooth
offset internally because that can remove real breakpoint signal.

The scan computes a unit-level log-Bayes-factor score
`m_ir = log BF^+(bhat_ir^scan, shat2_ir, tau2_scan)`, where
`shat2_ir` is the working sampling variance of `bhat_ir^scan`, not a
score. The default pooled curve is the weighted sum of unit scores.
Candidate regions are seeded from local maxima, expanded by a fixed
half-window `h`, and converted into candidate positions. The default
`h = 100` nt is a permissive region-building default, not a
cleavage-site precision claim. Nearby positions within
`merge_distance`, default exactly 25 nt, are collapsed by a documented
representative rule as microheterogeneity unless users disable merging
with `merge_distance = NULL` or `merge_distance <= 0`.

The returned candidate set must be the union of scan-window
positions, annotation-supported sites, tail-read sites, motif sites,
and fine-scale local background positions inside scan regions. Within
each scan-expanded region, keep flanking low-score positions, not
only the local maximum, so the final model retains nearby competing or
null-like variables.

### Step-basis module

`R/apa_step_basis.R` defines a basis object and fast matrix-free
operations for upstream step functions:

```text
S[j, t] = 1(x_j <= c_t)
```

and, optionally, downstream step functions:

```text
S[j, t] = 1(x_j > c_t).
```

The implementation must not materialize a dense `J x T` matrix in the
default path. It must provide cumulative-sum operators:

- `apa_step_Xb(basis, b)`
- `apa_step_Xtr(basis, r, weights = NULL)`
- `apa_step_colnorm(basis, weights = NULL)`
- `apa_step_explicit(basis, sparse = TRUE)` for tests, debugging, and
  optional initialization only

For upstream steps:

```text
(S' r)_t = sum_{j: x_j <= c_t} r_j
(S b)_j = sum_{t: c_t >= x_j} b_t
S_t' S_t = number of positions with x_j <= c_t
```

Weighted versions replace each summand by `weights[j] * ...`.

### Annotation prior module

`R/apa_prior_annotation.R` converts candidate PAS evidence to
location prior weights:

```text
omega_t = h_eta(F_t)
pi_t    = omega_t / sum_s omega_s.
```

`F_t` may include known PAS annotations, motif evidence, direct
tail-read support, cleavage cluster support, internal priming risk,
mappability, local sequence covariates, or user-defined numeric
features. These features are prior evidence only. The module must not
hard-filter candidates by annotation and must not encode effect signs
in annotation features.

The first implementation accepts either user-supplied `prior_weights`
or a simple fixed score from user-supplied or default feature weights.
It exposes a documented annotation-strength parameter so strong
annotations can produce strong prior odds without hard-filtering other
candidates. The APA fit treats these weights as fixed external inputs.

### Positive prior module

`R/apa_prior_positive.R` provides closed-form single-effect updates
for a one-sided Gaussian slab. Given

```text
bhat_t | beta_t ~ N(beta_t, shat2_t),    beta_t ~ N_+(0, tau2),
```

define

```text
mu_t = tau2 / (tau2 + shat2_t) * bhat_t
v_t  = tau2 * shat2_t / (tau2 + shat2_t).
```

The positive-prior log Bayes factor is

```text
log BF_t^+ = log BF_t^N + log(2 Phi(mu_t / sqrt(v_t))).
```

The posterior moments are those of
`N(mu_t, v_t)` truncated to `(0, Inf)`. The module must expose
helpers for logBF, posterior mean, and posterior second moment.

For a downstream-step orientation, the wrapper may internally sign-flip
the marginal estimates and use the same positive-prior functions. The
annotation prior remains sign-free.

### L0 initialization module

`R/apa_l0_init.R` provides an optional warm start, following susieR's
`L0Learn` initialization pattern. This module does not define the
posterior target. It only provides `model_init` for the main APA
step-SuSiE fit.

The implementation should fit an L0-regularized regression on the same
step-basis predictors, extract nonzero candidate PAS coefficients, and
convert them to a SuSiE-style initialization:

```text
apa_l0learn_init(Y, basis, L, unit_weights = NULL,
                 positive = TRUE, lambda_choice = "cv_min",
                 max_nonzero = L, ...)
```

For each analysis unit `i`, run `L0Learn::L0Learn.cvfit(S_w, y_w,
penalty = "L0")` only when the explicit sparse step matrix size is
below a documented threshold. Here `y_w = sqrt(w_i) * (Y_i - offset_i)`
and `S_w = sqrt(w_i) * S` when row weights are supplied; otherwise
`y_w = Y_i - offset_i` and `S_w = S`. Select the lambda by the
requested rule, remove the intercept, and keep nonzero coefficients.
Under the default positive upstream parameterization, discard
non-positive coefficients. Combine candidates selected across units by
the weighted support score
`score_t = sum_i unit_weights[i] * max(beta_hat_it, 0)^2`, keep at
most `L`, and create one one-hot SuSiE effect per retained candidate.
For selected candidates missing in a unit, initialize that unit's
coefficient by the positive marginal least-squares estimate at that
candidate, or zero if there is no positive marginal evidence.

If `L0Learn` is not installed, the sparse explicit matrix would exceed
the safety threshold, no positive candidate is selected, or the
L0Learn fit fails, the wrapper must fall back to the standard
initialization with a clear diagnostic.

### Cross-unit sharing module

`R/apa_cross_unit.R` adds a new `combine_outcome_lbfs()` method for
weighted analysis units:

```text
log BF_t^joint = sum_i w_i log BF_it.
```

The existing `prior_weights` mechanism supplies `log pi_t` through
the usual SuSiE softmax path. The weights `w_i` capture effective
information content and dependency among user-supplied units. If not
provided, all weights default to 1.

### APA phenotype module

`R/apa_phenotype.R` converts the posterior on step-basis effects into
unit-level APA phenotypes. It must output model-derived quantities,
not only legacy PDUI-like ratios:

- candidate-level posterior PAS inclusion probabilities;
- credible sets for PAS / breakpoints;
- unit-by-candidate posterior expected usage;
- proximal-distal balance for user-specified reporting cutpoints;
- expected 3'UTR length;
- posterior uncertainty summaries where available;
- optional uncertainty and effective-count fields that make the output
  compatible with downstream mixed-model software, without implementing
  the association model in this change.

For candidate usage:

```text
A_it = sum_l alpha_lt E[beta_il | Z_l = t, Y]
U_it = A_it / sum_s A_is.
```

If the denominator is below a numerical tolerance, the usage for that
unit/gene is returned as missing with a diagnostic flag.

Cutpoints are not model parameters. They are post-processing reporting
boundaries used to collapse a multi-PAS usage vector into a
proximal-versus-distal balance. If no cutpoint is supplied, the module
must still return usage, PIPs, credible sets, and expected length.

## Impact

- New specs: `inst/openspec/changes/add-unit-level-apa-step-susie/specs/apa-step-susie/spec.md`.
- New code files only; existing mfsusieR behavior is unchanged unless
  users call the new APA wrapper.
- Existing wavelet/DWT paths remain untouched.
- New S3 registrations are needed only for the APA data subclass and
  the weighted cross-unit combiner.

## Out of scope

- Any QTL association model.
- ISSAC-style measurement-error or uncertainty-aware association
  extensions.
- A built-in library-size, GC-bias, or 5'/3' bias correction pipeline.
- Scale-mixture priors, grid priors, or full EB mixture updates.
- A dense-matrix-only implementation as the default path.
- Any assumption about how analysis units are constructed.
