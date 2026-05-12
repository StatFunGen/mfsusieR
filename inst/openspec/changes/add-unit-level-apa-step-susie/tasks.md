# Tasks

## 1. Step-basis module

- [ ] 1.1 Add `R/apa_step_basis.R` with `apa_step_basis()`,
      `apa_step_Xb()`, `apa_step_Xtr()`, `apa_step_colnorm()`,
      `apa_step_center_cache()`, and `apa_step_explicit()`.
- [ ] 1.2 Use manuscript notation in docs and comments:
      `i` = analysis unit, `j` = observed position/bin,
      `t` = candidate PAS / breakpoint.
- [ ] 1.3 Support upstream orientation
      `S[j, t] = 1(x_j <= c_t)` using prefix / reverse-prefix
      cumulative sums.
- [ ] 1.4 Support downstream orientation
      `S[j, t] = 1(x_j > c_t)` using suffix cumulative sums or an
      exactly equivalent transformation.
- [ ] 1.5 Validate ordered positions, candidate mapping, duplicate
      candidates, missing values, and candidate names.
- [ ] 1.6 Add tests comparing every operator to `apa_step_explicit()`
      for upstream/downstream, weighted/unweighted, and irregular
      candidate spacing.
- [ ] 1.7 Mirror the susieR trend-filtering helper pattern:
      `compute_tf_Xb()` -> `apa_step_Xb()`,
      `compute_tf_Xty()` -> `apa_step_Xtr()`, and
      `compute_tf_d()` / `compute_colstats()` ->
      `apa_step_colnorm()` plus any weighted centered denominators.
- [ ] 1.8 Add a full-candidate unweighted parity test showing that
      the APA upstream/downstream basis is equivalent to an explicit
      0th-order changepoint design up to the documented orientation
      convention.
- [ ] 1.9 Implement and test candidate coordinate mapping: integer
      candidates use checked indices; coordinate candidates map to the
      documented upstream/downstream position index; boundary
      candidates with insufficient information are flagged.
- [ ] 1.10 Cache weighted step means, centered denominators, and
      low-information flags per analysis unit and candidate without
      materializing dense `S`.

## 1A. Candidate PAS pre-scan module

- [ ] 1A.1 Add `R/apa_prescan.R` with `apa_prescan()` as a separate
      coverage-based candidate construction function. It may accept
      `initial_candidates`, but it must not implement external PAS
      discovery.
- [ ] 1A.2 Default to a dense scan with `R = J` and `d_r = x_r`;
      support a custom coarser grid as a speed option.
- [ ] 1A.3 Treat `Y` as the working coverage outcome; if `offset` is
      supplied, subtract it, and if `offset = NULL`, use zero
      position-specific offset.
- [ ] 1A.4 Always fit a weighted intercept in each marginal scan by
      weighted centering; do not require or compute an average offset.
- [ ] 1A.5 Do not estimate a flexible smooth GC, mappability, or
      5'/3' bias offset inside `apa_prescan()`; accept such correction
      only as preprocessing or as a supplied fixed offset.
- [ ] 1A.6 If `initial_candidates` is supplied, validate and map them
      with the same basis rules as the final fit, retain them in the
      output, fit an additive weighted-intercept step model on those
      candidates per unit, and scan the residuals. If
      `initial_candidates = NULL`, scan offset-adjusted working
      coverage directly.
- [ ] 1A.7 Treat initial-candidate residualization as candidate
      construction only: it must not return usage, must not define the
      final Bayesian posterior, and must not turn coverage scan scores
      into strong final prior weights.
- [ ] 1A.8 Compute `bhat_ir^scan`, `shat2_ir`, and
      `m_ir = log BF^+(bhat_ir^scan, shat2_ir, tau2_scan)`, and
      document that `shat2_ir` is a sampling variance, not a score.
- [ ] 1A.9 Skip or flag scan positions with too little effective
      weight on either side or with a near-zero centered denominator.
- [ ] 1A.10 Compute the pooled score curve by weighted sums across
      analysis units after normalizing positive `unit_weights` to
      mean one.
- [ ] 1A.11 Seed candidate regions from local maxima in the pooled
      score curve; by default retain finite local maxima with
      positive pooled score, and use `max_peaks` only as an optional
      high-score cap.
- [ ] 1A.12 Expand retained peaks by `region_half_width`, default
      exactly 100 nt, and retain fine-scale plus flanking
      local-background candidates inside expanded regions.
- [ ] 1A.13 Collapse nearby candidates within `merge_distance`,
      default exactly 25 nt, using a documented representative rule
      while preserving distinct peaks outside that distance; allow
      users to disable merging with `merge_distance = NULL` or
      `merge_distance <= 0`.
- [ ] 1A.14 Return candidates as the union of retained
      `initial_candidates`, scan-region positions, and local
      background positions, with source metadata and region
      diagnostics.
- [ ] 1A.15 Implement the dense scan with prefix sums or equivalent
      linear-time operations; do not fit `R` separate dense
      regressions per unit.
- [ ] 1A.16 Add tests comparing prefix-sum scan estimates to explicit
      weighted least-squares with an intercept, including supplied
      offsets, no-offset behavior, low-side-weight filtering, local
      background retention, multiple nearby peaks, merge behavior, and
      source-union metadata.
- [ ] 1A.17 Add a conditional-scan test where an initial candidate is
      retained, residualization removes its fitted contribution, and
      the residual scan adds a second true candidate.
- [ ] 1A.18 If `tau2_scan = NULL`, estimate it once from scan
      marginal statistics using the positive-excess rule and record
      scan-scale diagnostics.

## 2. Annotation prior module

- [ ] 2.1 Add `R/apa_prior_annotation.R` with
      `apa_prior_uniform()`, `apa_prior_from_scores()`, and
      `apa_prior_from_annotations()`.
- [ ] 2.2 Ensure every prior helper returns a strictly positive,
      finite vector summing to one.
- [ ] 2.3 Ensure annotation features are only used as location
      evidence and never encode effect signs.
- [ ] 2.4 Add `prior_strength` and `prior_floor` controls to score
      and annotation helpers; do not add a `temperature` alias in the
      first implementation.
- [ ] 2.5 Reject invalid prior floors, including negative values and
      values with `T * prior_floor >= 1`.
- [ ] 2.6 Normalize direct `prior_weights` after enforcing the same
      positive floor.
- [ ] 2.7 Store prior diagnostics, including prior odds summaries such
      as `max(pi) / min(pi)`.
- [ ] 2.8 Add tests for normalization, invalid inputs,
      monotone prior-strength behavior, positive floor behavior, and
      explicit fallback to uniform.
- [ ] 2.9 Keep the first annotation implementation simple: support
      direct `prior_weights`, direct numeric `scores`, or a supplied
      small linear score from feature columns; do not train a
      genome-wide prediction model or classifier.
- [ ] 2.10 If `apa_prior_from_annotations()` is called with
      `coefficients = NULL`, require a numeric column named `score`;
      otherwise error rather than inventing feature weights. If
      coefficients are supplied, require names matching feature
      columns or length equal to the number of columns.

## 3. Positive prior module

- [ ] 3.1 Add `R/apa_prior_positive.R` with
      `apa_halfnormal_lbf()`, `apa_halfnormal_moments()`, and
      `apa_halfnormal_ser()`.
- [ ] 3.2 Use `beta` in function documentation for effect magnitude;
      reserve `alpha` for posterior location probabilities.
- [ ] 3.3 Compute `log(2 * Phi(x))` stably for large negative `x`.
- [ ] 3.4 Return finite posterior mean and second moment under
      machine-small `shat2` and `tau2`.
- [ ] 3.5 Add numerical-integration tests for logBF and posterior
      moments.

## 4. Cross-unit module

- [ ] 4.1 Add `R/apa_cross_unit.R` with
      `cross_unit_prior_weighted()`.
- [ ] 4.2 Add S3 method
      `combine_outcome_lbfs.mf_prior_cross_unit_weighted()`.
- [ ] 4.3 Validate `unit_weights`: finite, non-negative, length
      matching the number of outcomes; allow zero weights.
- [ ] 4.4 Normalize positive `unit_weights` to mean one before use,
      leaving zeros as zeros; store raw and normalized weights in
      diagnostics.
- [ ] 4.5 Add tests against manual weighted logBF sums and tests
      showing that multiplying all positive weights by a constant does
      not change the combined logBF.

## 5. Initialization module

- [ ] 5.1 Add `R/apa_init.R` with `apa_marginal_init()`,
      `apa_l0learn_init()`, `apa_init_from_step_coef()`,
      `apa_estimate_tau2()`, and
      `apa_estimate_residual_variance()`.
- [ ] 5.2 Add wrapper controls `init_method = c("marginal",
      "none", "l0learn")` and `l0_control = list()`. The default
      must be `"marginal"`.
- [ ] 5.3 Implement `apa_marginal_init()` using matrix-free marginal
      step statistics from `apa_step_Xtr()` and `apa_step_colnorm()`;
      do not create a dense or sparse explicit step matrix.
- [ ] 5.3A Use the same weighted-centering convention as the main fit:
      subtract offsets, center `Y - offset` per unit, and compute
      marginal `bhat_it` and `shat2_it` using cached centered
      denominators.
- [ ] 5.4 Pool default marginal support across units using a
      documented positive-evidence score such as
      `score_t = sum_i unit_weights_norm[i] * logBF_it^+`, retain at
      most `min(L, max_nonzero, T)` candidates, and build a
      SuSiE-compatible one-hot initialization.
- [ ] 5.5 Under the default upstream-positive model, initialize
      unit-specific starting effects from positive marginal
      estimates and use zero for non-positive marginal evidence.
- [ ] 5.5A Add `apa_estimate_tau2()` for the first implementation's
      fixed slab scale. If the user supplies a positive scalar
      `tau2`, keep it unchanged and record `tau2_source = "user"`.
- [ ] 5.5B When `tau2 = NULL`, reuse marginal `bhat_it` and
      `shat2_it` from `apa_marginal_init()`, compute
      `e_it = max(bhat_it, 0)^2 - shat2_it`, and keep finite positive
      excess estimates.
- [ ] 5.5C Estimate the default `tau2` as a robust weighted quantile
      of retained positive excess estimates, defaulting to the
      weighted median with normalized positive-mean-one
      `unit_weights`, and enforce a positive `tau2_floor` on the
      working-outcome scale.
- [ ] 5.5D If no finite positive excess estimates are available, set
      `tau2 = tau2_floor`, record
      `tau2_source = "floor_no_positive_excess"`, and store `tau2`,
      `tau2_floor`, quantile, and retained-excess count in diagnostics.
- [ ] 5.5E Define `y_scale` for `tau2_floor` as the median across
      units of finite weighted variances of centered `Y - offset`,
      with a documented positive fallback.
- [ ] 5.5F Implement `apa_estimate_residual_variance()` as a separate
      APA helper. Use supplied `residual_variance` if present;
      otherwise use centered working coverage and a first-difference
      MAD estimate. For varying row precision weights, standardize
      adjacent differences by
      `sqrt(1 / w_ij + 1 / w_i,j-1)`. Use weighted-variance and
      positive-floor fallbacks, and return diagnostics per unit.
- [ ] 5.5G Keep initialized `sigma2_i` fixed by default. Only update
      residual variance when `estimate_residual_variance = TRUE`, and
      never update `tau2` in the residual-variance update path.
- [ ] 5.6 In `apa_l0learn_init()`, check that `L0Learn` is installed;
      otherwise return `NULL` with a diagnostic and let the wrapper
      use marginal initialization.
- [ ] 5.7 Guard optional explicit sparse step-matrix construction with
      `max_explicit_entries`; do not create `S` except for
      `init_method = "l0learn"`.
- [ ] 5.8 For each analysis unit `i`, when `init_method = "l0learn"`,
      allow the optional L0 path only when row weights are absent or
      constant within every unit; otherwise return `NULL` with a
      diagnostic and fall back to marginal initialization.
- [ ] 5.9 For the allowed L0 path, run
      `L0Learn::L0Learn.cvfit(S, ytilde_i, penalty = "L0", ...)`,
      rely on L0Learn's intercept, and remove that intercept after
      fitting.
- [ ] 5.10 Select the lambda rule, initially `cv_min` using
      `which.min(fit$cvMeans[[1]])`, matching the susieR vignette.
- [ ] 5.11 Extract coefficients, remove the intercept, and map nonzero
      coefficients to candidate PAS indices.
- [ ] 5.12 Under the default upstream-positive model, discard
      non-positive L0 coefficients; do not truncate them to positive
      values by default.
- [ ] 5.13 Combine L0-selected candidates across units by weighted support,
      `score_t = sum_i unit_weights_norm[i] * max(beta_hat_it, 0)^2`,
      and keep at most `min(L, max_nonzero, T)` candidates.
- [ ] 5.14 Create a SuSiE-compatible initialization with one-hot
      `alpha[l, t_l] = 1` and mfsusieR-shaped unit-specific starting
      effects: `mu[[l]][[i]][t_l, 1] = beta_init[i, l]` and
      `mu2[[l]][[i]][t_l, 1] = beta_init[i, l]^2`.
- [ ] 5.15 For selected L0 candidates missing in a unit, fill starting
      effects using the positive marginal WLS estimate at that
      candidate, or zero if no positive evidence exists.
- [ ] 5.16 Store initialization metadata on the fit: method used,
      selected candidates, scores, optional lambda rule, package
      availability, fallback reason, and explicit-matrix size when
      relevant.
- [ ] 5.17 Add tests for marginal initialization, no explicit-matrix
      use in the marginal path, positive filtering, `L` cap,
      default `tau2` estimation and floor fallback, optional
      no-`L0Learn` fallback, optional threshold fallback, and
      residual-variance initialization fallbacks. Include a test that
      non-constant row weights make L0 initialization fall back to
      marginal initialization.
- [ ] 5.18 Validate the warm-start object before fitting: finite
      `alpha`, `mu`, `mu2`, non-negative `V` if present, normalized
      alpha rows, dimensions matching `L`, `T`, and the number of
      analysis units, and no negative starting `beta` under the
      upstream-positive model.
- [ ] 5.19 Add direct unit tests for
      `apa_estimate_residual_variance()` covering unweighted,
      weighted, zero-MAD, non-finite, and user-supplied
      `residual_variance` cases.

## 6. APA data subclass and wrapper

- [ ] 6.1 Add `R/apa_mfsusie.R` with an APA data constructor
      `create_mf_apa_individual()`.
- [ ] 6.2 Add public wrapper `apa_susie()`, accepting
      `Y`, `pos`, `candidates`, `L`, `offset`, row precision
      `weights`, `prior_weights`, `annotations`, `unit_weights`,
      `tau2`, `orientation`, `prior_strength`, `prior_floor`,
      `init_method`, `l0_control`, `residual_variance`,
      `estimate_residual_variance = FALSE`, and standard SuSiE
      control arguments.
      Do not expose
      `estimate_prior_variance` in the first APA wrapper; `tau2`
      controls the slab scale.
- [ ] 6.3 Map the input to mfsusieR's multi-response shape:
      rows are positions, predictors are candidate PAS step
      functions, outcomes are analysis units, and each outcome has
      one response column.
- [ ] 6.4 Ensure the wrapper uses the matrix-free basis object in
      the default path and does not materialize a dense step matrix.
- [ ] 6.4A After candidate validation, reduce `L` to `min(L, T)` with
      diagnostics; if no candidate maps to an informative step or no
      unit-candidate pair is informative, stop with an informative
      error before IBSS.
- [ ] 6.5 Initialize scalar `tau2` when omitted using
      `apa_estimate_tau2()` from marginal positive step coefficients;
      do not update `tau2` inside IBSS.
- [ ] 6.6 Pass annotation-derived or user-supplied `prior_weights`
      into the existing SuSiE variable-selection prior path.
- [ ] 6.7 Pass `unit_weights` through the weighted cross-unit
      combiner after normalization to positive mean one.
- [ ] 6.8 Unless `init_method = "none"` or a user-supplied
      `model_init` is provided, call `apa_marginal_init()` before the
      IBSS fit and pass its result as `model_init`.
- [ ] 6.9 If `init_method = "l0learn"` is explicitly requested, call
      `apa_l0learn_init()` before the IBSS fit; if it returns `NULL`,
      fall back cleanly to `apa_marginal_init()` with diagnostics.
- [ ] 6.10 Treat `Y` as a working coverage outcome; accept `offset`,
      `weights`, and already-corrected coverage, but do not implement
      a built-in library-size, GC-bias, or 5'/3' bias correction
      pipeline in this change. Document that default inputs should
      preserve linear coverage contrasts; do not add rank-normalized
      or quantile-normalized coverage preprocessing inside the APA
      wrapper.
- [ ] 6.10B Initialize per-unit residual variance with
      `apa_estimate_residual_variance()` when omitted. Keep it fixed
      by default and update it only when
      `estimate_residual_variance = TRUE`.
- [ ] 6.10A Validate row precision weights. Accept `NULL`, a length-`J`
      vector recycled across units, a `J x n` matrix, or a list of
      `n` length-`J` vectors. Reject missing, negative, non-finite,
      or zero-total weights for any unit.
- [ ] 6.11 Do not fit an explicit intercept inside IBSS. Store
      weighted response centers and step centers so the centered fit
      is equivalent to a unit-specific baseline model.
- [ ] 6.12 Reconstruct fitted working coverage as
      `offset_i + ybar_i + sum_l sum_t alpha_lt * mu_lit * Sc_it`.

## 7. Required S3 method overrides

- [ ] 7.1 Implement `compute_residuals.mf_apa_individual()` using
      `apa_step_Xb()` and `apa_step_Xtr()`.
- [ ] 7.2 Implement `update_fitted_values.mf_apa_individual()` using
      `apa_step_Xb()`.
- [ ] 7.3 Implement `loglik.mf_apa_individual()` using
      `apa_halfnormal_lbf()` and the existing SuSiE softmax over
      candidate locations.
- [ ] 7.4 Implement `calculate_posterior_moments.mf_apa_individual()`
      using `apa_halfnormal_moments()`.
- [ ] 7.5 Implement `compute_ser_statistics.mf_apa_individual()`
      explicitly using cached `apa_step_Xtr()` residual inner
      products and `apa_step_colnorm()` denominators.
- [ ] 7.5A Compute SER statistics from centered partial residuals:
      `x_it = sum_j w_ij Sc_it(j) Rc_il(j)`,
      `bhat_it = x_it / d_it`, and `shat2_it = sigma2_i / d_it`,
      with low-information pairs returning uninformative logBFs.
- [ ] 7.6 Implement no-op
      `pre_loglik_prior_hook.mf_apa_individual()` and
      `post_loglik_prior_hook.mf_apa_individual()` so `tau2` is fixed
      during IBSS.
- [ ] 7.7 Implement `SER_posterior_e_loglik.mf_apa_individual()`
      without dense design multiplication.
- [ ] 7.8 Implement `update_variance_components.mf_apa_individual()`
      and `Eloglik.mf_apa_individual()` without calling helpers that
      assume dense wavelet coefficient arrays.
- [ ] 7.9 Implement `update_model_variance.mf_apa_individual()` for
      residual-variance updates and derived-cache refresh, mirroring
      `update_model_variance.mf_individual` but using APA residual
      structures. This method must not update `tau2`.
- [ ] 7.9A Update residual variance per unit as
      `sigma2_i = max(sigma2_floor_i, ERSS_i / sum_j w_ij)`, where
      `ERSS_i` is the weighted expected squared residual including
      posterior second-moment terms, not just squared posterior mean
      residuals.
- [ ] 7.10 Register all `mf_apa_individual` S3 methods using the
      same mechanism as the corresponding `mf_individual` method:
      `.onLoad` for susieR generics already registered there, and
      roxygen/NAMESPACE `S3method(...)` entries for hooks such as
      `pre_loglik_prior_hook` and `post_loglik_prior_hook`.

## 8. Phenotype module

- [ ] 8.1 Add `R/apa_phenotype.R` with `apa_phenotype()`.
- [ ] 8.2 Return candidate PIP, credible sets, unit-by-candidate
      effect means, usage, optional balance, optional expected
      length, and uncertainty summaries.
- [ ] 8.2A Compute candidate PIP with susieR semantics,
      `1 - prod_l(1 - alpha_lt)`, preferably by reusing
      `susieR::susie_get_pip()`.
- [ ] 8.2B Build default credible sets with
      `susieR::susie_get_cs(fit, coverage = coverage)` without
      passing dense `X`; skip correlation-purity filtering by
      default. Optional purity requires a bounded-size step-basis
      `Xcorr` path.
- [ ] 8.2C Add credible-set diagnostics: size, `cs_fraction = size/T`,
      claimed coverage, and `cs_is_diffuse` using default threshold
      `cs_fraction > 0.5`. Diffuse CSs are reported but not treated
      as confident breakpoint calls.
- [ ] 8.3 Normalize usage by unit and return missing usage with a
      diagnostic when total inferred abundance is below tolerance.
- [ ] 8.3B Compute usage for all candidates before any PIP or
      credible-set reporting filter. PIP and CS summaries are
      uncertainty diagnostics; they must not be prerequisites for
      constructing the usage matrix.
- [ ] 8.3A Define first-version usage as relative modeled
      step-derived candidate usage. Do not include the centered
      baseline/intercept component as a distal PAS or in the usage
      denominator.
- [ ] 8.4 Carry optional effective numerator-denominator summaries or
      precision weights when supplied, so downstream mixed generalized
      linear models can use the phenotype without redefining it.
- [ ] 8.5 Treat `cutpoints` as optional reporting boundaries, not
      model parameters; if `cutpoints = NULL`, return usage, PIPs,
      credible sets, and expected length without requiring balance.
- [ ] 8.6 For one or more cutpoints, validate coordinate scale, return
      one balance column per cutpoint, and record cutpoint values in
      metadata.
- [ ] 8.7 When `export_for_mixed_model = TRUE`, return uncertainty
      fields suitable for downstream mixed models where available,
      such as `usage_var` or `usage_se`, `logit_usage` and
      `logit_usage_se`, and optional `effective_success` /
      `effective_n`.
- [ ] 8.8 Do not implement QTL or ISSAC-style association testing in
      the phenotype module.
- [ ] 8.9 If uncertainty summaries are not implemented by a documented
      approximation, return them as `NA` with diagnostics rather than
      using an undocumented delta-method shortcut.
- [ ] 8.10 Add `usage_has_mass` and `usage_reportable` diagnostics.
      Low denominator still makes usage missing. Other rules based on
      maximum PIP, entropy, or CS diffuseness are reporting flags only
      and must not alter the returned usage matrix.
- [ ] 8.11 Add deterministic tests for usage normalization, zero/one/
      multiple cutpoints, expected length, missing denominator flags,
      mixed-model export fields, and dimensions.

## 9. Documentation and exports

- [ ] 9.1 Add roxygen docs for exported functions.
- [ ] 9.2 Update `NAMESPACE` via roxygen, including the new
      `combine_outcome_lbfs` S3 method.
- [ ] 9.3 Add a vignette or article stub describing the unit-level
      APA model and its modules.
- [ ] 9.4 Explain the 0th-order trend-filtering / Haar-like
      interpretation and why this is a dedicated step transform rather
      than the existing DWT path.
- [ ] 9.5 Ensure docs do not name the model after a specific
      unit-construction strategy.
- [ ] 9.6 Explain that `Y_i(j)` is a working coverage outcome and
      that bias correction is an upstream preprocessing assumption,
      not a novel APA-model component.
- [ ] 9.7 Explain how strong annotation is represented through prior
      odds using `prior_strength` and `prior_floor`, without hard
      filtering or forcing per-unit usage.
- [ ] 9.8 Explain cutpoint `r` as a reporting choice for optional
      proximal-distal summaries and document alternatives when no
      single contrast is well motivated.
- [ ] 9.9 Document uncertainty exports for downstream mixed models
      while making clear that association testing is out of scope.
- [ ] 9.10 Document `apa_prescan()` as optional input-variable
      screening, including `initial_candidates`, conditional
      residualization before scan, the offset/intercept rule, the
      definition of `m_ir` and `shat2_ir`, the default
      `region_half_width`, the merge distance, local-background
      candidates, and why coverage-derived scan scores are normally
      not used as strong final priors.
- [ ] 9.11 Add package documentation and code comments under the
      ASCII-only rule. User-facing documentation should be concise;
      implementation comments around nontrivial math must spell out
      formulas in LaTeX-style ASCII text, including `Y_i(j)`,
      `S[j,t]`, `sigma2_i`, `tau2`, `bhat_it = x_it / d_it`, and
      `shat2_it = sigma2_i / d_it`.
- [ ] 9.12 Code comments must explicitly distinguish row precision
      weights `w_ij`, cross-unit weights `unit_weights_norm[i]`,
      residual variance `sigma2_i`, and slab scale `tau2`.

## 10. Integration tests

- [ ] 10.1 Simulate one gene with two shared PAS locations and
      unit-specific positive effects; verify `apa_susie()` recovers
      high PIP near the true candidates.
- [ ] 10.2 Verify annotation prior upweights the intended candidate
      through `prior_weights` without filtering other candidates.
- [ ] 10.3 Verify unit weights change joint evidence in the expected
      direction on a two-unit toy example.
- [ ] 10.3A Verify normalized unit weights are invariant to global
      rescaling of all positive weights.
- [ ] 10.4 Verify default fitting path does not call
      `apa_step_explicit()`.
- [ ] 10.5 Verify marginal initialization selects plausible
      positive candidates on a simple simulated breakpoint example
      without calling `apa_step_explicit()`.
- [ ] 10.6 Verify `init_method = "l0learn"` is optional, preserves
      final output shape when available, and falls back to marginal
      initialization when unavailable or above the explicit-matrix
      threshold.
- [ ] 10.7 Verify increasing `prior_strength` increases prior odds for
      higher-scored candidates while `prior_floor` preserves positive
      probability for all candidates.
- [ ] 10.8 Verify `apa_phenotype()` exports mixed-model compatible
      uncertainty summaries without running any association model.
- [ ] 10.9 Verify `apa_prescan()` output can be passed to
      `apa_susie()` as the candidate set on a simulated multi-PAS
      gene.
- [ ] 10.10 Simulate a no-signal gene and verify the model does not
      report a confident APA phenotype without evidence.
- [ ] 10.11 Simulate duplicated or overlapping analysis units and
      verify downweighting duplicates reduces overconfident PIPs.
- [ ] 10.12 Simulate broad positional bias with a supplied offset and
      verify false scan peaks are reduced while a true sharp
      breakpoint remains detectable.
- [ ] 10.13 Verify the centered no-explicit-intercept implementation
      matches an explicit weighted-intercept WLS calculation on small
      fixed examples.
- [ ] 10.14 Verify default APA credible sets match
      `susieR::susie_get_cs()` called without `X` or `Xcorr`, and
      candidate PIPs match `susieR::susie_get_pip()`.
- [ ] 10.15 Verify row precision weights accept all documented shapes
      and reject missing, negative, non-finite, or zero-total weights
      for any analysis unit.
- [ ] 10.16 Verify `L > T` is reduced to `L = T` with diagnostics and
      no informative candidate stops before IBSS with an informative
      error.
- [ ] 10.17 Verify residual-variance updates match an explicit small
      weighted ERSS calculation that includes posterior second-moment
      terms.
- [ ] 10.18 Verify fitted coverage reconstruction restores offset and
      stored baseline while keeping the baseline outside the IBSS
      variable set.
- [ ] 10.19 Verify first-version usage excludes the baseline from the
      usage denominator and returns missing usage when the modeled
      step-derived denominator is below tolerance.
- [ ] 10.20 Verify diffuse credible sets and low maximum PIP can make
      `usage_reportable = FALSE` without changing non-missing usage
      values when the usage denominator is adequate.
- [ ] 10.21 Verify `init_method = "l0learn"` falls back to marginal
      initialization when any unit has non-constant row precision
      weights.
- [ ] 10.22 Simulate endpoint-only units: one unit uses only a
      proximal endpoint and another uses only a distal endpoint under
      the same candidate basis. Verify unit-specific usage
      concentrates on the corresponding candidate.
- [ ] 10.23 Verify conditional coverage scan: supplied
      `initial_candidates` are retained, the residual scan is run
      after conditioning on them, and additional true candidates can
      be added.

## 11. Performance tests

- [ ] 11.1 Add a small benchmark proving matrix-free operators match
      explicit operator results on a small example.
- [ ] 11.2 Add a skipped larger performance check documenting the
      expected `O(n * (J + T))` per-effect scaling.
- [ ] 11.3 Add an assertion or test hook ensuring no dense `J x T`
      design is created in the default fit path.
- [ ] 11.4 Add a pre-scan performance smoke test showing dense
      `R = J` scanning uses prefix sums and scales like `O(n * J)`.
- [ ] 11.5 Record fit dimensions, IBSS iteration count, initializer
      used, explicit-matrix use, and fallback diagnostics so users can
      evaluate empirical scaling on their data. Do not impose fixed
      hard size limits in the default matrix-free fitting path.

## 11A. Code-generation and output-quality checks

- [ ] 11A.1 Add test helpers that can monkey-patch
      `apa_step_explicit()` and fail if the default fit path calls it.
- [ ] 11A.2 Add a fit-object validator used in tests to check finite
      `alpha`, `mu`, `mu2`, `lbf`, `pip`, normalized priors, and
      recorded diagnostics.
- [ ] 11A.3 Add an `apa_phenotype()` output validator checking
      required fields, dimensions, usage row sums for reportable
      units, missing flags for low-information units,
      credible-set diffuseness flags, `usage_reportable`, and finite
      uncertainty summaries where implemented.
- [ ] 11A.4 Ensure optional dependencies (`L0Learn`) are guarded by
      `requireNamespace()` and tested through skip/fallback paths;
      the default marginal path must work without optional
      dependencies.
- [ ] 11A.5 Keep APA code in APA-specific files except for minimal
      S3 registration, `NAMESPACE`, `DESCRIPTION`, and documentation
      edits.

## 12. Validation commands

- [ ] 12.1 From `/Users/gw2411/Documents/GIT/mfsusieR/inst`, run
      `openspec validate add-unit-level-apa-step-susie --type change --strict --no-interactive`.
- [ ] 12.2 Run `devtools::test()`.
- [ ] 12.3 Run `devtools::document()`.
- [ ] 12.4 Run `R CMD check` or the repository's standard check
      command.
- [ ] 12.5 Update this OpenSpec change and any package documentation
      if API names change during coding.
