# Refactor exceptions ledger

This file records every line of the two port sources that mfsusieR
intentionally does NOT carry over. Doctrine lives in
`inst/openspec/changes/add-mfsusier-s3-architecture/design.md`
D11e. Reviewer pass on every Phase 3 numerical PR confirms the
ledger is up to date; PRs that omit lines without a corresponding
entry SHALL be blocked.

## Port sources in scope

- `mvf.susie.alpha/R/multfsusie.R` and supporting
  `operation_on_multfsusie_*.R`, `EM.R`, `ELBO_mutlfsusie.R`,
  `computational_routine.R`, `utils_wavelet_transform.R`, and the
  parts of `utils.R` and `utils_formatting.R` that the IBSS path
  touches.
- `fsusieR/R/susiF.R` and `susiF_workhorse.R` and supporting
  `operation_on_susiF_obj.R`, `wavelet_utils.R`, the parts of
  `operation_on_prior.R` that init the scale-mixture prior, the
  per-method smoothers in `operation_on_susiF_obj.R`
  (`TI_regression`, `smash_regression`, `HMM_regression`), and
  the parts of `utils.R` needed by the susiF path.

## Out of scope (file-level entries; no per-line walk required)

- fsusieR/R/EBmvFR.R
  Behavior: EB multivariate functional regression, a different
    model with no SuSiE structure.
  Decision: out-of-scope-EBmvFR
  Reason: design.md D5 and CLAUDE.md hard rule #2 record EBmvFR as
    an entirely separate algorithm. Not in v1; possibly addressed
    in a follow-up OpenSpec change after the SuSiE-track ports
    stabilize.
- fsusieR/R/EBmvFR_workhorse.R
  Behavior: Inner loop for EBmvFR.
  Decision: out-of-scope-EBmvFR
  Reason: As above.
- fsusieR/R/operation_on_EBmvFR_obj.R
  Behavior: Operations on the EBmvFR fit object.
  Decision: out-of-scope-EBmvFR
  Reason: As above.

## Entry schema for in-scope omissions

```
- <port_source>/R/<file>.R:<lo>-<hi>
  Behavior: <one-line description of what the original lines do>
  Decision: omit | replaced-by-<upstream> | deferred-to-<phase> |
            out-of-scope-EBmvFR
  Reason: <one-paragraph justification, citing OpenSpec change or
          manuscript section>
```

## In-scope omissions

### PR group 2 (R/utils_wavelet.R, 2026-04-25)

- fsusieR/R/utils.R:113-119
  Behavior: `for (col in zero_var_cols) { csd[col] <- 1 }`. Loop
    over zero-variance columns to replace `csd[col]` with 1.
  Decision: replaced-by-vectorized (col_scale uses
    `csd[csd == 0] <- 1`).
  Reason: Trivially equivalent to the vectorized form; no
    behavioural change. Cosmetic style cleanup logged here per
    D11e.

- fsusieR/R/utils.R:135-137
  Behavior: `d = n * cm^2 + (n - 1) * csd^2; d = (d - n * cm^2) /
    csd^2`. Two-step computation that algebraically simplifies to
    a constant `(n - 1)` per column since `csd[csd == 0] <- 1`
    was applied upstream.
  Decision: replaced-by-simplified-math (port computes
    `attr(x, "d") <- rep(nrow(x) - 1, ncol(x))` directly).
  Reason: The two-step form is wasteful compute and introduces
    ULP-level rounding noise via the `(a + b) - a` pattern. The
    simplification is correct math (algebraically equivalent),
    and the C2 fidelity test on the `d` attribute is relaxed to
    `tolerance = 1e-12` (well within contract C2's `<= 1e-8`
    floor) to accept the resulting machine-epsilon difference.
    The test file documents the fsusieR-bug nature inline. Net
    effect: cleaner mfsusieR code, identical observable behaviour
    on the `d` attribute up to machine epsilon.

- fsusieR/R/wavelet_utils.R:295-297
  Behavior: `verbose` message text uses raw newlines:
    `"Response matrix dimensions not equal to nx 2^J \n or
    unevenly spaced data \n interpolation procedure used"`.
  Decision: replaced-by-cleaner-message-text.
  Reason: The message text is user-facing diagnostic output, not
    part of any numerical contract. mfsusieR replaces with a
    semicolon-separated single line for readability. No effect on
    C2 / C3 contracts.

- fsusieR/R/wavelet_utils.R:148
  Behavior: `for (s in 1:(lev_res-1))` in `gen_wavelet_indx`.
    For `lev_res = 1` this expands to `1:0 = c(1, 0)`, executing
    the loop twice with `s = 1` and `s = 0`, producing a
    nonsensical result.
  Decision: replaced-by-seq_len (port uses `seq_len(lev_res - 1)`
    which yields `integer(0)` for `lev_res = 1`).
  Reason: `lev_res = 1` corresponds to a 2-sample input which is
    too small for a meaningful wavelet decomposition and is not
    part of the C2 contract floor (`T_1 in c(64, 128)` per
    design.md D11b). The fsusieR pattern is a real bug at this
    edge case; the port adopts the safe `seq_len` idiom. The C2
    test range (lev_res 6-7) is bit-identical between fsusieR and
    mfsusieR; the bug-fixed edge case is below the contract.

- fsusieR/R/wavelet_utils.R:18-37 + R/wavelet_utils.R:284-297
  Behavior: `interpol_mat` returns `grid` rescaled to a `[0,
    length(grid)]` integer-like scale, and `remap_data` then
    rescales it back to position units via
    `outing_grid <- start_pos + (end_pos - start_pos) /
    length(grid) * grid`.
  Decision: replaced-by-single-step (port folds the rescaling
    into `interpol_mat`, returning `outing_grid` in position
    units directly).
  Reason: The two-step composition is mathematically identical to
    the one-step linear interpolation, but the order of floating-
    point operations differs, producing a ULP-level difference
    (~3e-15) in `outing_grid`. Y output is bit-identical (max
    diff = 0). The reordering is well within contract C2's `<=
    1e-8` floor. The simplification eliminates redundant
    `min(pos)`, `max(pos)`, `min(raw_grid)`, `max(raw_grid)`
    passes and is locked by a regression test in
    `tests/testthat/test_utils_wavelet.R`.

- fsusieR/R/wavelet_utils.R:102-105 (DWT2's `min.scale = 10`)
  Behavior: `wavethresh::wd(..., min.scale = 10)` hardcoded inside
    the DWT helper.
  Decision: replaced-by-parameter-unification (port adds
    `max_scale = 10` argument to `dwt_matrix`, forwarded as
    wavethresh's `min.scale`).
  Reason: The `10` is the same conceptual quantity as `max_scale`
    elsewhere in `R/utils_wavelet.R`. Hardcoding it inside the
    DWT helper is a magic-number style issue the audit flagged.
    The port lifts it to a function-level parameter named
    consistently with the rest of the file. No behavioural change
    when `max_scale = 10` (the default); behaviour preserved.

### PR group 4 (R/prior_scale_mixture.R + R/prior_cross_modality.R, 2026-04-25)

- fsusieR/R/computational_functions.R:69-155 (cal_Bhat_Shat,
  full-Y case)
  Behavior: per-position marginal OLS regression (Bhat, Shat)
    used to seed the ash mixture grid in init_prior.default.
  Decision: replaced-by-`susieR::compute_marginal_bhat_shat`
    (new upstream helper landed in this session).
  Reason: per the upstream-first principle in
    `inst/notes/refactor-discipline.md` section 0, the per-
    position OLS regression belongs in susieR where multiple
    downstream packages can share it. mfsusieR calls
    `susieR::compute_marginal_bhat_shat(X, Y_m)` from
    `init_scale_mixture_prior_default`. Bit-identical to fsusieR
    output (max diff = 0 verified in
    `tests/testthat/test_prior_scale_mixture.R` C2 fidelity test).

- fsusieR/R/operation_on_prior.R:42-185 (init_prior.default,
  mixture_normal and mixture_normal_per_scale branches)
  Behavior: Per-modality data-driven prior init via ash on a
    sample of (Bhat, Shat); replicate per scale. Hardcodes the
    init pi vector to `c(0.8, 0.2/(K-1), ...)` regardless of
    user-side prior weight parameters.
  Decision: replaced-by-`R/prior_scale_mixture.R::init_scale_mixture_prior_default`,
    with the init pi parameterised on `null_prior_weight`
    (default 2) so the same parameter drives both data-driven
    and user-supplied-grid paths.
  Reason: behaviour-preserving port for the sd-grid (same
    `set.seed(1)` sequence, sample-size caps 5000 / 50000, ash
    invocation). The hardcoded `pi_null = 0.8` upstream constant
    is replaced by `pi_null = null_prior_weight / (K + 1)`,
    matching the formula in `distribute_mixture_weights`.
    Default `null_prior_weight = 2` gives `pi_null ≈ 2 / (K+1)`,
    a weaker sparsity init than upstream's 0.8 but consistent
    with the M-step null:data balance set by
    `mixsqp_null_penalty = 0.1`. The C2 fidelity test recovers
    upstream bit-identity by passing
    `null_prior_weight = 0.8 * (K + 1)` (test hack only;
    production uses the parameterised default).

- fsusieR/R/operation_on_prior.R:79-165 (lowc_wc filtering branch)
  Behavior: when `lowc_wc` (low-count wavelet coefficient
    indices) is non-NULL, drop those columns from the Bhat/Shat
    sample before the ash fit.
  Decision: deferred-to-Phase-3-future-PR (not implemented in PR
    group 4).
  Reason: the ledger comment block says low-count filtering is a
    smashr / signal-quality concern downstream. mfsusieR's v1
    init does not filter low-count coefficients. If the FDR
    investigation in Phase 5 surfaces a need, a follow-up PR can
    add the option. The C2 fidelity tests run with `lowc_wc =
    NULL` on both sides, so the deferral does not affect the
    contract.

- mvf.susie.alpha/R/operation_on_multfsusie_prior.R:17-138
  (Prior object construction for `multfsusie`)
  Behavior: per-(modality, scale) prior init in the multi-
    functional case.
  Decision: replaced-by-`R/prior_scale_mixture.R::mf_prior_scale_mixture`
    (per-modality wrapper).
  Reason: mfsusieR's wrapper iterates over modalities and calls
    the per-modality init. The original mvf.susie.alpha file is
    a thicker wrapper that includes EM update logic that lives
    in PR group 6 in mfsusieR. Per-modality prior init line
    range mapped over the new file structure.

### PR group 2 (R/dwt.R + R/data_class.R, 2026-04-25)

- mvf.susie.alpha/R/utils_wavelet_transform.R:24-50 (DWT2)
  Behavior: Row-wise wavelet transform of a matrix of curves.
    Identical in shape to `fsusieR::DWT2` but without the
    `min.scale = 10` argument; the wavethresh default applies
    (decompose to scale 0).
  Decision: replaced-by-`R/utils_wavelet.R::dwt_matrix` (already
    landed in PR group 2 task 2.2).
  Reason: A single canonical DWT helper is enough; mvf.susie.alpha
    and fsusieR each shipped a near-duplicate. The C3 fidelity
    test in `tests/testthat/test_dwt.R` ("matches mvf.susie.alpha
    DWT2 row-wise") verifies bit-identity on the same input
    (after the same `col_scale` step).

- mvf.susie.alpha/R/utils_wavelet_transform.R:64-68 (pack_dwt)
  Behavior: One-line helper computing `cbind(W$D, W$C)` where W
    is `DWT2(Y)`.
  Decision: replaced-by-inline (the `cbind(D, C)` pack step is
    inlined into `R/dwt.R::mf_dwt`).
  Reason: The packing convention is exposed in `mf_dwt` directly;
    a separate one-line helper is unnecessary.

- mvf.susie.alpha/R/multfsusie.R:255-318 (per-modality DWT pipeline
  inside `multfsusie()`)
  Behavior: Inline loop in the `multfsusie()` body that calls
    `fsusieR::remap_data` -> `fsusieR::colScale` -> `DWT2` ->
    `cbind(D, C)` -> `fsusieR::gen_wavelet_indx` per modality.
  Decision: replaced-by-`R/dwt.R::mf_dwt` (per-modality wrapper)
    invoked from `R/data_class.R::create_mf_individual` (per-fit
    constructor).
  Reason: The original buries the DWT pipeline inside the public
    fitting function; mfsusieR factors it into a reusable
    per-modality helper (`mf_dwt`) and a separate constructor
    (`create_mf_individual`). The pipeline order and arguments
    are bit-identical to mvf.susie.alpha's; the C3 DWT fidelity
    test verifies this on the matrix output.

- fsusieR/R/operation_on_susiF_obj.R:1785-1796 (inverse-DWT
  reconstruction skeleton)
  Behavior: Constructs a zero `wavethresh::wd` object, sets `$D`
    and the last entry of `$C`, applies `wavethresh::wr` to
    reconstruct the position-space curve. Used by `update_cal_fit_func`
    and downstream coef / predict.
  Decision: replaced-by-`R/dwt.R::mf_invert_dwt` (PR group 2 task
    2.4).
  Reason: Same skeleton, factored into a reusable inverse-DWT
    helper. The 2e-9 wd/wr roundtrip noise inherent to
    `filter.number = 10` is documented in
    `tests/testthat/test_dwt.R` ("mf_dwt + mf_invert_dwt is
    identity to within wd/wr precision"), with tolerance set to
    `1e-8` (the contract C2 floor).

### PR group 6 (R/individual_data_methods.R variance + ELBO, 2026-04-25)

- mvf.susie.alpha/R/operation_on_multfsusie_obj.R:1046-1096
  (`get_ER2.multfsusie`)
  Behavior: ER2 per modality is computed as
    `sum((Y - X*postF)^2) - sum(postF^2) + sum(postF2_sd2)`,
    where `postF = sum_l alpha_l * mu_l` and
    `postF2_sd2 = sum_l alpha_l * mu2_l`. The bias correction
    omits the predictor-weights factor and aggregates effects
    via `sum(postF^2)` (sum-then-square) instead of
    `sum_l ||X * (alpha_l * mu_l)||^2` (per-effect square-then-sum).
  Decision: replaced-by-susieR-formula in
    `R/individual_data_methods.R::mf_get_ER2_per_position`.
  Reason: susieR's `get_ER2.individual` uses the correct
    `E_q[||Y - X*beta||^2]` decomposition for the SER posterior,
    `||Y - X*postF||^2 + sum_l ((alpha_l . pw)^T mu2_l - ||X (alpha_l . mu_l)||^2)`,
    derived from `Cov(beta_lj, beta_lk) = -alpha_lj alpha_lk mu_lj mu_lk`
    for j != k under the single-effect mixture posterior. The port
    source's formula deviates in two ways: (1) missing the
    `predictor_weights` factor (~ O(n) for scaled X), (2)
    sum-then-square vs square-then-sum across effects. Magnitude:
    sigma^2 underestimated by up to ~n on scaled X, which inflates
    Bayes factors and PIPs. Plausible mechanism for the FDR
    miscalibration observed in `inst/notes/paradigms/mvf-original.md`;
    Phase 5 will verify with controlled simulation. This is a
    Pattern A port-source-bug fix per refactor-discipline section 3.
    The C3 test (PR 7e) shall assert the deviation explicitly with
    a reference to this entry.

- mvf.susie.alpha/R/computational_routine.R:395-431
  (`estimate_residual_variance.multfsusie`)
  Behavior: divides `get_ER2` output by `n * T_m` (legacy
    `shared_per_modality` mode) or NA for per-(scale, modality).
  Decision: replaced-by-`R/individual_data_methods.R::update_variance_components.mf_individual`.
  Reason: the divisor is correct; the bug is upstream in `get_ER2`
    (separate entry above). Our `update_variance_components`
    branches on `params$residual_variance_method`:
    `shared_per_modality` divides `sum(ER2)` by `n * T_padded[m]`
    (matching the port-source divisor with the corrected ER2);
    `per_scale_modality` (default) divides `sum(ER2[idx_s])` by
    `n * |idx_s|` per scale.

- fsusieR/R/computational_functions.R `fit_hmm` (HMM-mu-subset)
  Behavior: after the pre-EM forward-backward, `idx_comp` selects
    a subset of the K_full mixture states whose average posterior
    occupancy clears `thresh`. Upstream subsets `prob[, idx_comp]`
    and `P[idx_comp, idx_comp]`, but does NOT subset `mu`. Inside
    the EM loop (`for (k in 2:K)`), the line `mu_ash <- mu[k]`
    therefore reads the k-th value of the un-subset grid rather
    than `mu[idx_comp[k]]`, breaking the state-emission alignment
    required by the forward-backward equations whenever
    `idx_comp` is not the contiguous prefix `1..K_full`.
  Decision: fixed-in `R/post_smooth_hmm.R::mf_fit_hmm` by
    appending `mu <- mu[idx_comp]` to the subset block, alongside
    `prob[, idx_comp]` and `P[idx_comp, idx_comp]`. The returned
    `mu` therefore has length K (= length(idx_comp)) rather than
    K_full.
  Reason: the upstream author identified this as a missing fix
    in commit 9f89333 (2026-03-12 14:30) with explicit comment
    `# *** FIX: subset mu to match the reduced state space ***`
    and `# was missing in previous version`. Four hours later in
    commit fc806a5 (2026-03-12 18:46, "Bump up version") the line
    was deleted during a Baum-Welch restructure that preserved
    the sibling subsettings on `prob` and `P`. There is no
    algorithmic justification for the deletion (Rabiner 1989
    §III.B, Cappé-Moulines-Rydén 2005 §6: state reindexing in an
    HMM is consistent only if all state-indexed quantities are
    reindexed together).
  Resolution (2026-04-26): fsusieR pull request
    stephenslab/fsusieR#31 restored the `mu <- mu[idx_comp]`
    line. Both the contiguous-prefix and non-contiguous
    regimes now bit-match upstream at `tolerance = 0`; the
    Pattern A deviation-asserting test in
    `tests/testthat/test_post_smooth_HMM.R` is replaced by
    a uniform bit-identity assertion across both regimes.


### Cross-package audit fix-now bundle (2026-04-26)

Three ledger entries from the 2026-04-26 audit. Code edits and
math derivations are colocated in
`inst/notes/cross-package-audit-derivations.md` and
`inst/notes/cross-package-audit-summary.md`.

- mvf.susie.alpha/R/operation_on_multfsusie_obj.R credible-band
  formula (position-space SD via `abs(invert_dwt(sqrt(var_w)))`)
  Behavior: the upstream credible band derives a position-space
    standard deviation by inverse-DWT-ing the per-coefficient
    standard deviations (square roots of the wavelet-domain
    variances) and taking the absolute value. This treats a
    linear combination of standard deviations as if it were the
    standard deviation of the corresponding linear combination
    of independent variables, which is wrong: the correct
    identity is `Var(pos[t]) = sum_k W^T_{t,k}^2 * Var(w[k])`,
    not `SD(pos[t]) = |sum_k W^T_{t,k} * SD(w[k])|`.
  Decision: replaced-by-`R/utils_wavelet.R::mf_invert_variance_curve`
    plus the call site in `R/mfsusie_methods.R::.post_smooth_scalewise`.
    `mf_invert_variance_curve(var_w, T_basis, ...)` builds the
    column-by-column squared inverse-DWT matrix `W_sq` and returns
    `W_sq %*% var_w`; the position-space SD is then
    `sqrt(mf_invert_variance_curve(var_w, ...))`. Pattern A
    port-source-bug fix per `inst/notes/refactor-discipline.md`
    section 3.
  Reason: the closed-form variational posterior on a wavelet
    basis is Gaussian per coefficient; mapping the variance
    through `W` requires squaring each row of W, not taking
    `|W * sqrt(var)|`. The upstream formula systematically
    miscalibrates credible bands at positions where multiple
    basis functions have non-trivial support (i.e., everywhere
    except the trivial T_basis = 1 case). Audit ID B-1.1.

- mvf.susie.alpha/R/EM.R:58-65 (nullweight scaled by `K_f + K_u`)
  Behavior: upstream multiplies the user-facing `nullweight`
    argument by `K_f + K_u` (number of functional + scalar
    outcomes) before passing it to the per-outcome mixsqp M
    step. mfsusieR's per-(outcome, scale) M step lives in
    `R/individual_data_methods.R::optimize_prior_variance.mf_individual`
    and originally consumed `params$mixture_null_weight %||% 0.05`
    unscaled. (Public-API rename from `mixsqp_null_penalty` to
    `mixture_null_weight`; default lowered from 0.7 to 0.05 as
    part of the per-(outcome, scale) calibration done in PR group
    6 because the per-(m, s) M step splits the regularization
    pressure across `M * S` independent solves rather than mvf's
    single joint solve.)
  Decision: scale-mfsusieR-by-M in
    `R/individual_data_methods.R::optimize_prior_variance.mf_individual`:
    `mixture_null_weight <- (params$mixture_null_weight %||% 0.05) *
                             max(1L, data$M)`.
  Reason: derivation in
    `inst/notes/cross-package-audit-derivations.md` section 1.
    With M outcomes, the joint per-effect SNP posterior alpha_l
    concentrates by an O(M) factor in the strong-signal regime
    (variance of joint lbf scales linearly with M, and the
    softmax over j amplifies the relative weight of top SNPs by
    a factor whose log scales linearly with M). The mixsqp
    first-order condition is `LHS = lambda / (1 - pi_1)` with
    LHS dominated by `alpha_top * G(pi_1, BF_top)`. To keep
    the equilibrium pi_1 invariant in M, lambda must scale
    linearly in M. mfsusieR's per-outcome M step uses the
    same alpha as mvf.susie.alpha's joint M step (the joint
    posterior is invariant to whether the M step is per-outcome
    or joint), so the same M-fold scaling applies. M = 1 fits
    are unchanged (`max(1L, 1L) = 1`); M >= 2 fits now apply
    M-fold stronger regularization on the mixture null. Audit
    ID C-1.3.

- mvf.susie.alpha/R/ELBO_mutlfsusie.R:269-282 (entropy term `o`
  added to the upstream ELBO)
  Behavior: upstream adds `o = sum_l sum_j alpha_lj * log(p /
    pmax(alpha_lj, 1e-6))` to its ELBO formula. The upstream
    per-effect KL is computed as
    `KL_l^{mvf} = -loglik_SFR(l) - loglik_SFR_post(l)`
    where `loglik_SFR_post(l) = E_q[log p(Y | b_l)]` is a
    typically large negative number; this differs in sign from
    the susieR convention by `-2 * E_q[log p(Y | b_l)]`. The
    `+o` term is bounded by `L * log p` and does not offset the
    sign error, so the upstream ELBO is mathematically incorrect.
  Decision: omit-the-o-term-in-mfsusieR. mfsusieR's
    `R/ibss_methods.R::get_objective.mf_individual` returns
    `Eloglik - sum(KL)` per the susieR convention. No code change
    is needed; this entry documents the deliberate omission.
  Reason: derivation in
    `inst/notes/cross-package-audit-derivations.md` section 2.
    The susieR per-effect KL formula
    `KL_l = -lbf_l - L_null + E_q[log p(Y | b_l)]`
    is the standard variational decomposition of
    `KL(q(b_l) || p(b_l))` into the categorical KL plus the
    weighted Gaussian KL. `Eloglik - sum_l KL_l` is therefore
    the complete ELBO. The watertight check is monotone non-
    decreasing ELBO across IBSS iterations at machine precision
    (`tests/testthat/test_variance_and_elbo.R`, "ELBO is non-
    decreasing post-iter-1 at machine precision"). Audit ID C-4.4.

### PR group ?? (R/adjust_covariates.R + R/em_helpers.R, 2026-04-26)

- fsusieR/R/EBmvFR.R, fsusieR/R/EBmvFR_workhorse.R,
  fsusieR/R/operation_on_EBmvFR_obj.R (covariate-adjustment
  pathway only; the `adjust = TRUE` exit branch).
  Behavior: wavelet-domain empirical-Bayes regression of a
    functional response Y on a covariate matrix Z, returning
    `Y_adjusted = Y - Z %*% fitted_func`.
  Decision: in-scope-as-of-round-2 for the covariate-adjustment
    use case. Assembled from existing primitives
    (`init_scale_mixture_prior_default`, `compute_marginal_bhat_shat`,
    `ashr::postmean/postsd`, `mixsqp::mixsqp`, `dwt_matrix`,
    `mf_invert_dwt`) plus a new `R/adjust_covariates.R` module.
    No port of the EBmvFR.R / EBmvFR_workhorse.R files
    themselves; the algorithm is reassembled. The full EB-
    multivariate-functional-regression model as a SuSiE
    alternative remains out of scope (the original
    out-of-scope-EBmvFR ruling stands for that purpose).
  Reason: round 2 vignette work needs a public covariate-
    adjustment utility (`mf_adjust_for_covariates`). The
    upstream `EBmvFR(adjust = TRUE)` exit pathway is a
    well-defined subset of the EBmvFR machinery and assembling
    it from existing mfsusieR primitives stays within the
    scope-discipline rule for the SuSiE-track ports.

- fsusieR/R/operation_on_EBmvFR_obj.R:351 (update_effect.EBmvFR
  storing `MLE_wc$Shat` into `obj$MLE_wc2`)
  Behavior: the slot named `MLE_wc2` (suggesting "second
    moment", i.e., squared standard error) is upstream
    populated with `MLE_wc$Shat` (the standard error itself,
    not the second moment). The downstream `update_prior_EBmvFR`
    then calls `L_mixsq(..., Shat = sqrt(obj$MLE_wc2[[1]]))`,
    which evaluates to `sqrt(Shat)`, NOT `Shat`. The
    "Shat" argument fed to L_mixsq is therefore in units of
    `sqrt(Shat)`, off by a square root.
  Decision: fixed-in `R/adjust_covariates.R::mf_residualize_wavelet_eb`.
    `MLE_wc2[j, ] <- shat_j^2` stores the actual second moment.
    The L_mixsq call passes `sqrt(MLE_wc2[, idx_s, drop = FALSE]) = Shat`,
    which is the correct standard error. Pattern A bug fix per
    refactor-discipline section 3.
  Reason: this is a slot-name-vs-content bug, not a defensible
    algorithmic choice. The author labels the slot `_wc2`
    explicitly to indicate a second moment; the storage line
    contradicts the label and the downstream consumer treats
    the value as if it were the second moment. The mathematical
    consequence is wrong-units arguments to mixsqp.

- fsusieR/R/operation_on_EBmvFR_obj.R:124-132 (get_ER2.EBmvFR
  formula)
  Behavior: `sum(t(R) %*% R) - sum(t(postF) %*% postF) + sum(postF2)`
    where R = Y - X*postF. The first two terms compute the FULL
    Gram-matrix sum of (T x T) cross-products, which equals
    `sum_i rowSums(R[i, ])^2` and `sum_j rowSums(postF[j, ])^2`
    respectively, NOT the diagonal traces `sum(R^2)` and
    `sum(postF^2)`. The third term `sum(postF2)` is the sum of
    the posterior variances (since `postF2 = postsd^2`) but is
    missing the `predictor_weights` factor `pw[j] = colSums(X^2)`
    that the variational decomposition requires.
    Denault flagged this himself with the inline TODO at line
    74: `"not correct, need to correct bottom part"`.
  Decision: replaced-by-correct-formula in
    `R/adjust_covariates.R::mf_residualize_wavelet_eb`:
    `er2 = sum((Y - X*postF)^2) + sum_j pw[j] * sum_t Var(beta[j,t])`
    where `pw[j] = colSums(X^2)` and `Var(beta[j,t]) = fitted_wc2[j,t]`.
    Pattern A bug fix per refactor-discipline section 3.
  Reason: the variational decomposition for the posterior-
    expected squared error under a factorized posterior is
    well-known; the upstream formula is dimensionally wrong.
    Fixing it in our port is the correct call per the
    "no-buggy-code" rule.

  Joint impact of the two bug fixes above: our
    `mf_adjust_for_covariates(method = "wavelet_eb")` output
    is no longer bit-identical to upstream `EBmvFR(adjust = TRUE)`.
    The deviation is the magnitude of the upstream bug. The
    `tests/testthat/test_adjust_covariates.R` Pattern A test
    asserts (i) closed-form closeness for the OLS path, (ii)
    structural correctness for the wavelet-EB path (positive
    sigma2, posterior shrinkage on near-null effects), and
    (iii) numerical agreement with upstream within the bug-
    magnitude tolerance documented per-scenario.

- mvf.susie.alpha/R/multfsusie.R::multfsusie.obj
  Behavior: legacy mid-run resume via `multfsusie.obj = m_prev`,
    coupled with a separate `max_step` per-call iteration cap.
  Decision: replaced-by-`R/ibss_methods.R::ibss_initialize.mf_individual`
  Reason: mfsusieR delegates to susieR's IBSS workhorse, which
    already supports warm-start through `model_init`. The
    mfsusieR S3 override now honors `params$model_init` by
    copying `alpha`, `mu`, `mu2`, `V`, `pi_V`, `G_prior`,
    `sigma2`, `fitted`, and `intercept` from the supplied
    fit before the IBSS loop runs, and resetting `KL`/`lbf`
    to NA. The legacy `multfsusie.obj` plus `max_step`
    workaround is retired in favor of the single
    `model_init` argument and the single `max_iter` budget.
    The L_greedy cross-round path threads `model_init` from
    one round to the next; `expand_model_init_to_L` mirrors
    `prune_single_effects` from susieR for the mfsusieR
    list-of-list `mu` / `mu2` shape, appending zero-state
    effect slots up to the requested L.

- mvf.susie.alpha/R/multfsusie.R::cor_small +
  mvf.susie.alpha/R/computational_routine.R::log_BFu (df branch)
  Behavior: Johnson 2005 scaled Student-t marginal Bayes factor
    for the per-variable SER, used as a small-sample
    correction to the Wakefield Normal marginal. df = n - 1
    per outcome.
  Decision: replaced-by-`mfsusie(small_sample_correction =
    TRUE)` -> `R/posterior_mixture.R::mixture_log_bf_per_scale_johnson`.
  Reason: the upstream `cor_small` argument name is opaque;
    the renamed argument describes the purpose. The kernel is
    ported R-only (no C++ change) and gates on
    `requireNamespace("LaplacesDemon")` via DESCRIPTION
    Imports. The fidelity test asserts per-variable LBFs from
    the ported kernel summed across scales match
    `fsusieR::log_BF(..., df = n - 1)` at tolerance <= 1e-12.

  Alternative considered: porting susieR's NIG path
    (`estimate_residual_method = "NIG"`). The Stage 4a audit
    established that every NIG hook in susieR is shadowed by
    an `mf_individual` S3 override, so plumb-through would
    have no effect, and that a wavelet-mixture-prior NIG
    marginal BF would be a separate research contribution
    (mixture x NIG composition; no manuscript derivation
    exists). Johnson-t is the implementable port that
    addresses the same use case.

### Cross-package audit (post hooks, 2026-04-30)

Phase A reports: `inst/notes/sessions/2026-04-30-cross-package-audit-{a,b,c}-*.md`.
Triage summary: `inst/notes/cross-package-audit-summary-posthooks.md`.
OpenSpec follow-ups: see summary for the V-semantics-cluster and
audit-followups change names. The three add-to-ledger entries from
this round follow.

- fsusieR/R/EM.R:67-71 (`max_SNP_EM` top-K thinning of M-step input)
  vs mvf.susie.alpha/R/EM.R:67-71
  Behavior: upstream caps the M-step input to the top
    `max_SNP_EM = 100` SNPs by `lBF`. Under heavy LD the cutoff
    varies per IBSS iter and modality.
  Decision: replaced-by-alpha-threshold in
    `R/individual_data_methods.R::optimize_prior_variance.mf_individual`
    via `keep_idx <- which(zeta_l > params$alpha_thin_eps)`
    (`alpha_thin_eps = 5e-5` default).
  Reason: alpha-driven thinning is invariant to lBF-scale shifts
    across (m, s) groups, so the per-(outcome, scale) M step in
    mfsusieR sees a consistent variant set per outer iter regardless
    of which modality or scale is being solved. Top-K thinning by
    lBF is incompatible with the per-(m, s) split because each (m,
    s) has its own lBF magnitude; running top-K independently per
    (m, s) would produce inconsistent variant sets across (m, s)
    within one effect's SER step. The alpha threshold is set well
    below SuSiE's `prior_tol = 1e-9` (the V-on-effect filter), so
    no truly-signal variant is dropped: an alpha < 5e-5 means the
    variant contributes < 0.005% to the SER posterior. Audit ID
    C-2 (2026-04-30).

- fsusieR/R/computational_functions.R:1604-1693 (`TI_regression.susiF`,
  lead-variant column) and 771-875 (`HMM_regression.susiF`)
  Behavior: upstream picks the lead variable via
    `which.max(obj$alpha[[l]])` and regresses on the single
    column `X[, idx[l]]` for the per-effect TI / HMM smoother.
  Decision: replaced-by-alpha-weighted-aggregate in
    `R/mfsusie.R:419-427` via
    `X_eff[[l]] <- X %*% (fit$alpha[l, ] * data$csd)`, consumed
    by `R/mfsusie_methods.R:738` (TI) and `R/mfsusie_methods.R:914`
    (HMM).
  Reason: the SuSiE variational posterior on the variable index
    is a categorical distribution with mass `alpha[l, ]`. The
    lead-variant approach commits to one variant; the alpha-
    weighted aggregate respects the posterior coverage uncertainty
    and gives the same answer when alpha is concentrated on a
    single variant. The two coincide on tightly-mapped SuSiE fits
    and diverge when the credible set is broad; the latter is
    where preserving posterior uncertainty matters. Audit ID B-8
    (2026-04-30). Phase 5 FDR investigation should compare both
    on the same fixtures before this design is locked.

- mvf.susie.alpha/R/EM.R:33 (`espsilon = 0.0001` separate inner-EM tol)
  Behavior: upstream maintains a separate `espsilon` for the inner
    EM convergence cutoff (`abs(newloglik - oldloglik) < espsilon`),
    distinct from the outer ELBO-change tol.
  Decision: shared-tol-in mfsusieR. The post-loglik hook reads
    `inner_tol <- params$tol %||% 1e-4`
    (`R/individual_data_methods.R:702`), reusing the IBSS outer
    convergence tol for the inner per-effect lbf-change cutoff.
  Reason: the user-facing surface for tolerance control is a
    single `tol` argument on `mfsusie()` (mfsusie.R:267); a
    separate inner_em_tol would proliferate knobs for a
    rarely-tuned semantic. Both checks compare per-iteration
    log-likelihood-proxy changes at the same scale, so a shared
    threshold is reasonable. If a Phase 5 / Phase 7 investigation
    surfaces a need to tune the two independently, the shared-tol
    decision can be split into a public `inner_em_tol` argument.
    Audit ID C-6 (2026-04-30).
