# Cross-package audit Agent B — mfsusieR vs fsusieR (2026-04-30)

## Summary
- Total findings: 9 (by severity: fix-now=0, track-for-later=4, needs-trace=2, accept-pre-existing=3)
- Pathways audited: data prep / prior init / SER step / M-step / inner-EM / variant thinning / residual variance / ELBO / smoothers / CS / PIP
- Agent confidence: data prep = high; prior init = high; SER step = high; M-step = high; inner-EM = high; variant thinning = high; residual variance = medium (mfsusieR's per-scale path has no fsusieR analogue, only smoke comparison possible); ELBO = high; smoothers = high (TI / smash / HMM all read carefully); CS = high; PIP filter = high.

## Findings

### Finding B-1: mixsqp control defaults differ between packages
- **Severity**: track-for-later
- **Pathway**: M-step
- **Location**: mfsusieR/R/em_helpers.R:143-146 ; fsusieR/R/susiF.R:287-290
- **Verified by reading**: yes
- **Description**: fsusieR's default `control_mixsqp = list(verbose=FALSE, eps=1e-6, numiter.em=40)` (susiF.R:287-290). mfsusieR's `mf_em_m_step_per_scale` default `default_ctrl <- list(verbose=FALSE, tol.svd=0, numiter.em=10L, convtol.sqp=1e-6)` (em_helpers.R:143-146). The two differ on three axes: (a) mfsusieR sets `tol.svd = 0` to skip an irlba SVD that helps low-rank L matrices but adds cost on full-rank ones (mfsusieR's L is full-rank by construction); (b) mfsusieR uses 4x fewer EM-inside-mixsqp iterations (10 vs 40); (c) mfsusieR uses `convtol.sqp = 1e-6` explicitly while fsusieR relies on mixsqp's default. mfsusieR amortizes the lower per-call iteration count by warm-starting from `pi_prev` (em_helpers.R:153-159, 179-189), an option fsusieR does not use. The tighter cold-start budget is intentional but is a noticeable single-shot cost difference at iter 1 (no warm start available).
- **Triage hint**: track-for-later
- **Suggested action**: Add a comment to `mf_em_m_step_per_scale` explaining that `numiter.em = 10` is paired with the warm-start; if the warm-start branch is ever bypassed (e.g., a future refactor of the iter cache), revisit the iteration cap.

### Finding B-2: mfsusieR's PIP filter now has a meaningful gate that fsusieR lacks
- **Severity**: track-for-later
- **Pathway**: PIP filter and refinement
- **Location**: mfsusieR/R/individual_data_methods.R:725, 734-748, 1001-1003 ; mfsusieR/R/mfsusie.R:358 ; susieR/R/susie_get_functions.R:494-498 ; fsusieR/R/operation_on_susiF_obj.R:1562-1575
- **Verified by reading**: yes
- **Description**: After commit 605f395 the post-loglik hook returns `V = .effective_V_l(...)` (the per-effect mean over `(m, s)` of `sum_k pi[l, m, s, k] * var_k`, individual_data_methods.R:725, 734-748). `cleanup_extra_fields.mf_individual` returns `character(0)` (1001-1003) so V survives the cleanup. mfsusieR's `params$prior_tol = 1e-9` (mfsusie.R:358) feeds `susie_get_pip(model, prior_tol = 1e-9)` (susieR/R/susie_get_functions.R:494-498), which drops effects with `V[l] < 1e-9`. Effects whose mixture pi has collapsed to all-null see `var_k = 0` for the only positive-weight component and `V[l] = 0`, so they are excluded from the PIP product. fsusieR's `update_cal_pip.susiF` (operation_on_susiF_obj.R:1562-1575) computes `1 - prod_l(1 - alpha[l, j])` over every effect with no V-based filter; collapsed effects still contribute their (typically uniform 1/p) alpha to the PIP product. The mfsusieR behavior is the more principled choice and matches susieR's scalar-V convention; the difference will surface as different PIPs whenever any effect collapses, particularly in low-signal scenarios.
- **Triage hint**: track-for-later
- **Suggested action**: Document the gate explicitly in the `mfsusie()` return-value docstring (alongside the existing description of `pip`). Phase 5 FDR investigation should compare mfsusieR's filtered PIP to the pre-605f395 unfiltered PIP on a benchmark with collapsing effects.

### Finding B-3: Wakefield Bayes factor formulation diverges (joint sigma2 vs marginal residual variance)
- **Severity**: accept-pre-existing
- **Pathway**: SER step
- **Location**: mfsusieR/R/individual_data_methods.R:120-146 (joint sigma2 path) ; fsusieR/R/computational_functions.R:69-165 (marginal residual SE)
- **Verified by reading**: yes
- **Description**: This is the documented Divergence 3 from the prompt's KNOWN DIVERGENCES section. mfsusieR's `mf_per_outcome_bhat_shat` reads `shat2_m` from the per-outcome cache built by `refresh_iter_cache.mf_individual` as `outer(1/pw, sigma2_per_pos)` (individual_data_methods.R:512), the susieR convention. fsusieR's `cal_Bhat_Shat` (computational_functions.R:69-165, full-Y branch lines 85-101) computes `Shat[k, j] = sd(Y[, j] - X[, k] * Bhat[k, j]) / sqrt(n - 1)` per `(j, k)` pair using the residual after an OLS fit, the marginal-residual convention. mfsusieR uses the joint variance from the running SER residual; fsusieR uses the per-pair marginal residual. Per the prompt's KNOWN DIVERGENCES, this is intentional and not a fix-now finding; mentioned once here per the prompt's instruction.
- **Triage hint**: accept-known
- **Suggested action**: None.

### Finding B-4: `mfsusie` warns and emits a `status_ok` check that fsusieR's m_step omits
- **Severity**: track-for-later
- **Pathway**: M-step
- **Location**: mfsusieR/R/em_helpers.R:191-199 ; fsusieR/R/EM.R:333-343
- **Verified by reading**: yes
- **Description**: mfsusieR's `mf_em_m_step_per_scale` checks `res$status` after the mixsqp call and emits a `warning_message(...)` when mixsqp does not converge to optimal (em_helpers.R:191-199), with a hint suggesting a stronger `control_mixsqp` setting. fsusieR's `m_step.lik_mixture_normal` (EM.R:333-343) and `scale_m_step` (EM.R:429-441) do not inspect `res$status` and silently return `res$x`. mfsusieR's behavior is more conservative; the warning surfaces a real signal (mixsqp cold-start failure) but is per-(outcome, scale) per inner-EM cycle, which can be noisy on a hard problem. Not a divergence to fix; flagging because the warning behavior is a deliberate user-visible difference from fsusieR.
- **Triage hint**: track-for-later
- **Suggested action**: Consider a once-per-fit aggregation if mixsqp non-convergence proves a frequent vignette annoyance.

### Finding B-5: CS construction routes through susieR, not the fsusieR cumsum-based path
- **Severity**: needs-trace
- **Pathway**: CS construction
- **Location**: mfsusieR/R/mfsusie.R:393 (delegation to `susie_workhorse`) ; mfsusieR/R/mfsusie.R:361-364 (coverage / min_abs_corr forwarded) ; fsusieR/R/operation_on_susiF_obj.R:1601-1636
- **Verified by reading**: yes
- **Description**: mfsusieR's CS comes from `susie_workhorse(data, params)` -> `ibss_finalize.default` -> `get_cs(...)` -> `susie_get_cs(...)`, the susieR backbone path. It uses susieR's threshold + min-abs-corr purity filter. fsusieR's `update_cal_cs.susiF` (operation_on_susiF_obj.R:1601-1636) uses the simpler `cumsum(sort(alpha, decreasing=TRUE))` until coverage is met, with `check_cs.susiF` (operation_on_susiF_obj.R:196-211) running the purity filter after the fact. The threshold-based selection logic is the same (sort by alpha, accumulate until cumsum > coverage). mfsusieR's `n_purity = 100` (mfsusie.R:363) and `min_abs_corr = 0.5` (mfsusie.R:262) match fsusieR's `min_purity = 0.5` (susiF.R:283). What I cannot verify from reading alone: whether the susieR path applies the purity filter at the same point in the pipeline as fsusieR (susieR may filter first and not include purity-failed CSs at all; fsusieR keeps then drops). Smoke testing a fsusieR-vs-mfsusieR fit on a controlled scenario would resolve this.
- **Triage hint**: needs-trace
- **Suggested action**: A C2 fidelity test on a small fixture should compare `fit$sets$cs` between mfsusieR (`fsusie(Y, X, pos)`) and `fsusieR::susiF(Y, X, pos)` and document the per-CS purity / coverage rounding.

### Finding B-6: HMM `mu` subset patch is correctly applied in mfsusieR
- **Severity**: accept-pre-existing
- **Pathway**: smoothers
- **Location**: mfsusieR/R/post_smooth_hmm.R:170-173
- **Verified by reading**: yes
- **Description**: mfsusieR's `mf_fit_hmm` correctly subsets `mu <- mu[idx_comp]` after `prob[, idx_comp]` and `P[idx_comp, idx_comp]` (post_smooth_hmm.R:170-173, with explicit comment citing the refactor-exceptions ledger). fsusieR's `fit_hmm` was missing the `mu <- mu[idx_comp]` line (computational_functions.R:600-604) at the time of the original audit; per the refactor-exceptions ledger, fsusieR PR #31 has since restored it. The bug fix is documented in the ledger and the mfsusieR test asserts uniform bit-identity across both regimes.
- **Triage hint**: accept-known
- **Suggested action**: None.

### Finding B-7: smash post-smoother credible band derives sd from the band half-width
- **Severity**: track-for-later
- **Pathway**: smoothers
- **Location**: mfsusieR/R/post_smooth_smash.R:56-61 ; fsusieR/R/computational_functions.R:1858-1920
- **Verified by reading**: yes
- **Description**: mfsusieR's `.post_smooth_smash` recovers per-position sd from `(out$cred_band[1L, ] - out$cred_band[2L, ]) / (2 * z_crit)` (post_smooth_smash.R:58-59) and uses it for lfsr derivation. The fsusieR equivalent (`smash_regression.susiF`, computational_functions.R:1858-1920) does not produce per-position lfsr — it only computes the credible band and stores `fitted_var` for later cred_band reconstruction. mfsusieR's lfsr derivation is a strict superset of fsusieR's smash output and uses the same posterior mean / sd algebra. The reverse derivation `(up - low) / (2 * z_crit) = sd` is correct algebraically only when the band is symmetric; for `smash.gaus` with `post.var = TRUE` the band is built with `coeff * sqrt(var)`, so this is exact. Not a divergence; flagging because the lfsr derivation is mfsusieR-only and downstream consumers may not realize fsusieR's smash output had no lfsr.
- **Triage hint**: track-for-later
- **Suggested action**: None other than ensuring the mfsusieR smash path's lfsr is documented as a feature add over fsusieR.

### Finding B-8: TI / HMM smoother regression input differs from fsusieR (alpha-weighted X aggregate vs lead-variant column)
- **Severity**: needs-trace
- **Pathway**: smoothers
- **Location**: mfsusieR/R/mfsusie.R:419-427 (X_eff construction) ; mfsusieR/R/mfsusie_methods.R:738 (TI consumption) ; mfsusieR/R/mfsusie_methods.R:914 (HMM consumption) ; fsusieR/R/computational_functions.R:1604-1693 (TI lead-variant) ; fsusieR/R/computational_functions.R:771-875 (HMM lead-variant)
- **Verified by reading**: yes
- **Description**: mfsusieR's post-fit `X_eff[[l]] = X %*% (fit$alpha[l, ] * data$csd)` (mfsusie.R:425-427) is an alpha-weighted aggregate across all variables for each effect. The TI / HMM smoothers regress the per-effect isolated response on `fit$X_eff[[l]]` (mfsusie_methods.R:738, 914). fsusieR's `TI_regression.susiF` and `HMM_regression.susiF` (computational_functions.R:1613-1618, 781-786) instead pick the lead variable by `which.max(obj$alpha[[l]])` and regress on `X[, idx[l]]` (a single column). The two approaches give identical results when `alpha[l, ]` concentrates near a Dirac on one variable; they diverge when the CS has spread mass. The alpha-weighted aggregate captures the SuSiE posterior coverage uncertainty; the lead-variant approach commits to one variant. Not a bug in either direction, but a distinct algorithmic choice. Phase 5 FDR investigation should compare both on the same fixtures.
- **Triage hint**: needs-trace
- **Suggested action**: Add a vignette-level note documenting the alpha-weighted regressor choice. Optionally provide a `lead_variant_only = TRUE` mode on `mf_post_smooth(method = "TI"/"HMM")` for backward-compatible comparisons.

### Finding B-9: ELBO short-circuits to NA on PIP convergence path
- **Severity**: accept-pre-existing
- **Pathway**: ELBO
- **Location**: mfsusieR/R/individual_data_methods.R:973-977 ; mfsusieR/R/mfsusie.R:403 ; fsusieR/R/ELBO.R:211-217 ; fsusieR/R/susiF.R:292
- **Verified by reading**: yes
- **Description**: mfsusieR's `get_objective.mfsusie` short-circuits to `NA_real_` when `params$convergence_method == "pip"` (individual_data_methods.R:973-975), the default for `mfsusie()` (`convergence_method = c("pip", "elbo")`, mfsusie.R:266). fsusieR's `get_objective.susiF` (ELBO.R:211-217) returns `sum(obj$KL)` (note: just sum of KL, not Eloglik - sum(KL); fsusieR's stated ELBO is incomplete per the documented Divergence 2). fsusieR's `cal_obj = FALSE` default (susiF.R:292) means fsusieR also skips ELBO evaluation by default. Both packages skip ELBO computation by default; mfsusieR is principled (full ELBO computed only on demand), fsusieR is incidental (the get_objective body is wrong anyway). mfsusieR's iter-1 ELBO is masked to NA (mfsusie.R:403) which is correct (the iter-1 sigma2 is the inflated `var(Y)` initial). No divergence to fix; mentioned for completeness.
- **Triage hint**: accept-known
- **Suggested action**: None.

## Pathways with no findings

- **Data prep** (silent because-clean). mfsusieR's `mf_dwt` (R/dwt.R:56-132) replicates fsusieR's `remap_data + colScale + DWT2 + cbind(D, C) + colScale + gen_wavelet_indx` pipeline (R/susiF.R:372-413). The differences are documented in the refactor-exceptions ledger PR group 2 entries (vectorized `csd[csd == 0] <- 1`, simplified `d` attribute, `seq_len(lev_res - 1)`, single-step interpol_mat, `min.scale` parameterization). Defaults match: `wavelet_standardize = TRUE` matches fsusieR's `colScale(Y_f)` (always on); `wavelet_qnorm = FALSE` matches fsusieR's `quantile_trans = FALSE` default; `max_padded_log2 = 10` matches fsusieR's `max_scale = 10`. The lowc_idx filtering deferral is documented in the ledger.

- **Prior init** (silent because-clean). The data-driven path divergence (`init_scale_mixture_prior_default` deterministic 10th-percentile-Shat anchor vs fsusieR's `init_prior.default` ash-driven init) is logged in the refactor-exceptions ledger PR group 4 entries with derivation. The cold-start `null_prior_init` parameter (mfsusieR's preferred-position name; the public arg in `mfsusie()` is `null_prior_init = 0` per mfsusie.R:251 and prior_scale_mixture.R:89) is parameterised whereas fsusieR hardcodes `pi[1] = 0.8`; per the ledger this is intentional. The new ebnm-backed prior classes (per_scale_normal / per_scale_laplace) are mfsusieR-only; they have no fsusieR analogue and the M-step kernel correctly locks `mode = 0` (individual_data_methods.R:858-861). Init pi via `lead_init_s = which.max(rowMeans(bhat^2))` (prior_scale_mixture.R:202-208) is sane.

- **SER step** (silent because-clean). mfsusieR's `compute_residuals.mf_individual` -> `compute_ser_statistics.mf_individual` -> `loglik.mf_individual` -> `calculate_posterior_moments.mf_individual` chain (individual_data_methods.R:28-308) replicates fsusieR's `cal_partial_resid -> cal_Bhat_Shat -> log_BF -> post_mat_mean / post_mat_sd` chain (operation_on_susiF_obj.R:37-82, computational_functions.R:69-1349). The Bhat / Shat divergence is the documented Divergence 3 (Finding B-3). `mixture_log_bf_per_scale` (posterior_mixture.R:35-65) implements the same mixture marginal as `log_BF.mixture_normal_per_scale` (computational_functions.R:1028-1135) using `logSumExp` over components. Both packages clamp `Shat[Shat<=0] <- 1e-32` (computational_functions.R:1038) — mfsusieR does this implicitly via the cached `shat2` floor. The cross-outcome combiner is mfsusieR-only (not relevant for M=1).

- **M-step** (silent because-clean other than B-1, B-4). mfsusieR's `mf_em_likelihood_per_scale` + `mf_em_m_step_per_scale` (em_helpers.R:48-202) replicate fsusieR's `cal_L_mixsq_s_per_scale` + `scale_m_step` (EM.R:157-185, 405-443). Both build the `(p * |idx_s| + 1) x K` mixsqp likelihood matrix with the same `c(100, 0, ..., 0)` penalty row prepend (em_helpers.R:83 vs EM.R:180-181). Both build the same weight vector `c(nullweight * |idx_s|, rep(zeta, |idx_s|))` (em_helpers.R:136 vs EM.R:418-420). The `tol_null_prior` collapse path is the same (em_helpers.R:200 vs EM.R:437-440). For M=1 the `mixture_null_weight * max(1L, data$M) = mixture_null_weight` is unchanged from the per-outcome design (this is the documented cross-package fix-now bundle in the ledger).

- **Inner-EM control surface** (silent because-clean). mfsusieR's `max_inner_em_steps = 5L` default (mfsusie.R:282) plus the off-by-one cap `inner_cap = max(0L, max_inner_em_steps) + 1L = 6L` (individual_data_methods.R:701) matches the prompt's description: 6 M-step cycles per effect per outer iter. Early exit on `abs(model$lbf[l] - prev_lbf) < inner_tol` (individual_data_methods.R:723) maps to fsusieR's `EM_pi` early exit on `abs(newloglik - oldloglik) >= espsilon` (EM.R:99). fsusieR's `max_step_EM = 1` (susiF.R:300) gives 2 cycles per outer iter (loop runs `k <= max_step` so 2 cycles). The semantics match; the default values differ as documented in the prompt.

- **Variant thinning** (silent because-clean). mfsusieR's `alpha_thin_eps = 5e-5` default (mfsusie.R:279) drops variables with `alpha[l, j] < 5e-5` from the M-step input (individual_data_methods.R:587-590). fsusieR's `max_SNP_EM = 100` default (susiF.R:299, EM.R:80-85) keeps the top-100 variables by lbf. Different mechanism (threshold vs rank), different defaults; both effective at capping the L-matrix row count. Neither drops truly-signal variants in expectation: the alpha threshold drops variables with alpha < 5e-5 which is below SuSiE's own pruning threshold (`prior_tol = 1e-9` is on V not alpha, but alpha < 5e-5 means a variable is contributing < 0.005% to the SER posterior). The rank-based filter keeps up to 100 variables which on a typical p ~ 1000 fine-mapping is generous. The defaults are not directly comparable.

- **Residual variance update** (silent because-clean for `per_outcome` mode; mfsusieR's `per_scale` mode has no fsusieR analogue). mfsusieR's `update_variance_components.mf_individual` (individual_data_methods.R:458-477) for `per_outcome` matches fsusieR's `estimate_residual_variance.susiF` (operation_on_susiF_obj.R:302-307): both divide `sum(ER2)` by `n * T_basis`. The `per_scale` mode (default in mfsusieR per mfsusie.R:255-256) divides per-scale `sum(ER2[idx_s]) / (n * |idx_s|)` and is mfsusieR-only. The ER2 computation in `mf_get_ER2_per_position` (individual_data_methods.R:422-443) uses the corrected susieR-style `||Y - X*postF||^2 + sum_l ((alpha_l . pw)^T mu2_l - ||X * (alpha_l . mu_l)||^2)` decomposition; fsusieR's `get_ER2.susiF` (operation_on_susiF_obj.R:525-534) uses the buggy `sum(t(R) %*% R) - sum(t(postF) %*% postF) + sum(postF2)` form. The mfsusieR fix is documented in the refactor-exceptions ledger as the Pattern A bug fix for the multivariate version (entry on `mvf.susie.alpha/R/operation_on_multfsusie_obj.R:1046-1096`); the same fix applies to the single-modality M=1 path because the algebra is identical.

- **ELBO** (silent because-clean other than B-9). mfsusieR's `get_objective.mfsusie` returns `Eloglik(data, model) - sum(model$KL)` after a `refresh_lbf_kl` pass (individual_data_methods.R:973-977). The susieR-convention KL is computed in `compute_kl.mf_individual` (individual_data_methods.R:368-385) as `KL_l = -lbf_l - L_null + SER_posterior_e_loglik(l)`, which matches the documented Divergence 2 (susieR's correct identity, not mvf's sign-flipped one). For M=1 the only difference from fsusieR is that fsusieR's `get_objective.susiF` returns `sum(KL)` only and is essentially unused (cal_obj = FALSE default).

- **Smoothers** (silent because-clean other than B-6, B-7, B-8). The four smoothers in mfsusieR (`scalewise`, `TI`, `HMM`, `smash`) replicate the algebra of fsusieR's three (`TI`, `HMM`, `smash`) plus an extra `scalewise` which fsusieR lacks. The TI's stationary-wavelet decomposition + scalewise ash + cycle-spinning av-basis pipeline (mfsusie_methods.R:802-880) is line-for-line a refactor of `univariate_TI_regression` (computational_functions.R:1756-1857) using `compute_marginal_bhat_shat` instead of `cal_Bhat_Shat`. Both use `nullweight = 30` for D coefficients and `nullweight = 3` for C. The variance basis uses Daubechies 10 / DaubLeAsymm in both. The smash path matches fsusieR closely; the lfsr derivation is mfsusieR-only (B-7).

- **CS construction** (silent because-clean other than B-5).

- **PIP filter** (silent because-clean other than B-2).

## Limits of this audit

- I did not run any tests. All findings are from code reading.
- The covariate-adjustment pathway (`adjust_covariates.R`) is documented in the refactor-exceptions ledger and out of scope for this audit per the EBmvFR ruling.
- I did not audit the susieR backbone (Agent A territory) or the multi-modality combiner (Agent C territory). Cross-outcome combiner code in `prior_cross_outcome.R` is M >= 2 only and was not traced.
- Smoke testing fsusieR vs mfsusieR on the same fixture would resolve B-5 (CS purity / coverage rounding) and B-8 (alpha-weighted vs lead-variant smoother regressor) more conclusively. I read both implementations and confirmed the algebra; observable divergence depends on alpha concentration which I cannot quantify from code alone.
- I did not trace the L-greedy outer loop interaction with the new V-as-effective-slab-variance behavior. The `expand_model_init_to_L` helper (ibss_methods.R:120) initializes new effects' V to 1; on greedy expansion, the new effects' V will then be overwritten to the effective slab variance after their first post-loglik hook. This is not a bug, but the prompt mentioned PIP filter behavior under L-greedy is worth a smoke test.
