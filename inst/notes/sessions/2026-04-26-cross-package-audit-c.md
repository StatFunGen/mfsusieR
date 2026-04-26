# Cross-package audit C: mfsusieR vs mvf.susie.alpha (multi-modality functional)

Date: 2026-04-26
Agent: Explore (opus), read-only.
Scope: M >= 1, T_m can vary, legacy variance mode (C3 contract). Numerical equivalence at the documented tolerances modulo divergences listed in `inst/notes/refactor-exceptions.md`.

## 1. Multi-modality SER kernel

### 1.1 Outcome-only multifsusie SER kernel matches at 1e-12 (accept-as-documented)
- mfsusieR: `R/individual_data_methods.R:159-212` (`loglik.mf_individual`).
- upstream: `mvf.susie.alpha/R/EM.R:42-50`, `R/computational_routine.R:259-334`.
- Reference test `tests/testthat/test_mvf_alpha_fidelity.R:23-76` asserts at `1e-12`.

### 1.2 `max_SNP_EM` truncation absent (track-for-later)
- mvf.susie.alpha truncates the SER step to top `max_SNP_EM = 100` SNPs by lBF (`R/EM.R:67-71`).
- mfsusieR runs mixsqp on all p SNPs unconditionally (`R/individual_data_methods.R:436-473`).
- Defensible (mixsqp is linear in p) but no ledger entry.

### 1.3 K-scaled `nullweight` not propagated (fix-now)
- mvf.susie.alpha scales `nullweight` by `K_f + K_u` (`R/EM.R:58-65, 102`).
- mfsusieR passes `mixsqp_null_penalty %||% 0.7` unscaled by M.
- Upstream rationale (inline comment): with M modalities, joint lBF concentrates ~M-fold, so null regularization must scale.
- Potential FDR-calibration source for M > 1.
- **Investigate**: does mfsusieR's per-(outcome, scale) M-step (independent solves, not joint) need the M-scaling? If yes, fix; if no, document the deliberate omission.

### 1.4 Inner-EM loop collapsed to single-pass (track-for-later)
- mvf.susie.alpha has an inner EM loop in `EM_pi_multsusie` (`R/EM.R:81-120`) with `max_step` iterations.
- mfsusieR collapses to one mixsqp solve per IBSS iteration. Consistent with susieR pattern.
- Add ledger entry confirming this is intentional.

## 2. Cross-outcome combine

### 2.1 Default combiner matches `apply(., 2, sum)` semantics (accept-as-documented)
- `R/prior_cross_outcome.R:34-37` defines `combine_outcome_lbfs.mf_prior_cross_outcome_independent` as `Reduce("+", outcome_lbfs)`.
- Equivalent to upstream `colSums(do.call(rbind, per_trait_lbfs))` up to FP summation order. Test `1e-15`.

### 2.2 u-trait/f-trait merge unified into single M-outcome path (accept-as-documented)
- mvf.susie.alpha sums `f_logBF + u_logBF`.
- mfsusieR collapses both into "M outcomes, T_m can vary"; T_m == 1 reproduces a univariate trait. Consistent with `inst/CLAUDE.md` hard rule #2.

### 2.3 Audit prompt mentioned combiner in prior init; misread
- `combine_outcome_lbfs` is only called from the SER step. Prior init runs per-modality and does not combine.

## 3. Per-(outcome, scale) prior

### 3.1-3.3: `mf_prior_scale_mixture` matches upstream (accept-as-documented)
- Replicated single ash fit per modality, hardcoded `(0.8, 0.2/(K-1), ...)` overwrite, parameter rename `gridmult` -> `grid_multiplier`. All matches upstream behaviour or covered by ledger.
- Test `1e-12` against `fsusieR::init_prior.default`.

## 4. Variance / ELBO

### 4.1 ER2 Pattern A bug fix - load-bearing (accept-as-documented)
- mvf.susie.alpha `R/operation_on_multfsusie_obj.R:1075-1090` has a dimensionally-wrong formula (missing `predictor_weights` factor; aggregates effects via `sum(postF^2)` sum-then-square).
- mfsusieR `R/individual_data_methods.R:360-380` (`mf_get_ER2_per_position`) uses the correct decomposition.
- Ledger entry `inst/notes/refactor-exceptions.md:251-276` documents this as Pattern A. Verified current.
- Pattern A regression guard: `test_mvf_alpha_degeneracy.R:114-129` `expect_gt(max(abs(s2_v - s2_m)), 1e-4)`.

### 4.2 Per-scale residual variance is new mode (accept-as-documented)
- mfsusieR adds `per_scale` variance mode. `per_outcome` reproduces upstream legacy mode.
- Ledger entry `:278-289` (PR group 6) documents.

### 4.3 `Eloglik` aggregate matches at `1e-12` (accept-as-documented)
- Per-position generalization. Test `test_variance_and_elbo.R:153-166`.

### 4.4 `get_objective` omits upstream `o` entropy term (fix-now)
- mvf.susie.alpha `R/ELBO_mutlfsusie.R:269-282`: `out <- Eloglik - sum(KL) + o`, where `o = Reduce("sum", ... alpha_l * log(p / pmax(alpha_l, 1e-6)))` — categorical-posterior entropy correction.
- mfsusieR `R/individual_data_methods.R:548-551`: `Eloglik - sum(KL)`. No `o` term.
- mfsusieR follows susieR's clean ELBO decomposition (`compute_kl.mf_individual` folds the categorical KL into the lbf-derived loglik aggregate).
- Either upstream's `o` is double-counting (Pattern A) or it represents a different KL convention (Pattern B).
- **Add ledger entry** documenting the deliberate omission with the susieR-derivation justification, citing the manuscript ELBO derivation.

### 4.5 Stale docstring fragment in `update_model_variance.mf_individual` (track-for-later)
- `R/individual_data_methods.R:480-483` mentions `params$residual_variance_lowerbound` "for refinement" but the field is not consulted.

## 5. Post-processing

### 5.1 HMM-mu-subset Pattern A fix preserved (accept-as-documented)
- `R/post_smooth_hmm.R:166` `mu <- mu[idx_comp]`. Ledger entry `:291-323` current.
- Test `test_post_smooth_HMM.R:56-94` at `tol = 0` against `fsusieR:::fit_hmm`.

### 5.2 Lazy IDWT vs eager (track-for-later)
- mfsusieR makes IDWT reconstruction lazy (called via `coef`/`predict`/`plot`); upstream `multfsusie` does it eagerly inside `out_prep`. No numerical impact.

### 5.3 Smash branch uses local helper instead of `compute_marginal_bhat_shat` (fix-now)
- `R/post_smooth_smash.R:108, 153-160` uses local `mf_per_position_bhat_shat`.
- HMM branch correctly uses `compute_marginal_bhat_shat` (`R/post_smooth_hmm.R:286`).
- Replace with the susieR helper; same closed form, removes a duplicate.

### 5.4 `mf_col_scale` duplicates `col_scale` (track-for-later)
- `R/post_smooth_smash.R:140-149` duplicates `R/utils_wavelet.R::col_scale`. Bundle with 5.3 in same edit pass.

### 5.5 Smash `mixVBEM` ashparam matches reference at `1e-12` (accept-as-documented)
- Test `test_post_smooth_smash.R:12-32`.

### 5.6 Cross-modality residual reconstruction algebraically equivalent (accept-as-documented)
- mfsusieR caches `D = residuals + fitted` in wavelet domain and inverts; upstream sums in position space. Same answer up to wavelet roundtrip noise.

## 6. Reference test fidelity

### 6.1 Tolerances match contract (accept-as-documented)
- `test_mvf_alpha_fidelity.R`: `1e-12` per-(modality) log-BF, `1e-15` cross-modality combine, `1e-12` posterior. Within C3 floor (`<= 1e-8`) and Pattern A/B convention.
- `test_mvf_alpha_degeneracy.R`: structural tests with documented bug-magnitude bounds.

### 6.2 No end-to-end IBSS-step fidelity tests (track-for-later)
- Per-component tests pass; full `alpha`/`mu`/`mu2`/`pip` end-to-end is blocked by sigma2 divergence (Pattern A).
- Could extend to per-IBSS-step bit-identity (first SER step before any sigma2 update).

### 6.3 PIP-within-0.05 structural test (accept-as-documented)
- `test_mvf_alpha_degeneracy.R:102-110`. Graduated tolerance is acceptable as a smoke threshold per `refactor-discipline.md` section 3 (apple-to-orange).

## Final triage

| ID | Severity | Decision |
|---|---|---|
| C-1.1 | n/a | accept-as-documented |
| C-1.2 | low | track-for-later |
| C-1.3 | medium | fix-now (potential FDR impact) |
| C-1.4 | n/a | track-for-later (doc) |
| C-2.1 | n/a | accept-as-documented |
| C-2.2 | n/a | accept-as-documented |
| C-2.3 | n/a | accept-as-documented |
| C-3.1 | n/a | accept-as-documented |
| C-3.2 | n/a | accept-as-documented |
| C-3.3 | n/a | accept-as-documented |
| C-4.1 | n/a | accept-as-documented |
| C-4.2 | n/a | accept-as-documented |
| C-4.3 | n/a | accept-as-documented |
| C-4.4 | medium | fix-now (ledger entry) |
| C-4.5 | low | track-for-later (doc) |
| C-5.1 | n/a | accept-as-documented |
| C-5.2 | low | track-for-later |
| C-5.3 | low | fix-now |
| C-5.4 | low | track-for-later |
| C-5.5 | n/a | accept-as-documented |
| C-5.6 | n/a | accept-as-documented |
| C-6.1 | n/a | accept-as-documented |
| C-6.2 | medium | track-for-later |
| C-6.3 | n/a | accept-as-documented |

**Three `fix-now` items: C-1.3 (K-scaled nullweight investigation), C-4.4 (`get_objective` ledger entry), C-5.3 (smash Bhat/Shat helper delegation).**
