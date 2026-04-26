# Cross-package audit A: mfsusieR vs susieR

Date: 2026-04-26
Agent: Explore (opus), read-only.
Scope: numerical divergence + delegation misses, mfsusieR public API and S3 overrides against `susieR/R/`.

## 1. Public API surface

| Public symbol | Argument-name observation | Triage |
|---|---|---|
| `mfsusie()` | `prior_variance_grid`, `prior_weights`, `residual_variance`, `coverage`, `min_abs_corr`, `max_iter`, `tol`, `L`, `L_greedy`, `lbf_min`, `model_init`, `verbose`, `track_fit`, `intercept`, `standardize` consistent with susieR `susie()` | accept-as-documented |
| `mfsusie()` divergent names | `null_prior_weight` (vs susieR `null_weight`), `mixture_weight_method` (vs `estimate_prior_method`), `prior_variance_scope`, `residual_variance_scope`, `wavelet_filter_number`, `wavelet_family`, `low_count_filter`, `quantile_norm`, `mixsqp_null_penalty`, `small_sample_correction`, `max_padded_log2`. Wavelet-specific or scope-specific; no susieR equivalent. Translation layer maps `mixture_weight_method` to susieR's `estimate_prior_method` internally. | accept-as-documented |
| `fsusie()` | Single-outcome wrapper. Argument order `(Y, X, pos, ...)` swaps `(X, y)` of susieR; intentional fsusieR-port compatibility. | accept-as-documented |
| `predict.mfsusie` | `newx` matches susieR `predict.susie`. Returns `list[M]` of `n x T_m`. | accept-as-documented |
| `coef.mfsusie` | Adds `smooth_method` not in susieR. Documented opt-in. | accept-as-documented |
| `summary.mfsusie` | Returns `n_effects`/`n_variables`/`n_outcomes`/`T_basis`. susieR returns `vars`/`cs` data-frames. Different schema. | track-for-later (Phase 8 polish) |
| `print.mfsusie`, `fitted.mfsusie` | Custom format / list-valued. No susieR conformance issue. | accept-as-documented |
| `mf_post_smooth`, `mf_adjust_for_covariates`, `mf_simu_ibss_per_level`, `mfsusie_plot` | mfsusieR-specific, no susieR equivalents. | accept-as-documented |

## 2. S3 dispatch overrides on `mf_individual` / `mfsusie`

### 2.1 `ibss_initialize.mf_individual` (`R/ibss_methods.R:140-171`)
- Skips susieR's scalar-only residual_variance check (justified by per-(scale, outcome) shape).
- Custom `expand_model_init_to_L` instead of susieR's `validate_init` + `prune_single_effects` + `adjust_L`. Justified by list-of-list `mu`/`mu2` shape; ledgered.
- Skips slot_prior init; not a numerical divergence on the supported path.
- Triage: shape divergences accept-as-documented; slot_prior gap track-for-later.

### 2.2 `validate_prior.mf_individual` (`R/ibss_methods.R:184-186`)
- Body is `invisible(TRUE)`, byte-equivalent to `validate_prior.default`. Redundant override.
- Triage: track-for-later (cosmetic).

### 2.3 `track_ibss_fit.mf_individual` (`R/ibss_methods.R:200-209`)
- Records `pi_V` instead of scalar `V`, skips slot-activity. Justified.
- Triage: accept-as-documented.

### 2.4 `check_convergence.mf_individual` (`R/ibss_methods.R:220-227`)
- Tiny ELBO-only check; silently disables susieR's PIP-convergence and verbose tabular output.
- Triage: track-for-later (feature gap, not a bug).

### 2.5 `trim_null_effects.mf_individual` (`R/ibss_methods.R:241-243`)
- Pass-through. Justified (mfsusieR holds `V[l] = 1`). Documented inline.
- Triage: accept-as-documented.

### 2.6 `cleanup_model.mf_individual` (`R/ibss_methods.R:263-270`)
- Re-implements field-stripping rather than calling `NextMethod()`. Differences are inert (the un-stripped fields aren't populated in the mfsusieR path) but the duplication is a delegation miss.
- Triage: track-for-later.

### 2.7 `configure_data.mf_individual` (`R/ibss_methods.R:283-285`)
- Pass-through, byte-equivalent to default. Redundant.
- Triage: track-for-later (cosmetic).

### 2.8 `compute_residuals.mf_individual` (`R/individual_data_methods.R:28-51`)
- Uses `crossprod(X, R_m)` instead of susieR `compute_Xty`. Forecloses sparse-X support.
- Triage: track-for-later.

### 2.9 `compute_ser_statistics.mf_individual` (`R/individual_data_methods.R:71-93`)
- Same arithmetic as susieR up to per-outcome list shape. No undocumented divergence.
- Triage: accept-as-documented.

### 2.10 `loglik.mf_individual` (`R/individual_data_methods.R:159-212`)
- Genuinely different algorithm (mixture vs single-Gaussian). Stable-softmax block re-implements `lbf_stabilization` + `compute_posterior_weights`; comment explains the inlining (per-position zero-pw logic).
- Triage: accept-as-documented.

### 2.11 `calculate_posterior_moments.mf_individual` (`R/individual_data_methods.R:230-254`)
- Forwards to `mixture_posterior_per_scale` cpp11 kernel. Different algorithm; not a delegation miss.
- Triage: accept-as-documented.

### 2.12 `SER_posterior_e_loglik.mf_individual` (`R/individual_data_methods.R:275-297`)
- Per-position sigma2 generalization of susieR scalar form.
- Triage: accept-as-documented.

### 2.13 `compute_kl.mf_individual` (`R/individual_data_methods.R:312-328`)
- Per-position generalization of susieR Gaussian KL.
- Triage: accept-as-documented.

### 2.14 `update_variance_components.mf_individual` (`R/individual_data_methods.R:395-414`)
- Genuinely different shape; documented in refactor-exceptions (`get_ER2.multfsusie`).
- Triage: accept-as-documented.

### 2.15 `update_model_variance.mf_individual` (`R/individual_data_methods.R:488-492`)
- Roxygen truncated mid-sentence (`#' refinement.`).
- Triage: track-for-later (doc fix).

### 2.16 `update_derived_quantities.mf_individual`
- Rebuilds `model$fitted[[m]]` from scratch. Justified by mfsusieR's running-fit storage.
- Triage: accept-as-documented.

### 2.17 `Eloglik.mf_individual`, `get_objective.mf_individual`
- Per-position generalization. No slot_prior / omega.
- Triage: accept-as-documented.

### 2.18 `optimize_prior_variance.mf_individual`
- Mixsqp per (outcome, scale). Different algorithm; correct generic contract.
- Triage: accept-as-documented.

### 2.19 `update_fitted_values.mf_individual`
- Uses `X %*% b_lm` instead of susieR `compute_Xb`. Same delegation issue as 2.8.
- Triage: track-for-later.

## 3. Helpers vs susieR equivalents

### 3.1 `R/utils_wavelet.R`
- All helpers (`col_scale`, `power_of_two_gridn`, `interpol_ks`, `interpol_mat`, `remap_data`, `gen_wavelet_indx`, `dwt_matrix`, `wd_variance`, `wst_variance`, `av_basis_variance`, `mf_invert_variance_curve`, `mf_low_count_indices`, `mf_quantile_normalize`) are local ports; susieR has no equivalents.
- Triage: accept-as-documented.

### 3.2 `R/post_smooth_*.R`, `R/mfsusie_methods.R`
- `.post_smooth_smash` uses local `mf_per_position_bhat_shat` instead of `compute_marginal_bhat_shat`. Delegation miss.
- `mf_col_scale` (in `R/post_smooth_smash.R`) duplicates package-level `col_scale`. Delegation miss within mfsusieR.
- Triage: both track-for-later.

### 3.3 `R/em_helpers.R`, `R/posterior_mixture.R`, `R/prior_scale_mixture.R`, `R/prior_cross_outcome.R`
- All mfsusieR-specific. Calls to `compute_marginal_bhat_shat` are correct.
- Triage: accept-as-documented.

### 3.4 `R/adjust_covariates.R`
- Closed-form FWL has no susieR equivalent. Wavelet-EB uses `compute_marginal_bhat_shat`.
- Triage: accept-as-documented.

### 3.5 `R/simulation.R`
- mfsusieR-specific. accept-as-documented.

### 3.6 zzz.R cached internals
- `compute_Xb` and `compute_Xty` are not cached, despite being needed for sparse-X support (root cause of 2.8 / 2.19).
- Triage: track-for-later.

## 4. Refactor-exceptions ledger cross-check

All 24 ledger entries verified against current code. Two stale:

- **`lowc_wc filtering deferred`** entry says "v1 init does not filter low-count coefficients" and that C2 fidelity tests run with `lowc_wc = NULL`. Stale: `low_count_filter` is now a public `mfsusie()` argument and the prior init excludes flagged columns from the ash sampling pool (`R/prior_scale_mixture.R:69-84`).
- **`fit_hmm` HMM-mu-subset** entry's "Resolution (2026-04-26): fsusieR PR #31 restored the line" suggests upstream caught up. Test-status ambiguous (does test still assert deviation, or now bit-identity?).

Both: track-for-later.

## Final triage

| ID | Severity | Decision |
|---|---|---|
| F1: validate_prior redundant override | low | track-for-later |
| F2: check_convergence disables PIP-convergence | medium | track-for-later |
| F3: cleanup_model duplicates susieR field-stripping | low | track-for-later |
| F4: configure_data redundant override | low | track-for-later |
| F5: ibss_initialize skips slot_prior | medium | track-for-later |
| F6: compute_residuals/update_fitted_values skip susieR helpers | low | track-for-later |
| F7: update_model_variance roxygen truncated | low | track-for-later |
| F8: mf_per_position_bhat_shat duplicates compute_marginal_bhat_shat | low | track-for-later |
| F9: mf_col_scale duplicates col_scale | low | track-for-later |
| F10: summary.mfsusie schema diverges | low | track-for-later (Phase 8) |
| F11: lowc_wc filtering ledger entry stale | low | track-for-later |
| F12: HMM-mu-subset ledger entry status ambiguous | low | track-for-later |
| F13: compute_Xb/Xty not cached in zzz.R | low | track-for-later |

**No `fix-now` items. No undocumented numerical divergences.**
