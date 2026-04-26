# Cross-package audit B: mfsusieR vs fsusieR (M = 1 functional case)

Date: 2026-04-26
Agent: Explore (opus), read-only.
Scope: M = 1, T_1 > 1; numerical equivalence at the documented tolerances modulo divergences listed in `inst/notes/refactor-exceptions.md`.

## 1. Wavelet pipeline

| mfsusieR routine | mfsusieR location | fsusieR equivalent | Verdict |
|---|---|---|---|
| `is_wholenumber` | `R/utils_wavelet.R:17-19` | `is.wholenumber` (`R/utils.R:21-22`) | identical (cosmetic rename) |
| `mf_low_count_indices` | `R/utils_wavelet.R:40-46` | `which_lowcount` (`R/utils.R:9-18`) | empty result `integer(0)` vs `NULL` (cosmetic) |
| `mf_quantile_normalize` | `R/utils_wavelet.R:63-76` | `Quantile_transform` | bit-identical |
| `col_scale` | `R/utils_wavelet.R:96-137` | `colScale` (`R/utils.R:85-141`) | values + center/scale `tol = 0`; `d` attr `1e-12` (Pattern B; ledgered) |
| `interpol_ks` | `R/utils_wavelet.R:171-173` | `interpolKS2` | identical |
| `interpol_mat` | `R/utils_wavelet.R:199-213` | `interpol_mat` | Pattern B simplification; ledgered |
| `remap_data` | `R/utils_wavelet.R:233-262` | `remap_data` | Pattern B; ledgered |
| `gen_wavelet_indx` | `R/utils_wavelet.R:287-298` | `gen_wavelet_indx` | bit-identical for `lev_res >= 2`; ledgered |
| `dwt_matrix` | `R/utils_wavelet.R:331-365` | `DWT2` | bit-identical at default `max_scale = 10`; ledgered |
| `mf_dwt` | `R/dwt.R:51-120` | inline pipeline in EBmvFR / susiF_workhorse | bit-identical |
| `mf_invert_dwt` | `R/dwt.R:148-184` | `wd(...)` skeleton + `wr` | `1e-8` roundtrip; ledgered |
| `wd_variance` | `R/utils_wavelet.R:389-446` | `wd.var` | bit-identical (test asserts `tol = 0`) |
| `wst_variance` | `R/utils_wavelet.R:454-488` | `convert.var` | bit-identical |
| `av_basis_variance` | `R/utils_wavelet.R:497-520` | `AvBasis.var` (real branch) | bit-identical |
| `mf_invert_variance_curve` | `R/utils_wavelet.R:552-574` | n/a (Pattern A correction of upstream `(invert_dwt(sqrt(var_w)))^2`) | **No ledger entry** |
| `create_mf_individual` | `R/individual_data_class.R:37-189` | `multfsusie()` body and `EBmvFR()` setup | reorders quantile-norm and low-count-mask vs upstream |

### Findings

**B-1.1 (fix-now): `mf_invert_variance_curve` is a Pattern A bug fix without a ledger entry.** Add entry to `refactor-exceptions.md` covering `R/mfsusie_methods.R::.post_smooth_scalewise` and the variance-propagation derivation.

**B-1.2 (track-for-later): `create_mf_individual` reorders quantile-normalize and low-count-mask vs upstream.** Default-equivalent (only matters when both flags are nondefault). Add ledger entry explaining the rationale.

**B-1.3 (accept-as-documented): `mf_low_count_indices` returns `integer(0)` not `NULL`.** Cosmetic.

## 2. Prior init / SER / EM

| Routine | Verdict |
|---|---|
| `init_scale_mixture_prior_default` (`R/prior_scale_mixture.R:64-111`) vs `init_prior.default` (`R/operation_on_prior.R:42-185`) | bit-identical (test `1e-12`); ledgered |
| `mf_prior_scale_mixture` (`R/prior_scale_mixture.R:137-231`) | per-outcome wrapper; no algorithmic deviation |
| `mixture_log_bf_per_scale` (cpp11) vs `log_BF.mixture_normal_per_scale` | C++ port; test `1e-12` |
| `mixture_log_bf_per_scale_johnson` vs `log_BF` `df` branch | Student-t marginal; ledgered (`cor_small`) |
| `mixture_posterior_per_scale` (cpp11) vs `post_mat_mean` + `post_mat_sd` per scale | combined cpp11; test `1e-12` |
| `mf_em_likelihood_per_scale` vs `cal_L_mixsq_s_per_scale` | bit-identical structure |
| `mf_em_m_step_per_scale` vs `m_step.lik_mixture_normal_per_scale` -> `scale_m_step` | bit-identical mixsqp call |
| `univariate_smash_regression` vs upstream | test `1e-12` |

### Findings

**B-2.1 (accept-as-documented): mfsusieR omits the upstream `EM_pi` outer loop.** Consistent with the IBSS-loop delegation pattern (design.md D7).

**B-2.3 (accept-as-documented): `cal_lik` per-scale convergence not ported.** Legacy `EM_pi` `max_step_EM` is dead code at the upstream default of 1.

## 3. Post-smoothers

| Routine | fsusieR equivalent | Notes |
|---|---|---|
| `.post_smooth_scalewise` (`R/mfsusie_methods.R:457-501`) | n/a | mfsusieR-specific; uses Pattern A variance correction (B-1.1) |
| `.post_smooth_ti` -> `.univariate_ti_regression` (`R/mfsusie_methods.R:505-549, 591-673`) | `univariate_TI_regression` (`R/computational_functions.R:1756-1857`) | bit-identical kernel; test `tol = 0` |
| `.post_smooth_hmm` -> `mf_univariate_hmm_regression` -> `mf_fit_hmm` | `univariate_HMM_regression` -> `fit_hmm` | bit-identical post fsusieR PR #31; ledgered |
| `.post_smooth_smash` -> `univariate_smash_regression` | `smash_regression.susiF` -> `univariate_smash_regression` | kernel `1e-12`; outer-loop divergence (B-3.1) |

### Findings

**B-3.1 (track-for-later): `.post_smooth_smash` does not implement the upstream `n_iter = 5` outer Gauss-Seidel loop.** Single pass for L > 1; upstream alternates updates between effects. No L > 1 test exposes this.

**B-3.2 (track-for-later): `.post_smooth_ti` lacks the same `n_iter = 5` outer loop.** Same as B-3.1.

**B-3.3 (accept-as-documented): HMM cred bands return NULL by design.**

## 4. Covariate adjustment

| Routine | Verdict |
|---|---|
| `mf_adjust_for_covariates` dispatcher (`R/adjust_covariates.R:100-149`) | adds `"ols"` branch; upstream is wavelet-EB only |
| `mf_residualize_ols` | closed-form FWL; no upstream equivalent |
| `mf_residualize_wavelet_eb` | Pattern A bug fixes per ledger entries `EBmvFR-ER2`, `EBmvFR-MLE_wc2`; both load-bearing |

### Findings

**B-4.1 (accept-as-documented): EBmvFR Pattern A entries are still load-bearing.**

**B-4.2 (track-for-later): wavelet-EB ash refit uses R-level `ashr::*` calls instead of cpp11 `mixture_posterior_per_scale_cpp` kernel.** Cosmetic / performance.

**B-4.3 (fix-now): `mf_residualize_wavelet_eb` defaults disagree with the wrapper.** `wavelet_filter_number = 1L, wavelet_family = "DaubExPhase"` (inner) vs `wavelet_filter_number = 10L, wavelet_family = "DaubLeAsymm"` (wrapper). Wrapper passes its defaults through, so the inner defaults only matter on direct calls. Either align defaults or document.

## 5. Reference test fidelity

All fidelity tests against fsusieR pass at the documented tolerances:

- `test_fsusier_fidelity.R`: `1e-12` on per-modality log-BF, posterior moments; `0` on cal_zeta.
- `test_utils_wavelet.R`: `0` on values; `1e-12` on attribute deviations (Pattern B).
- `test_post_smooth_HMM.R`: `0` on bit-identity post fsusieR PR #31.
- `test_post_smooth_TI.R`: `0`.
- `test_post_smooth_smash.R`: `1e-12`.
- `test_post_smooth_wavelet_var.R`: `0` on the wd_variance / wst_variance / av_basis_variance pipeline.
- `test_post_smooth_scalewise.R`: `1e-14` on the corrected variance curve (closed-form).
- `test_prior_scale_mixture.R`: `1e-12` on prior init.
- `test_adjust_covariates.R`: `< 0.1` (Y_adjusted) and `< 0.05` (fitted_func) — bounded by the documented bug magnitude (Pattern A).

### Findings

**B-5.1 (track-for-later): `test_adjust_covariates.R:101-102` uses graduated tolerances forbidden by `refactor-discipline.md` section 3.** Justified inline by the Pattern A bug-magnitude framing, but the form clashes with the binary apple/orange philosophy. Rewrite to assert structural properties.

**B-5.2 (track-for-later): No top-level `fsusie(Y, X)` vs `fsusieR::susiF(Y, X)` end-to-end fidelity test confirmed.** Verify scope of `test_fsusier_degeneracy.R`.

## Final triage

| ID | Severity | Decision |
|---|---|---|
| B-1.1 | medium | fix-now (add ledger entry) |
| B-1.2 | low | track-for-later |
| B-1.3 | trivial | accept-as-documented |
| B-2.1 | n/a | accept-as-documented |
| B-2.3 | n/a | accept-as-documented |
| B-3.1 | medium | track-for-later (L > 1 path only) |
| B-3.2 | medium | track-for-later (L > 1 path only) |
| B-3.3 | n/a | accept-as-documented |
| B-4.1 | n/a | accept-as-documented |
| B-4.2 | low | track-for-later |
| B-4.3 | low | fix-now |
| B-5.1 | low | track-for-later |
| B-5.2 | medium | track-for-later |

**Two `fix-now` items. No undocumented numerical correctness divergences on the M = 1 path.**
