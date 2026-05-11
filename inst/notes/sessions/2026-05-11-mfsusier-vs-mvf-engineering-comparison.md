# Engineering comparison: mfsusieR HEAD vs mvf.susie.alpha

Date: 2026-05-11  
Branch: fix-mu-storage  
Companion MWE: `tests/testthat/test_mfsusier_vs_mvf_mwe.R`

Every section is derived from reading all source files in both packages.
Files read: `mfsusieR` — mfsusie.R, individual_data_class.R, individual_data_methods.R,
ibss_methods.R, dwt.R, mfsusie_methods.R, prior_scale_mixture.R, em_helpers.R,
save_mu_method.R, model_class.R, zzz.R, utils_wavelet.R, posterior_mixture.R,
post_smooth_hmm.R, post_smooth_smash.R.
`mvf.susie.alpha` — multfsusie.R, EM.R, computational_routine.R,
operation_on_multfsusie_obj.R, operation_on_multfsusie_prior.R,
ELBO_mutlfsusie.R, utils_wavelet_transform.R, utils.R.

---

## §1. Architecture

### mfsusieR

Entry point: `mfsusie()` (`R/mfsusie.R`). The function:
1. Calls `create_mf_individual()` to build the `mf_individual` data object
   (wavelet-domain `D[[m]]` + metadata).
2. Calls `mf_prior_scale_mixture()` to build `G_prior` (per-outcome mixture-normal prior).
3. Assembles `params` (a named list matching what `susie_workhorse` expects).
4. Calls `susieR::susie_workhorse(data, params)`, which runs the IBSS loop,
   dispatching every per-effect and per-iteration step to mfsusieR via S3 on
   `mf_individual`.
5. Attaches `dwt_meta` (inverse-DWT metadata) and smoothing inputs.
6. Applies `save_mu_method` storage trim.

S3 methods registered in `R/zzz.R` `.onLoad`:
- Per-effect: `compute_residuals`, `compute_ser_statistics`, `loglik`,
  `calculate_posterior_moments`, `optimize_prior_variance`,
  `update_fitted_values`, `compute_kl`, `SER_posterior_e_loglik`,
  `pre_loglik_prior_hook`, `post_loglik_prior_hook`.
- Per-iteration: `update_variance_components`, `update_derived_quantities`,
  `update_model_variance`, `Eloglik`, `get_objective`, `track_ibss_fit`,
  `trim_null_effects`, `ibss_initialize`, `cleanup_extra_fields`.
- Finalize: `get_scale_factors`, `get_intercept`, `get_fitted`,
  `get_variable_names`, `get_zscore`.

The IBSS loop itself runs entirely inside `susie_workhorse`; mfsusieR adds
no loop code.

### mvf.susie.alpha

Entry point: `multfsusie()` (`R/multfsusie.R`). The function runs its own
`while (check > tol)` IBSS loop (lines 488-590):

```
while (check > tol && iter < maxit) {
  for (l in 1:L) {
    R_l  <- cal_partial_resid(obj, l-1, X, Y, list_indx_lst)
    bs   <- cal_Bhat_Shat_multfsusie(R_l, X, obj, ind_analysis)
    EM   <- EM_pi_multsusie(bs, obj, G_prior, ...)  # mixsqp M-step
    obj  <- update_multfsusie(obj, l, EM, bs, ...)  # updates alpha, fitted_wc
  }
  obj   <- greedy_backfit(obj, ...)
  obj   <- estimate_residual_variance(obj, Y, X, ind_analysis)
  ELBO  <- get_objective(obj, Y, X, ind_analysis)
  check <- |ELBO[t] - ELBO[t-1]| / |ELBO[t]|   (if convergence_method="elbo")
}
```

No S3 delegation. All functions called directly.

---

## §2. Y input format

| Package | Expected format |
|---------|----------------|
| mfsusieR | `list(Y_1, ..., Y_M)` — flat list of M matrices, each `n x T_m` |
| mvf | `list(Y_f = list(Y_1, ..., Y_M), Y_u = matrix n x K_u)` — separates functional (`Y_f`) and univariate (`Y_u`) modalities |

mfsusieR supports scalar outcomes via `T_m = 1`: the DWT step is skipped and
the single column is treated as the wavelet-domain representation directly
(`individual_data_class.R:147`). Mixed ragged lists such as
`list(n x 1, n x 1, n x 64)` are valid inputs. There is no separate `Y_u`
channel; scalar outcomes share the same S3 code path as functional outcomes.

mvf has a two-channel design: `Y_f` (functional, DWT applied) and `Y_u`
(univariate, no DWT). The `Y_u` path calls a dedicated `m_step_u` and stores
results separately in `fitted_u` / `fitted_u2` / `sigma2$sd_u`. The two
channels combine at the log-BF level inside `EM_pi_multsusie`.

---

## §3. Wavelet preprocessing

### mfsusieR (`R/dwt.R`, `R/utils_wavelet.R`)

Per outcome `m`:
1. `remap_data(Y_m, pos_m)`: pads `T_m` to next power of 2.
2. `col_scale(Y_padded, center = TRUE, scale = TRUE)`: column-center + scale.
   `column_center` and `column_scale` stored in `dwt_meta` for inverse.
3. `dwt_matrix()`: row-by-row `wd()` from `wavethresh`; packs `[D | C]` into an
   `n x T_basis` matrix.
4. Optional `wavelet_standardize`: scales wavelet coefficients by their
   across-row standard deviation.
   When `wavelet_qnorm = TRUE` (default FALSE), applies rank-based quantile
   normalization inside the wavelet domain column by column (`mf_quantile_normalize`).
   NA bug present here: `rank(x, ties.method = "random")` without `na.last = "keep"`
   assigns NA entries large positive ranks instead of NA, deflating sigma2 by
   28-48% on real data. See §10.
5. NA rows: `na_idx[[m]] <- which(complete.cases(Y[[m]]))`. Rows with any NA
   are excluded from all residual computations for outcome `m`.

`D[[m]]` (wavelet coefficients) stored as `n x T_basis[m]` on `mf_individual`.
`scale_index[[m]]` maps each wavelet scale to its column indices in packed D.

### mvf (`R/utils_wavelet_transform.R`, `R/multfsusie.R`)

Per outcome `k`:
1. `DWT2(Y_k)`: row-by-row `wd()`, stores `D` (n x T-1) and `C` (length n).
   NA rows are imputed to 0 before DWT, then reinstated as NA after. Imputing
   0 changes wavelet coefficients in neighboring rows through filter support.
2. `pack_dwt(Y_k)`: `cbind(D, C)` — C in the last column.
3. Optional `quantile_trans = FALSE` (default): when TRUE applies quantile
   normalization in POSITION space BEFORE the wavelet transform. This is different
   from mfsusieR's `wavelet_qnorm` which operates in the wavelet domain.

---

## §4. Bhat / Shat (D2 divergence)

### mfsusieR (`R/individual_data_methods.R:mf_per_outcome_bhat_shat`)

```r
bhat_m[j, t] <- XtR_m[j, t] / xtx[j]
shat2_m[j, t] <- sigma2_per_pos[t] / xtx[j]
```

where `xtx[j] = sum_i X_{ij}^2` and `sigma2_per_pos` is derived from the
global per-outcome sigma2. The same sigma2 is shared across all variables.

### mvf (`R/computational_routine.R:cal_Bhat_Shat_multfsusie`)

Calls `fsusieR:::cal_Bhat_Shat(R_l, X, v1, ...)` — fixed, no parameter to
select a different formula. The function computes the marginal per-variable
effect estimate from univariate regression of each wavelet column:

```
Bhat[j, t] = (X_j^T R_t) / ||X_j||^2
Shat[j, t] = sqrt(Var(residual_t) / ||X_j||^2)   # per-(j, t) marginal SE
```

This uses a per-(j, t) residual variance rather than a global sigma2. At a true
causal SNP `j*`, `||X_{j*}||^2` is large and the per-variable marginal SE is
smaller than mfsusieR's global-sigma2 SE, giving mvf higher log-BF at signals.

This is NOT a bug in either package. mfsusieR's construction matches the SuSiE
variational derivation (posterior moments and likelihood use consistent sigma2);
mvf's construction follows fsusieR's marginal approach with potentially higher power
at signals but breaks the Gaussian-likelihood variational consistency property.

---

## §5. Prior structure

### mfsusieR (`R/prior_scale_mixture.R`)

Four `prior_variance_scope` modes:

| Mode | G_prior structure | M-step solver |
|------|------------------|--------------|
| `per_outcome` | One group per outcome, covers all T_basis columns | mixsqp (one solve per outcome per effect) |
| `per_scale` | One group per (outcome, wavelet scale) | mixsqp (S_m solves per outcome per effect) |
| `per_scale_normal` | One ebnm_point_normal per (outcome, scale) | `ebnm::ebnm_point_normal` |
| `per_scale_laplace` | One ebnm_point_laplace per (outcome, scale) | `ebnm::ebnm_point_laplace` |

Grid init (mixsqp paths): marginal Bhat/Shat via `compute_marginal_bhat_shat`;
`sd_min = quantile(Shat, 0.1) / 10`, `sd_max = 2 * sqrt(max(Bhat^2 - Shat^2))`,
`gridmult = sqrt(2)` (following Stephens 2017 but with 10th-percentile lower end
to avoid degenerate tiny-grid components).

### mvf (`R/operation_on_multfsusie_prior.R`)

`init_prior_multfsusie()` calls:
- For functional modalities: `fsusieR:::init_prior.default()` which calls
  `ashr::ash()` on marginal Bhat/Shat per outcome per scale, returning a
  `mixture_normal_per_scale` ash object.
- For univariate modalities: `ashr::ash()` per column of `Y_u`.

---

## §6. M-step (mixture-weight update)

### mfsusieR (`R/em_helpers.R`, `R/individual_data_methods.R`)

`mf_em_likelihood_per_scale()` builds the mixsqp `L` matrix from the `(bhat, shat)`
slice for a (outcome, scale, keep_idx) rectangle, plus a `(100, 0, ..., 0)` penalty row.

`mf_em_m_step_per_scale()` calls `mixsqp`:
- Weight vector: `w = c(mixture_null_weight * idx_size, rep(zeta_keep, idx_size))`.
- `mixture_null_weight` scaled by M inside `.opv_mixsqp`: `mnw * max(1, M)`.
- Warm-started from `model$G_prior[[m]][[s]]$fitted_g$pi` (prior iter pi).
  `control_mixsqp = NULL` (default) uses warm-start with fast defaults.
  Cold-start equivalent: `list(convtol.sqp = 1e-8, numiter.em = 20, tol.svd = 1e-10)`.

`max_SNP_EM` does NOT exist in mfsusieR; all p variables enter the M-step.

### mvf (`R/EM.R`)

`EM_pi_multsusie()`: for functional modalities, calls `fsusieR::m_step` per scale;
for univariate, calls `m_step_u`. Both use cold start every outer iteration:
```r
x0 <- c(init_pi0_w = 0.9, rep(1e-12, K - 1))
```

`nullweight` in mvf is divided by K (number of modalities) inside `EM_pi_multsusie`:
`nullweight_m = nullweight / K`.

`max_SNP_EM = 100` (default): before the M-step, limits the input to the top-100
SNPs by log-BF. This is a computational shortcut with no mfsusieR equivalent; it
means the mixture weights are learned from at most 100 data points per (effect, scale).

---

## §7. Inner EM loop (mfsusieR only)

`post_loglik_prior_hook.mf_individual` (`R/individual_data_methods.R:672`):

```
inner_cap = max(0, max_inner_em_steps) + 1   # default max_inner_em_steps=5 -> cap=6
for k in 1..inner_cap:
    M-step: opv_fn(...)     # updates G_prior + fitted_g_per_effect[[l]]
    if k == inner_cap: break
    re-run loglik / moments / KL to update alpha[l, ] against new pi
    if |lbf[l] - lbf_prev| < inner_tol: break
```

After `inner_cap` M-step cycles, the per-effect `(alpha[l, ], pi_V[[l]])` end
in lockstep within one outer IBSS iteration.

`max_inner_em_steps = 0` (one M-step, no re-evaluation) approximates mvf's
one-step behavior. `max_inner_em_steps = 5` (default) runs up to 6 cycles.

The inner loop breaks the IBSS ELBO monotone guarantee for the outer loop.
A coherent ELBO is restored at iter end via `refresh_lbf_kl.mf_individual`
(inside `get_objective.mfsusie`).

mvf has no inner EM loop; it runs one M-step per effect per outer iteration.

---

## §8. Residual variance update (D1 divergence)

### mfsusieR (`R/individual_data_methods.R:400-421`)

`mf_get_ER2_per_position(data, model, m)`:
```r
res_m  <- (D[[m]] - fitted[[m]])[na_idx[[m]], ]
rss_t  <- colSums(res_m^2)
bias_t <- sum over l of:
    colSums((alpha_l * pw) * mu2_l_m)     # pw = xtx_diag_list[[m]], length p
    - colSums((X %*% (alpha_l * mu_lm))^2)
sigma2_m <- sum(rss_t + bias_t) / (n_m * T_basis[m])   # per_outcome
```

The bias correction uses `predictor_weights = xtx_diag` in `colSums((alpha_l * pw) * mu2)`.
This is the correct derivation.

### mvf (`R/operation_on_multfsusie_obj.R:get_ER2.multfsusie`)

```r
ER2$f[k] <- sum((Y_f[[k]] - X %*% postF$post_f[[k]])^2)
             - sum(postF$post_f[[k]]^2)
             + sum(postF2$post_f_sd2[[k]])
```

The second term uses `sum(E[b_l]^2)` (without `xtx_diag`) rather than the correct
`sum(xtx_diag * E[b_l^2]) - sum((X * E[b_l])^2)`. The missing `xtx_diag` factor
(O(n) on average) deflates the bias correction term, making ER2 smaller and sigma2
smaller. This is the D1 bug.

Practical effect: mfsusieR sigma2 is consistently larger than mvf sigma2 at the same
data. The ratio is empirically 1.5-2x in perm runs. Under-estimated sigma2 in mvf
makes lbf values larger, inflating CS counts and FDR.

---

## §9. KL / ELBO (D3 divergence)

### mfsusieR

`compute_kl.mf_individual`: inherits from susieR's correct formula.
`get_objective.mfsusie` calls `refresh_lbf_kl.mf_individual` at iter end to
re-evaluate lbf/KL against the iter-final pi, returning a coherent ELBO.

### mvf (`R/ELBO_mutlfsusie.R:55`)

`cal_KL_l.multfsusie`:
```r
out <- -loglik_SFR(obj, l, Y, X, ind_analysis)
       - loglik_SFR_post(obj, l, R_l, X, ind_analysis)
```

The correct formula is `KL = loglik_SFR_post - loglik_SFR` (posterior expected
log-likelihood minus prior marginal log-likelihood). mvf negates both terms.
Effect: `sum(KL)` in the ELBO is `-(sum_loglik_SFR + sum_loglik_SFR_post)`,
an O(L * n * T) inflation of the reported objective. Consecutive ELBO differences
still shrink as alpha stabilizes, so ELBO-convergence termination is unaffected.
PIP-convergence runs (default in both packages) are completely unaffected.

---

## §10. NA handling and qnorm bug (mfsusieR)

`mf_quantile_normalize` (`R/utils_wavelet.R`):
```r
rank(x, ties.method = "random")   # missing: na.last = "keep"
```
NA entries receive large positive ranks instead of NA. In real atac-seq data with
~5-10 samples missing for some cell types, this deflated sigma2 by 28-48%.

Fix: `rank(x, na.last = "keep", ties.method = "random")`. This bug is in the
`wavelet_qnorm = TRUE` path (default FALSE), so it only affects users who enable
wavelet quantile normalization.

mvf's `DWT2` imputes NA rows to 0 before DWT (a different trade-off that changes
neighboring wavelet coefficients via filter support). Neither solution is perfect;
the mfsusieR fix above is correct for the rank-normalization path.

---

## §11. LFSR

### mfsusieR

`lfsr_from_gaussian(mean, sd) = pnorm(-|mean| / sd)` (`R/mfsusie_methods.R:460`).
Computed during `mf_post_smooth()` on the smoothed per-effect curve in position space.
Stored at `fit$smoothed[[method]]$lfsr_curves[[m]][[l]]` (a list of length-T_m
numeric vectors, one per (outcome, effect) pair).

A per-variant clfsr at the wavelet-coefficient level is also stored at
`clfsr_curves[[m]][[l]]` (p x T_basis matrix) when `save_mu_method = "complete"`,
derived from `(mu[[l]][[m]], mu2[[l]][[m]])` via the Gaussian approximation.

All four post-smooth methods (TI, HMM, smash, scalewise) produce lfsr.

### mvf

`HMM_regression.multfsusie` (called inside `out_prep` when `post_processing = "HMM"`)
delegates to `fsusieR::HMM_regression.susiF` per functional outcome. HMM-derived
lfsr is stored at `fit$lfsr[[l]]$est_lfsr_functional[[k]]`.

For `post_processing` choices other than `"HMM"` (smash, TI, none), mvf does NOT
compute lfsr. The lfsr output is thus tightly coupled to the post-processing
choice baked into `out_prep`.

---

## §12. Posthoc trait configuration probabilities

### mfsusieR

Implemented via `susieR::susie_post_outcome_configuration()`. mfsusieR stores a
per-(effect, variant, outcome) log Bayes factor array `lbf_variable_outcome`
(L x p x M) during the IBSS sweep (`R/individual_data_methods.R:231`). Users call
this susieR function after fitting:

```r
susieR::susie_post_outcome_configuration(fit, by = "outcome", method = "susiex")
```

Two modes via `by`:
- `by = "outcome"` (for mfsusie): expands the single fit into M per-outcome views
  using `lbf_variable_outcome[, , m]` slices. Each view has `alpha` (the shared
  L x p matrix) and per-outcome `lbf` (L x p). Runs `susiex_configurations()` on
  the M views.
- `by = "fit"`: treats the whole fit as one trait, using the joint composite
  `lbf_variable` (L x p). Used for pairwise comparison across separate fits.

Two computation methods via `method`:
- `"susiex"`: enumerates `2^M` configurations per CS tuple. For each tuple of
  L-indices (one per outcome), computes
  `logBF_trait[m] = sum_j alpha[l_m, j] * lbf[l_m, j, m]`, then
  `prob_conf = normalize(exp(configs %*% logBF_trait))`,
  `marginal_prob = crossprod(configs, prob_conf)`.
  Same enumeration logic as mvf's `posthoc_multfsusie`.
- `"coloc"`: pairwise Bayes-factor coloc (`coloc.abf` style) for N=2 comparisons
  across separate fits. Not applicable to a single mfsusie fit with `by = "outcome"`.

Output: a `susie_post_outcome_configuration` object with either `$susiex` (list of
CS tuples, each with `cs_indices`, `logBF_trait`, `configs`, `config_prob`,
`marginal_prob`, `active`) or `$coloc_pairwise`.

NOT called automatically; the user must invoke it after `mfsusie()`.

### mvf

`posthoc_multfsusie()` (`R/operation_on_multfsusie_obj.R:1514`), called inside
`out_prep` when `posthoc = TRUE` (default). Per credible set `l`:
1. Compute per-trait log BF: `logBF_l[k] = get_cs_logBF_multfsusie(alpha_l, lBF_per_trait_l)`.
2. Enumerate `2^S` trait configurations (S = number of modalities; capped at S <= 20).
3. `prob_conf = normalize(exp(configs %*% logBF_l))`.
4. Marginal per-trait probability: `posthoc_trait = colSums(configs * prob_conf)`.

Output at `fit$posthoc[[l]]`: list with `logBF_trait`, `posthoc`, `active`
(posthoc >= 0.8), `configs`, `config_prob`. Called automatically by default.

Reference: Yuan et al., Nat Genet 2024.

### Key differences

| Aspect | mfsusieR | mvf |
|--------|---------|-----|
| Where implemented | `susieR::susie_post_outcome_configuration()` | `posthoc_multfsusie()` in mvf |
| When called | User calls explicitly after `mfsusie()` | Automatic inside `out_prep` |
| Storage | `fit$lbf_variable_outcome` (L x p x M) enables the call | `fit$posthoc[[l]]` stores results |
| Methods | `"susiex"` (2^M enum) and `"coloc"` (pairwise BF) | susiex-style only |
| Input flexibility | Accepts single fit (`by="outcome"`) or list of fits (`by="fit"`) | Per-CS on one fit only |

---

## §13. Output structure

### mfsusieR

| Field | Type | Content |
|-------|------|---------|
| `fit$alpha` | L x p **matrix** | variational SNP-level PIP per effect |
| `fit$mu[[l]][[m]]` | p x T_basis[m] matrix (or 1 x T if trimmed) | wavelet posterior mean |
| `fit$mu2[[l]][[m]]` | p x T_basis[m] matrix (or 1 x T if trimmed) | wavelet posterior second moment |
| `fit$sets$cs` | list of CS objects (susieR format), each with `$variables` and `$purity` | credible sets |
| `fit$pip` | length-p vector | aggregate PIP across L effects |
| `fit$sigma2` | list[M] of scalar or length-S_m vector | per-outcome residual variance (not SD) |
| `fit$pi_V[[l]][[m]]` | S_m x K matrix | mixture weights per (effect, outcome, scale) |
| `fit$dwt_meta` | list | inverse-DWT parameters |
| `fit$smoothed[[method]]` | list | populated after `mf_post_smooth()` |
| `fit$lbf_variable_outcome` | L x p x M array | per-(effect, variant, outcome) lbf; input to `susie_post_outcome_configuration()` |
| `fit$posthoc` | NULL | not stored; call `susieR::susie_post_outcome_configuration(fit)` explicitly |

### mvf

| Field | Type | Content |
|-------|------|---------|
| `fit$alpha` | **list of L plain numeric vectors**, each length p | NOT an L x p matrix |
| `fit$fitted_wc[[l]][[k]]` | p x T_basis[k] matrix | wavelet posterior mean per (effect, functional outcome) |
| `fit$fitted_wc2[[l]][[k]]` | p x T_basis[k] matrix | wavelet posterior second moment |
| `fit$cs` | list of L integer vectors (SNP indices) | flat list, NOT susieR format |
| `fit$purity` | matrix (one row per CS) | computed top-level by `fsusieR::cal_purity` |
| `fit$pip` | length-p vector | aggregate PIP |
| `fit$sigma2$sd_f` | length-M vector of **SDs (not variances)** | residual SD per outcome |
| `fit$fitted_func[[l]][[k]]` | length-T_basis vector | position-space curve; baked in by out_prep |
| `fit$cred_band[[l]][[k]]` | 2 x T_basis matrix | credible band; baked in |
| `fit$lfsr[[l]]$est_lfsr_functional[[k]]` | length-T_basis vector | only when post_processing = "HMM" |
| `fit$posthoc[[l]]` | list | config probabilities |

Key access differences for MWE code:
- Alpha matrix: mfsusieR `fit$alpha[l, ]` vs mvf `fit$alpha[[l]]`
- Convert mvf alpha to matrix: `do.call(rbind, fit_mvf$alpha)` (NOT `lapply(x, function(a) a$alpha_f)`)
- Sigma2: mfsusieR `unlist(fit$sigma2)` = variances; mvf `unlist(fit$sigma2$sd_f)^2` = variances
- CS purity: mfsusieR `cs[[i]]$purity` (nested); mvf `fit$purity` (top-level)

---

## §14. save_mu_method (mfsusieR only)

Three modes via `mfsusie(save_mu_method = ...)`:
- `"complete"` (default): full p x T_basis per (l, m). Supports warm-start
  (`model_init`), per-variant lfsr, `predict(newx)`.
- `"alpha_collapsed"`: replaces p x T by 1 x T alpha-weighted summary.
  A separate `coef_wavelet[[l]][[m]]` (1 x T) is precomputed before collapse
  (per-j csd_X scaling cannot be recovered post-collapse). Factor-p storage reduction.
- `"lead"`: keeps only the lead variable `j* = which.max(alpha[l, ])`;
  `top_index[l]` records `j*`. Cheapest but biased to lead.

`mf_thin(fit, method)` applies the trim post-fit without modifying the original.

mvf always stores full p x T in `fitted_wc`. No equivalent storage mode.

---

## §15. Warm start / model_init

mfsusieR: `mfsusie(..., model_init = prior_fit)` warm-starts the IBSS loop from
a prior fit's `fitted_g_per_effect`, `alpha`, `mu`, and `fitted` fields.
Requires `save_mu_method = "complete"` on the prior fit.

mvf: always cold start. No equivalent parameter.

---

## §16. Greedy / backfit

### mfsusieR

`L_greedy = NULL` (default): no greedy expansion. `L` effects are allocated at
init; `trim_null_effects` prunes effects whose effective slab variance falls
below `prior_tol` (susieR default 1e-9). Number of effects is fixed.

### mvf

`greedy = TRUE` (default): starts with `L_start` effects, adds blocks of 7 via
`expand_multfsusie_obj` when all effects have non-dummy CS. `backfit = TRUE`
(default): prunes low-purity effects after each outer iter via `discard_cs`.
Both are on by default, making the loop structure substantially different from
mfsusieR's fixed-L approach.

To achieve a comparable fixed-L run: `greedy = FALSE, backfit = FALSE`.

---

## §17. Convergence

| Parameter | mfsusieR default | mvf default |
|-----------|-----------------|-------------|
| Criterion | `convergence_method = "pip"` | ELBO-based (`check = |ELBO_diff| / |ELBO|`) |
| Tolerance | `tol = 1e-3` | `tol = 1e-3` |
| Max iter | `max_iter = 100` | `maxit = 100` |

mfsusieR's pip convergence: `max(|pip[t] - pip[t-1]|) < tol`. Unaffected by the
ELBO sign error in mvf (D3).

---

## §18. Nullweight parameterization

Both packages allow setting the null-component penalty for the mixture M-step.

mfsusieR: `mixture_null_weight` (default NULL, resolves to 0.05 internally).
Passed as the direct `mixture_null_weight` to `mf_em_m_step_per_scale`, then scaled
by M inside `.opv_mixsqp`: effective weight = `mixture_null_weight * max(1, M)`.

mvf: `nullweight = 0.7` (default). Divided by K (number of modalities) inside
`EM_pi_multsusie`: `nullweight_m = nullweight / K`.

For M = 3 (three functional outcomes), comparable settings are approximately:
mfsusieR `mixture_null_weight = 0.7 / 3 ≈ 0.23` vs mvf `nullweight = 0.7`
(which applies `0.7 / 3 ≈ 0.23` per outcome internally). The scaling strategy
is equivalent in intent; the default values differ.

---

## §19. Confirmed divergences

| ID | Component | mfsusieR | mvf | Impact |
|----|-----------|----------|-----|--------|
| D1 | ER2 bias correction | Correct (`xtx_diag` present in `colSums((alpha * pw) * mu2)`) | Missing `xtx_diag` factor; sigma2 deflated ~1.5-2x | FDR inflation in all real-data runs |
| D2 | Bhat/Shat | Global `sigma2/xtx` per outcome | Per-(j,t) marginal SE from univariate regression | Detection power difference at true signals; not a bug |
| D3 | KL/ELBO sign | Correct | `-loglik_SFR - loglik_SFR_post` instead of `loglik_SFR_post - loglik_SFR`; ELBO inflated by O(L*n*T) | ELBO value wrong; PIP and CS outputs unaffected |
| D4 | qnorm NA | `rank(x)` without `na.last = "keep"` deflates sigma2 28-48% on real data | DWT2 imputes NA to 0 (different trade-off) | Real-data FDR inflation when `wavelet_qnorm = TRUE` |

D1 and D4 are bugs. D2 is an intentional design choice in each package's variational
derivation. D3 is a derivation error in mvf that does not affect PIP convergence runs.

---

## §20. Parameter name mapping

| Concept | mfsusieR | mvf |
|---------|---------|-----|
| Number of effects | `L` | `L` |
| Null weight | `mixture_null_weight` (default NULL -> 0.05; scaled by M) | `nullweight` (default 0.7; divided by K) |
| Max iterations | `max_iter` | `maxit` |
| Convergence tol | `tol` | `tol` |
| Prior mode | `prior_variance_scope = "per_outcome"/"per_scale"/...` | `prior = "mixture_normal"/"mixture_normal_per_scale"` |
| Inner EM steps | `max_inner_em_steps = 5L` | no equivalent |
| Mixsqp control | `control_mixsqp = NULL` (warm) / `list(...)` (cold) | no exposed control |
| Warm start | `model_init = prior_fit` | no equivalent |
| Wave qnorm | `wavelet_qnorm = FALSE` (in wavelet domain) | `quantile_trans = FALSE` (in position domain, before DWT) |
| Small-sample BF | `small_sample_correction = FALSE` | `cor_small = FALSE` |
| Post-processing | `mf_post_smooth(fit, method)` (separate step after fit) | `post_processing = "smash"/"TI"/"HMM"/"none"` (baked into return) |
| Posthoc | not available | `posthoc = TRUE` (default) |
| Storage trim | `save_mu_method = "complete"/"alpha_collapsed"/"lead"` | not available |
| Greedy | `L_greedy = NULL` (disabled by default) | `greedy = TRUE` (default) |
| Backfit | implicit via `trim_null_effects` | `backfit = TRUE` (default) |
| Max SNPs for M-step | all p (no limit) | `max_SNP_EM = 100` (top-100 by log-BF) |

---

## §21. Next session

- Fix the NA bug in `mf_quantile_normalize` (`R/utils_wavelet.R`), pending user
  approval. Change: add `na.last = "keep"` to `rank()` call.
- Monitor SLURM job 32644042 (20 tasks, conditions A-D); aggregate results when
  all 20 CSVs appear.
- Verify the MWE tests run without error when `mvf.susie.alpha` is installed.
