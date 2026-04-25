# mvf.susie.alpha: the port source

This note maps `mvf.susie.alpha`, William Denault's original implementation of
multi-functional SuSiE. It is the source from which mfsusieR is being ported.
Anjing has reported FDR miscalibration (observed FDR exceeds nominal); part of
the purpose of this note is to identify candidate sites for that bug.

Source root: `~/GIT/mvf.susie.alpha`. All citations are relative to that root.

## 1. Public API

Entry function `multfsusie()` (`R/multfsusie.R:172-208`):

```
multfsusie(Y, X, L=2, pos=NULL, prior="mixture_normal",
           post_processing=c("smash","TI","HMM","none"), verbose=TRUE,
           maxit=100, tol=1e-6, cov_lev=0.95, min_purity=0.5, L_start=3,
           filter_cs=TRUE, init_pi0_w=1, nullweight=.7, max_scale=10,
           max_SNP_EM=100, gridmult=sqrt(2), thresh_lowcount, cal_obj=FALSE,
           greedy=TRUE, backfit=TRUE, max_step_EM=1, cor_small=FALSE,
           filter.number=10, family="DaubLeAsymm", e=0.001,
           tol_null_prior=0.001, lbf_min=0.1, posthoc=TRUE,
           max_step=100, multfsusie.obj=NULL)
```

`Y` is a list with two slots: `Y$Y_f` (a list of functional responses, one
matrix per trait, each n by T_k) and `Y$Y_u` (a matrix of univariate
responses, n by n_u). The package exports ~70 functions including the S3
methods on the `multfsusie` class (`update_alpha.multfsusie`,
`update_ELBO.multfsusie`, `estimate_residual_variance.multfsusie`,
`update_cal_pip.multfsusie`, `update_cal_cs.multfsusie`,
`check_cs.multfsusie`, `discard_cs.multfsusie`, ...).

Typical workflow: call `multfsusie(Y, X, pos, L=...)`, inspect `fit$cs`,
`fit$pip`, `fit$fitted_func`, `fit$lfsr`. Internally the data are remapped to
a 2^J grid, wavelet-decomposed row-by-row, run through IBSS with EM-fitted
mixture priors, and post-processed by inverse DWT plus an optional smoother
(SMASH, TI thresholding, or HMM).

## 2. Core data structures

The fit object (class `multfsusie`, initialized at
`R/operation_on_multfsusie_obj.R:508-644`):

- `alpha` - list of length L, each a length-p vector of posterior inclusion
  probabilities for one effect.
- `pip` - length-p vector, computed as `1 - prod_l (1 - alpha[[l]])`
  (`R/operation_on_multfsusie_obj.R:2001`).
- `cs` - list of length L of credible sets (variant indices).
- `fitted_wc` - list of L; for each effect, a list of K matrices p by n_wac[k]
  holding fitted wavelet coefficients per functional trait.
- `fitted_u` - list of L; for each effect, a p by n_u matrix of univariate
  fitted effects.
- `fitted_func` - list of L of reconstructed functional curves after inverse
  DWT and optional smoothing.
- `lBF` - list of L of length-p log Bayes factor vectors.
- `sigma2` - list with `sd_f` (per-trait functional residual SD) and `sd_u`
  (per-trait univariate residual SD), initialized by `init_var_multf()`
  (`R/operation_on_multfsusie_obj.R:606`), updated each iteration
  (`R/multfsusie.R:481-485, 562-567`). Note the misleading name: stores SD,
  not variance.
- `G_prior` - object of class `multfsusie_prior` containing `G_prior_u` (list
  of `ash` priors, one per univariate trait) and `G_prior_f` (list of
  `mixture_normal_per_scale` priors from `fsusieR`, one per functional
  trait).
- `n_wac` - list of K of wavelet-coefficient counts per trait (= 2^J after
  padding).
- `csd_X` - vector used to unstandardize X.

Wavelet representation:

- Each functional trait is padded to length 2^J via
  `fsusieR::remap_data()` (`R/multfsusie.R:284`).
- Columns are scaled with `fsusieR::colScale()` (`R/multfsusie.R:308`).
- Each row goes through `DWT2()` (`R/utils_wavelet_transform.R:24-50`), a
  thin wrapper around `wavethresh::wd()`. Daubechies least-asymmetric basis
  with 10 vanishing moments by default.
- Coefficients are stored as `cbind(D, C)` (`R/multfsusie.R:311`); the index
  list `list_indx_lst[[k]]` (from `fsusieR::gen_wavelet_indx()`) maps columns
  to scales.

`pos` is either NULL (assumed equally spaced, `pos[[i]] <- 1:ncol(Y$Y_f[[i]])`)
or user-provided; in either case it is consumed by `fsusieR::remap_data()`.

## 3. Wavelet pathway

Forward transform (`R/multfsusie.R:308-316`):

1. `fsusieR::remap_data()` pads to 2^J (`R/multfsusie.R:284`).
2. `fsusieR::colScale()` standardizes columns (`R/multfsusie.R:308`).
3. `DWT2()` applies `wavethresh::wd()` row by row
   (`R/utils_wavelet_transform.R:24-50`). Output has detail coefficients D
   ((J-1) by n) and the constant C (length n).
4. Stacked as `cbind(temp$D, temp$C)`.
5. `list_indx_lst[[k]] <- fsusieR::gen_wavelet_indx(...)` records which
   coefficient columns belong to which scale.

Scale-specific prior variance is inside `fsusieR::init_prior.default()`
(`R/operation_on_multfsusie_prior.R:81-96`). The functional prior is a
`mixture_normal_per_scale` from `fsusieR`; the per-scale variance comes
out of an empirical-Bayes fit `fsusieR` performs internally. There is no
explicit user-facing `max_scale` filtering of wavelet coefficients;
`max_scale` only controls the padded length used during data remap
(`R/multfsusie.R:287`).

Inverse transform (`R/operation_on_multfsusie_obj.R:2074-2099`,
`update_cal_fit_func()`): for the "none" post-processing path it constructs
a dummy wavelet object via `wavethresh::wd(rep(0, n_wac[[k]]))`, sets the
detail coefficients to `alpha[[l]] %*% fitted_wc[[l]][[k]][,-last_scale_idx]`
and the constant to the corresponding column, then calls `wavethresh::wr()`
to reconstruct the curve. Scaling is reversed by dividing by `csd_X`
(`R/operation_on_multfsusie_obj.R:2091-2092`).

For the smoother paths (`smash`, `TI`, `HMM`) the package delegates to
fsusieR's reconstruction routines (`R/operation_on_multfsusie_obj.R:2301-2477`).

Boundary handling is whatever `wavethresh` provides by default; padding is
done by `fsusieR::remap_data()`.

## 4. IBSS loop structure

Outer loop at `R/multfsusie.R:489-589`:

```
while (check > tol & iter <= maxit & n_step <= max_step) {
  for (l in 1:L) {
    update_Y       <- cal_partial_resid(...)        # line 495
    Bhat_Shat      <- cal_Bhat_Shat_multfsusie(...)  # line 519
    EM_out         <- EM_pi_multsusie(...)           # line 527
    multfsusie.obj <- update_multfsusie(...)         # line 541
  }
  multfsusie.obj <- greedy_backfit(...)              # line 555
  sigma2         <- estimate_residual_variance(...)  # lines 562-567
  multfsusie.obj <- update_residual_variance(...)
  multfsusie.obj <- test_stop_cond(...)              # line 571
}
```

Per-iteration order:

1. For l = 1..L: partial residuals minus all other effects, marginal
   sufficient statistics `Bhat`/`Shat`, EM update of mixture proportions,
   posterior moment update (`fitted_wc[[l]]`, `fitted_u[[l]]`, `alpha[[l]]`).
2. `greedy_backfit()` to add or refit effects
   (`R/operation_on_multfsusie_obj.R:1111-...`).
3. Residual variance update (per-trait scalar, see section 6).
4. ELBO check via `test_stop_cond()`.

The per-effect SER analogue:

- `EM_pi_multsusie` (`R/EM.R:31-120`) calls `log_BF`
  (`R/computational_routine.R:259-333`) which dispatches to
  `fsusieR::log_BF` per functional trait (`R/computational_routine.R:283`)
  and to `log_BFu` for univariate traits (`R/computational_routine.R:201`).
- alpha is computed from the aggregated lBF via `fsusieR::cal_zeta()`
  (`R/operation_on_multfsusie_obj.R:1739`).
- Posterior moments are aggregated across scales inside fsusieR.

Convergence: if `cal_obj = TRUE` ELBO is computed via `get_objective()` and
`test_stop_cond()` checks `|ELBO[iter] - ELBO[iter-1]| < tol` with default
`tol = 1e-6` (`R/multfsusie.R:178`).

## 5. Priors and posteriors

Prior object `multfsusie_prior` (`R/operation_on_multfsusie_prior.R:17-138`):

- `G_prior_u`: per-univariate-trait `ash` priors built via
  `ashr::ash(Bhat, Shat, mixcompdist="normal", gridmult=gridmult, ...)`
  (`R/operation_on_multfsusie_prior.R:40`).
- `G_prior_f`: per-functional-trait `mixture_normal_per_scale` priors from
  `fsusieR::init_prior.default(..., indx_lst=..., ...)`
  (`R/operation_on_multfsusie_prior.R:82`).
- Initialized once via `init_prior_multfsusie()` (`R/multfsusie.R:391`).

Posterior moments (`R/operation_on_multfsusie_obj.R:1705-1753`,
`get_post_effect_multfsusie()`): per effect, posterior mean and second moment
are computed and stored in `fitted_wc[[l]][[k]]` / `fitted_u[[l]]` and
`fitted_wc2[[l]][[k]]` / `fitted_u2[[l]]`.

The mixture prior over scales lives inside fsusieR; per scale, an empirical
Bayes mixture-of-normals prior is fit. lBF is aggregated across scales there
before being summed across traits in `EM.R:48-50`.

## 6. Residual variance and noise model

`estimate_residual_variance.multfsusie` (`R/computational_routine.R:395-431`):

```
R2 <- get_ER2(multfsusie.obj, Y, X, ind_analysis)   # expected squared residuals
est_sd2$sd_u <- R2$uni / nrow(Y$Y_u)                # line 402
est_sd2$sd_f <- R2$f   / (n * t)                    # line 407
```

Where `n = nrow(Y$Y_f[[k]])`, `t = ncol(Y$Y_f[[k]])`. Note again that the
field name is `sd` but the value is in fact a variance (R2 over n).

`get_ER2()` (`R/operation_on_multfsusie_obj.R:1046-1096`) computes
`sum((Y - X B)^2) - sum(B_hat^2) + sum(B2_hat)`, i.e. the expected squared
residual using posterior first and second moments.

The residual variance is *global per trait*. Even though the prior variance
is per scale (via fsusieR's `mixture_normal_per_scale`), a single
`sigma2$sd_f[k]` is shared across all wavelet scales of trait k. Cross-position
or cross-scale noise correlation is not modeled. This is a structural mismatch
that the FDR investigation in Phase 5 should pay attention to.

## 7. Credible sets and PIP computation

PIP at `R/operation_on_multfsusie_obj.R:1990-2003`:

```
update_cal_pip.multfsusie <- function(multfsusie.obj) {
  tpip <- list()
  for (l in 1:L) tpip[[l]] <- 1 - alpha[[l]]
  pip <- 1 - apply(do.call(rbind, tpip), 2, prod)   # line 2001
}
```

Credible sets (`update_cal_cs.multfsusie`, called from `out_prep`,
`R/operation_on_multfsusie_obj.R:1450-1454`):

- Sort alpha[[l]] descending, accumulate until cumulative probability
  reaches `cov_lev` (default 0.95).
- Apply purity threshold `min_purity` (default 0.5).
- For functional traits, scale-aware CS construction is delegated to fsusieR.

The post-processing path in `out_prep`:

1. `update_cal_pip.multfsusie()` - computes `pip` from current `alpha`.
2. If `filter_cs = TRUE`: `check_cs()` removes low-purity effects,
   `merge_effect()` merges low-signal effects, `discard_cs()` flags effects
   for removal (`R/operation_on_multfsusie_obj.R:1452-1454`).
3. `discard_cs()` removes the flagged effects from `alpha` and friends.

`update_cal_pip()` is *not* recalled after step 3. The `pip` vector therefore
includes contributions from effects that have been discarded. See the bug
candidates below.

## 8. Candidates for the FDR bug

Cited file:line and the suspect mechanism. Confidence flagged.

1. *PIP not refreshed after CS filtering* (high confidence).
   `update_cal_pip()` runs at `R/operation_on_multfsusie_obj.R:1447`, but
   `check_cs` / `merge_effect` / `discard_cs` run at `:1452-1454` and remove
   effects from `alpha` without recomputing `pip`. A user thresholding on
   `fit$pip` at e.g. 0.05 sees probabilities that include discarded
   contributions. Direct route to inflated discovery counts.

2. *Residual variance estimator divides by n only* (medium-high).
   `R/computational_routine.R:402, :407`: `est_sd2$sd_u <- R2$uni / nrow(...)`
   and `est_sd2$sd_f <- R2$f / (n * t)`. No degrees-of-freedom correction for
   the L estimated effects. Residual variance is biased low; lBFs are biased
   high; false discoveries grow with L. Standard variance-bias issue.

3. *Global residual variance per trait, scale-specific prior variance*
   (low-medium). Section 6 above. The mismatch between scale-specific signal
   modeling and trait-level noise modeling can cause systematic miscalibration
   when true signal is concentrated in a small number of scales.

4. *Empirical-Bayes prior variance with no df correction* (low-medium).
   `R/operation_on_multfsusie_prior.R:81-96` and the fsusieR side. ash
   estimates mixture proportions via empirical Bayes, which can underestimate
   the null component when sample size is large; the consequence is inflated
   non-null prior variance and Bayes factors.

5. *lBF aggregation across modalities without weighting* (medium). `R/EM.R:48-50`:
   ```
   f_logBF <- apply(lBF_per_trait$f_logBF, 2, sum)
   u_logBF <- apply(lBF_per_trait$u_logBF, 2, sum)
   lBF     <- f_logBF + u_logBF
   ```
   The functional and univariate log-BFs are summed with no per-modality
   variance weighting. A single high-signal univariate trait can dominate
   the EM and inflate inclusion probabilities for other modalities.

6. *ELBO convergence threshold may stop early* (low).
   `R/multfsusie.R:490`, `tol = 1e-6` and `cal_obj = FALSE` by default. With
   `cal_obj = FALSE`, the convergence check uses something other than ELBO
   (alpha-difference); even when on, the ELBO computation itself depends on
   the residual variance estimator from candidate 2.

The Phase 5 investigation will need scenario-by-scenario simulation to
attribute the observed miscalibration to one or several of these.

## 9. Code quality observations

Elegant:

- Wavelet scaling and DWT are cleanly separated (`R/multfsusie.R:308`).
- Functional prior delegated to fsusieR's `init_prior.default` so the
  per-scale machinery is not reimplemented (`R/operation_on_multfsusie_prior.R:82`).
- Inverse DWT is constructed step by step with explicit D / C assembly
  (`R/operation_on_multfsusie_obj.R:2087-2093`); the reconstruction logic is
  visible.
- Post-processing dispatch to smash / TI / HMM is a clean switch
  (`R/operation_on_multfsusie_obj.R:1456-1487, 2401-2477`).
- EM convergence is logged with old/new likelihood
  (`R/EM.R:81-120`), which makes inner-loop debugging tractable.

Warts:

- `pip` is not refreshed after CS filtering (section 8.1 above).
- `R2 / n` residual variance with no df correction (section 8.2).
- lBF aggregation across modalities lacks a weighting term (section 8.5).
- `init_var_multf()` averages per-column variance estimates across positions
  (`R/operation_on_multfsusie_obj.R:606`); per-position information is
  collapsed up front.
- Wavelet coefficients are `cbind`'d into a single matrix and indexed by an
  external list (`list_indx_lst`); structure mismatches between the two are
  silently invalid.
- `discard_cs()` removes from `alpha` but not from `pip` (`:268-291`),
  matching the bug above and adding to the inconsistency.
- Long block of commented-out lfsr code (`R/operation_on_multfsusie_obj.R:1716-1735`):
  `# fsusieR::cal_clfsr(...)`. The return claims to populate lfsr fields,
  but the path is dead.
- `nullweight = 0.7` (`R/multfsusie.R:184`) is a magic default unrelated to
  susieR's value, scaled later by `max(K)` in `R/EM.R:65`. No documentation
  for the choice.
- `sigma2$sd_f` / `sigma2$sd_u` store variances despite the name `sd`
  (`R/operation_on_multfsusie_obj.R:520, :525`); the naming inconsistency
  propagates throughout.

## 10. Evidence pointers

- `R/multfsusie.R:172` - public signature.
- `R/multfsusie.R:308-316` - wavelet pipeline.
- `R/multfsusie.R:489-589` - IBSS loop.
- `R/operation_on_multfsusie_obj.R:508-644` - object initialization.
- `R/operation_on_multfsusie_obj.R:1046-1096` - `get_ER2`.
- `R/computational_routine.R:395-431` - residual variance.
- `R/operation_on_multfsusie_obj.R:1990-2003` - PIP formula.
- `R/operation_on_multfsusie_obj.R:1447-1454` - PIP-then-filter ordering.
- `R/operation_on_multfsusie_obj.R:268-291` - `discard_cs`.
- `R/EM.R:48-50` - lBF aggregation across modalities.
- `R/operation_on_multfsusie_prior.R:81-96` - functional prior init.
- `R/utils_wavelet_transform.R:24-50` - DWT2.
- `R/operation_on_multfsusie_obj.R:2074-2099` - inverse DWT.
- `R/operation_on_multfsusie_obj.R:2301-2477` - post-processing dispatch.
- `R/operation_on_multfsusie_obj.R:1834-...` - convergence check.
- `R/utils.R:389-415` - `init_var_multf`.

## Implications for mfsusieR

The wavelet pathway is the part of mvf.susie.alpha that has to be ported.
Forward DWT, scale-indexed mixture priors, scale-aware CS construction, and
inverse DWT are all reusable concepts; the implementations should be
reorganized into the S3 dispatch shape of mvsusieR (data class
`multifunctional_individual` carries Y_f / Y_u and the DWT cache,
`compute_ser_statistics.*` reads the cache and calls per-scale prior code
delegated to fsusieR or reimplemented). The PIP-then-filter ordering bug is
a hard rule for the port: PIP must be a deterministic function of the
*final* alpha state, not the alpha state mid-pipeline. The residual-variance
estimator should be revisited under Phase 5 with at least the option of a
per-scale or DOF-corrected variant, since the global-per-trait variance is
the most plausible structural source of FDR miscalibration. The
multi-modality lBF aggregation needs a principled weighting before it can
be considered ported, not just rewritten.
