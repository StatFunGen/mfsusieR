# mf-post-processing capability

## ADDED Requirements

### Requirement: `mf_post_smooth()` is the unified post-processing entry point

`mf_post_smooth()` SHALL be the only public function for adding
smoothed effect curves and credible bands to an `mfsusie` /
`fsusie()` fit. The dispatch axis is the `method` argument.

#### Scenario: signature

`mf_post_smooth(fit, method, level = 0.95, threshold_factor = 1)`
SHALL accept exactly one fit object and one method string.

The function SHALL NOT take `X` or `Y` arguments. Every method
SHALL operate using only the fit and what the fit carries
(`fit$residuals`, `fit$mu`, `fit$mu2`, `fit$alpha`, `fit$dwt_meta`,
and the per-effect lead-variable column saved by the public
`mfsusie()` wrapper as `fit$lead_X`, see below).

#### Scenario: required fit fields

`mf_post_smooth()` SHALL emit a clear error when any of the
following fields is missing:

- `fit$residuals` -- a length-`M` list of `n x T_basis[m]`
  wavelet-domain residual matrices (already populated when
  `mfsusie(... save_residuals = TRUE)`, the default).
- `fit$lead_X` -- a length-`L` list; entry `l` is a length-`n`
  numeric vector containing the centred / scaled column of `X`
  at the lead variable of effect `l`. Populated by `mfsusie()`
  at finalize time when `save_residuals = TRUE`. Removed when
  `save_residuals = FALSE`.
- `fit$mu`, `fit$mu2`, `fit$alpha`, `fit$dwt_meta` -- standard
  fit fields.

### Requirement: methods

#### Scenario: `method = "scalewise"`

Per-scale soft-thresholding of the lead variable's wavelet
posterior mean. For each effect `l` and outcome `m`:

1. `lead <- which.max(fit$alpha[l, ])`.
2. `mean_w <- fit$mu[[l]][[m]][lead, ]` (length `T_basis[m]`).
3. `var_w  <- pmax(fit$mu2[[l]][[m]][lead, ] - mean_w^2, 0)`.
4. Per wavelet scale `s` in `fit$dwt_meta$scale_index[[m]]`,
   compute `sigma_s <- mean(sqrt(var_w[idx_s]))` and soft-
   threshold: `mean_w[idx_s] <- sign(mean_w[idx_s]) *
   pmax(abs(mean_w[idx_s]) - factor * sigma_s * sqrt(2 log T), 0)`.
5. Inverse-DWT the thresholded `mean_w` to position space ->
   `fit$effect_curves[[m]][[l]]`.
6. Inverse-DWT `sqrt(var_w)` (Parseval approximation on the
   orthonormal basis) -> pointwise SD; band is
   `mean_pos +/- z(level) * sd_pos`.

This method does not consume `fit$residuals` or `fit$lead_X`
beyond their existence check.

#### Scenario: `method = "TI"` (cycle-spinning translation-invariant denoising)

Faithful port of `fsusieR::TI_regression.susiF` /
`univariate_TI_regression`. For each effect `l` and outcome `m`:

1. Reconstruct the position-space "isolated" residual:
   `iso_pos <- mf_invert_dwt(fit$residuals[[m]] +
                            fit$lead_X[[l]] %*% mu_lead_w_l[[m]])`,
   where `mu_lead_w_l[[m]]` is the lead variable's posterior mean
   row of the wavelet representation.
2. Compute the **stationary** (i.e. translation-invariant) wavelet
   transform of each row of `iso_pos` via
   `wavethresh::wd(..., type = "station")`, collecting `D` and `C`
   coefficients.
3. Regress wavelet coefficients on `fit$lead_X[[l]]` via
   `crossprod(lead, D)/crossprod(lead) -> Bhat_w`, `Shat_w` from
   the residual variance.
4. Apply scalewise `ashr::ash(Bhat_w[idx_s], Shat_w[idx_s])`
   shrinkage; collect posterior mean and posterior SD per
   coefficient.
5. Reconstruct via `wavethresh::convert()` + `wavethresh::av.basis()`
   (cycle-spinning average across translates).
6. Pointwise band from `wd.var` + `AvBasis.var` per the upstream
   formula.

This method requires `ashr` (Imports) and `wavethresh` (Imports;
already present).

#### Scenario: `method = "HMM"` (HMM-based denoising)

Faithful port of `fsusieR::HMM_regression.susiF` /
`univariate_HMM_regression` / `fit_hmm`. For each effect `l` and
outcome `m`:

1. Reconstruct `iso_pos` as in TI step 1.
2. For each position `t`, regress `iso_pos[, t]` on
   `fit$lead_X[[l]]` to obtain `(Bhat_t, Shat_t)`.
3. Fit a hidden Markov model on `(Bhat_t, Shat_t)` per
   `fsusieR::fit_hmm(x, sd, halfK)` to obtain a posterior mean
   `x_post`, lfsr, and per-position log-Bayes-factor.
4. `effect_curves[[m]][[l]] <- x_post * (Y_csd / X_csd)` per
   the upstream formula. `credible_bands[[m]][[l]]` from the
   posterior SD.
5. Populate `fit$lfsr_curves[[m]][[l]] <- lfsr`. The plot
   functions read this when present.

This method requires `ashr` (already in Imports). The `fit_hmm`
helper itself ports cleanly without external dependencies.

### Requirement: lead-variable storage on fit

The public `mfsusie()` wrapper SHALL, after the IBSS workhorse
returns and when `save_residuals = TRUE`, populate
`fit$lead_X` with `length(L)` entries:

`fit$lead_X[[l]] <- data$X[, which.max(fit$alpha[l, ])]`

where `data$X` is the centred + scaled covariate matrix from the
`mf_individual` data class. The cost is `n * L` doubles, small
relative to the fit's other state.

### Requirement: optional `lfsr_curves` slot

When `method = "HMM"`, `mf_post_smooth()` SHALL populate
`fit$lfsr_curves[[m]][[l]]` (length `T_basis[m]` numeric vector)
in addition to `effect_curves` and `credible_bands`. The plot
functions SHALL detect this slot and draw the lfsr curve as a
secondary black trace with a dashed `0.05` reference line, per
`fsusieR::plot_susiF_effect`'s pattern.

## Output contract

After `mf_post_smooth(fit, method = ...)`:

- `fit$effect_curves[[m]][[l]]` SHALL be a length-`T_basis[m]`
  numeric vector on the post-remap grid.
- `fit$credible_bands[[m]][[l]]` SHALL be a `T_basis[m] x 2`
  numeric matrix; column 1 is `lower`, column 2 is `upper`;
  `lower <= mean <= upper` element-wise.
- `fit$lfsr_curves[[m]][[l]]` SHALL be either NULL or a
  length-`T_basis[m]` numeric vector in `[0, 1]`.

## Test obligations

The implementation SHALL include:

- `test_post_smooth_scalewise.R`: smoke test +
  `expect_equal(dim(fit$credible_bands[[1]][[1]]), c(T_basis, 2))`.
- `test_post_smooth_TI.R`: bit-identity vs
  `fsusieR::univariate_TI_regression` on a single-effect fit at
  `<= 1e-8`. Skip when `fsusieR` is not installed.
- `test_post_smooth_HMM.R`: bit-identity vs
  `fsusieR::univariate_HMM_regression` on a single-effect fit at
  `<= 1e-8`. Skip when `fsusieR` is not installed.
