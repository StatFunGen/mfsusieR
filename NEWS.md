# mfsusieR 0.0.1 (2026-04-25)

First release.

`mfsusieR` is now official implementation of the fSuSiE
and mfSuSiE models on the `susieR` backbone. It unifies
the previously separate `fsusieR` and `mvfsusieR` codebases under
one package; both methods, as described in the corresponding
manuscripts, are now distributed through `mfsusieR`.

What is in this release:

- `mfsusie()` for multi-outcome functional fine-mapping; each
  outcome can be a scalar or a curve sampled at a grid of
  positions.
- `fsusie()` for  single-outcome functional fine-mapping 
  (M = 1, T_1 > 1) with the conventional `(Y, X, pos, ...)`
  argument order.
- S3 methods: `predict.mfsusie()`, `coef.mfsusie()`,
  `fitted.mfsusie()`, `summary.mfsusie()`, `print.mfsusie()`,
  and `plot.mfsusie()`.
- `mf_post_smooth(fit, method = c("TI", "scalewise", "HMM",
  "smash"))` for posterior-smoothing of effect curves with
  credible bands. `"TI"` is the default (translation-
  invariant wavelet denoising via cycle spinning); `"HMM"`
  additionally returns per-position lfsr; `"smash"` delegates
  to `smashr::smash.gaus` (Suggests). The `"scalewise"`
  pointwise variance is computed via the squared inverse-DWT
  matrix (`Var(pos[t]) = sum_k W^T_{tk}^2 * var_w[k]`), which
  is the textbook variance propagation for an orthonormal
  linear operator with diagonal input covariance; the prior
  `(invert_dwt(sqrt(var_w)))^2` formula confused the linear
  combination of standard deviations for the position SD.
- `mf_adjust_for_covariates(Y, Z, X = NULL,
   method = c("wavelet_eb", "ols"))` for pre-fit covariate
  adjustment of a functional response. `method = "wavelet_eb"`
  fits a wavelet-domain empirical-Bayes regression of Y on Z;
  `method = "ols"` is the closed-form linear regression covariate
  adjustment using ordinary least square method.
- `mfsusie_plot()` and `mfsusie_plot_lfsr()` for visualizing
  PIPs, credible-set membership, per-CS effect curves with
  optional credible bands, and per-CS local false sign rates.
- Bundled simulated example datasets: `dnam_example`,
  `rnaseq_example`, `multiomic_example`, `gtex_example`. Each
  ships with a build script under `data-raw/`.