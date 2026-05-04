# mfsusieR (development version)

- susieR 0.16.1 compat: `track_ibss_fit.mf_individual` now delegates
  to `susieR:::make_track_snapshot` so the per-iteration snapshots
  satisfy the new `is_compact_track_snapshots` validator that
  `ibss_finalize -> make_susie_track_history` runs at finalize.
  Without this, `mfsusie(..., track_fit = TRUE)` errored with
  "fit$trace is not a compact SuSiE track" against susieR >= 0.16.1.
  The snapshot records mfsusieR's `sigma2` (a `list[M]`) as
  `NA_real_` because susieR's helper assumes a scalar; the real
  per-iteration values remain available on `fit$elbo` and
  `fit$sigma2`.
- `pip` now incorporates the `V[l] > prior_tol` filter. The
  per-effect effective slab variance (mean over (m, s) of
  `sum_k pi[l, m, s, k] * var_k`) populates `model$V[l]` after
  every IBSS sweep; effects whose mixture pi has collapsed onto
  the null spike are dropped from the PIP product and zeroed
  by `trim_null_effects.mf_individual` on the rest of the fit.
  Pass `prior_tol = 0` to disable the filter.

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
  "smash"), overwrite_previous = FALSE)` for posterior-
  smoothing of effect curves with credible bands. `"TI"` is
  the default (translation-invariant wavelet denoising via
  cycle spinning); `"HMM"` additionally returns per-position
  lfsr; `"smash"` delegates to `smashr::smash.gaus`
  (Suggests). The `"scalewise"` pointwise variance is
  computed via the squared inverse-DWT matrix
  (`Var(pos[t]) = sum_k W^T_{tk}^2 * var_w[k]`).
- Each `mf_post_smooth()` call adds an entry to
  `fit$smoothed[[method]]` rather than overwriting top-level
  slots, so multiple smoothers coexist on the same fit. With
  `overwrite_previous = FALSE` (default), re-applying the
  same smoother errors instead of clobbering. `coef(fit)`
  returns the raw inverse-DWT curves; pass
  `coef(fit, smooth_method = "<name>")` for a smoothed
  version. `mfsusie_plot()` picks a smoother by priority
  `TI > smash > HMM > scalewise` when several are present
  and emits a hint listing the alternatives; pass
  `smooth_method = "<name>"` to plot a specific one.
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