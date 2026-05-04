# 2026-05-03 PR plan: mu/mu2 storage + per_scale_normal benchmark

Date: 2026-05-03
Scope: design and PR plan for the `save_mu_method` storage policy
on `mfsusie()` (closes issue #7) and the surrounding 6-grid
benchmark over `prior_variance_scope` x `wavelet_qnorm` x
`mixture_null_weight`.

## Design

`save_mu_method = c("complete", "alpha_collapsed", "lead")` on
`mfsusie()` and the `fsusie()` wrapper, default `"complete"`.

| Mode | `mu[[l]][[m]]` shape | coef / mf_post_smooth numerics | predict(newx) | per-variant lfsr | model_init |
|---|---|---|---|---|---|
| complete           | p x T_basis[m] | exact | works | works | works |
| alpha_collapsed    | 1 x T_basis[m]  (= alpha %*% mu_full)         | numerically equivalent to complete (1e-12) | errors | errors | errors |
| lead               | 1 x T_basis[m]  (= mu_full[j*, ], j* = which.max(alpha[l, ])) | cheap lead-variable summary, biased toward j* | errors | errors | errors |

To recover raw-X coef under `alpha_collapsed` the fit also carries
`fit$coef_wavelet[[l]][[m]] = alpha %*% (mu_full / csd_X)` (1 x T)
because per-j csd_X scaling cannot be recovered after the
alpha-collapse. Under `lead` the fit carries `fit$top_index[l]`.

A post-fit helper `mf_thin(fit, method)` performs the same trim,
so a caller can keep a `complete` checkpoint for warm-start and a
thinned distribution copy.

## Benchmark grid

Six cells, intended to exercise `prior_variance_scope` x
`wavelet_qnorm` x `mixture_null_weight`:

```r
bench_grid <- rbind(
  expand.grid(wavelet_qnorm        = c(FALSE, TRUE),
              prior_variance_scope = "per_scale",
              mixture_null_weight  = c(0.05, 0)),
  data.frame  (wavelet_qnorm        = c(FALSE, TRUE),
               prior_variance_scope = "per_scale_normal",
               mixture_null_weight  = NA_real_)
)
```

Driver: `inst/bench/profiling/benchmark_per_scale_normal_6grid.R`
(Gaussian baseline, 30 fits) and
`inst/bench/profiling/benchmark_heavy_tailed_null_6grid.R`
(heavy-tailed signal + null no signal, 60 fits). Per-cell metrics:
empirical FDR, power, n_disc, cs_count, cs_purity, niter,
runtime, fit_size_mb. The benchmark fits use
`save_mu_method = "alpha_collapsed"` so each saved fit is small.

Results memo:
`inst/notes/sessions/2026-05-03-2314-per-scale-normal-baseline-results.md`.

## Default-value table

Defaults are unchanged in PR-1 (the storage policy is opt-in;
behaviour for users who do not pass `save_mu_method` is identical
to before).

| Arg | Default (mfsusieR 0.0.2) | Notes |
|---|---|---|
| `save_mu_method`         | `"complete"` (new arg) | opt-in to the 1D modes |
| `prior_variance_scope`   | `"per_outcome"`        | benchmark covers per_scale and per_scale_normal but does not motivate a switch |
| `wavelet_qnorm`          | `FALSE`                | rebenchmark with `TRUE` per scenario; results memo discusses |
| `mixture_null_weight`    | `NULL` (resolves to 0.05) | per_scale + 0.05 is the only well-calibrated `per_scale` cell across scenarios |
| `null_prior_init`        | `0`                    | only an init; the EM M-step overwrites within a few iterations |
| `small_sample_correction`| `FALSE`                | sensitivity only; see issue #8 |
| `L`, `L_greedy`, `greedy_lbf_cutoff` | `20`, `5`, `0.1` | unchanged |

## Three-layer warm start (terminology pin)

1. mixsqp internal warm start: on by default
   (`be2722e perf(ibss): mixsqp warm start + ser_cache`).
2. `L_greedy` ramp 5 -> L = 20: on by default.
3. operational warm start via `model_init` (cheap fit -> expensive
   refit): works only with `save_mu_method = "complete"`. Both 1D
   modes drop the per-variant axis the SER step needs at iter 0;
   the guard in `R/save_mu_method.R::mf_apply_save_mu_method`
   stop()s when a thinned fit is supplied as `model_init`.

## susieR 0.16.1 compat (PR-1 commit fix(track))

susieR 0.16.1 added `is_compact_track_snapshots(fit$trace)` inside
`ibss_finalize -> make_susie_track_history`. The validator
requires each tracking list element to be
`list(alpha = data.frame, effect = data.frame, iteration = data.frame)`.
The old `track_ibss_fit.mf_individual` wrote
`list(alpha = matrix, sigma2 = list, pi_V = list, elbo = scalar)`,
which the new validator rejects. Symptom on `track_fit = TRUE`:
"fit$trace is not a compact SuSiE track" against susieR >= 0.16.1.

Fix: `track_ibss_fit.mf_individual` delegates to
`susieR:::make_track_snapshot(model, iteration)` (cached as a
package-level binding via `R/zzz.R::.onLoad`). `model$sigma2` is
sanitised to `NA_real_` for the snapshot copy because mfsusieR's
`sigma2` is `list[M]` whereas `susieR:::track_scalar` assumes a
length-1 numeric. The real per-iteration `sigma2` stays on the
fit; only the trace records `NA`. `fit$elbo` is unaffected.

Pre-fix `tests/testthat/test_ibss_methods_branches.R:35:3` failed
on both `main` and `fix-mu-storage` against susieR 0.16.1. Post-fix
the test passes and the full suite is green.

## Vignette / test parameter audit (2026-05-03)

Confirmed all vignette and test code uses current public-API
parameter names (`mixture_null_weight`, `prior_variance_scope`,
`wavelet_qnorm`, `wavelet_magnitude_cutoff`, `null_prior_init`,
`alpha_thin_eps`), not legacy names from `mvf.susie.alpha` /
`fsusieR` (`null_weight`, `null_prior_weight`, `max_SNP_EM`,
`low_count_filter`, `cor_small`, `maxit`, `cov_lev`, `min_purity`,
`init_pi0_w`, `posthoc`).

Findings:
- `vignettes/post_processing.Rmd:269-276` uses `nullweight = 300`.
  This is `ashr::ash()`'s own `nullweight` argument forwarded
  through `mf_post_smooth(..., method = "ash", nullweight = ...)`'s
  `...`, not a mfsusieR null-prior knob. Action: tighten the
  surrounding prose to disambiguate from `mfsusie()`'s
  `mixture_null_weight` (PR-1 ships this prose).
- Legacy names in `test_fsusier_degeneracy.R`,
  `test_mvf_alpha_degeneracy.R`, `test_cs_parity_fsusier.R`
  (`backfit`, `maxit`) are arguments of the apple-to-apple
  comparison targets (`mvf.susie.alpha::multfsusie`,
  `fsusieR::susiF`). These must NOT be renamed; the equivalence
  contract requires calling the upstream signatures literally.
- `test_public_api_naming.R:15-25` already lists
  `max_SNP_EM`, `max_scale`, `init_pi0_w`, `min_purity` etc as
  forbidden on the mfsusieR public API.
- `test_per_scale_normal_degeneracy.R:684-706` exercises
  `mixture_null_weight` correctly and validates the warning
  emitted under `prior_variance_scope = "per_scale_normal"`.

No code changes to vignettes or tests for parameter names beyond
the PR-1 prose tightening.

## Status

- [x] Pull latest mfsusieR main (at `fa53ad2`)
- [x] Verify `mixsqp_null_penalty` is not in current main
- [x] Verify William's `fitted_func` storage pattern in
      `mvf.susie.alpha`
- [x] Update susieR from GitHub master in local R lib
      (0.16.1, commit `220a191`)
- [x] Confirm `r-ebnm` in pixi env (1.0.55)
- [x] PR-1 branch `fix-mu-storage`, full implementation, full test
      suite green (0 fail / 0 regression)
- [x] PR-1 6-cell baseline benchmark + heavy-tailed + null
      follow-up benchmarks complete (results memo separate)
