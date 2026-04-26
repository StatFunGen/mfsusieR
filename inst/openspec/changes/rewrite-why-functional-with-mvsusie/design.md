# Design — `fsusie_why_functional` rewrite

## D1. Simulation

```
n      = 100             # rows of N3finemapping$X[1:100, ]
T      = 128             # functional response length
p      = ncol(X)         # full SNP set after column-zero filter
pos    = 700             # causal SNP index (in unfiltered X)
effect = 1.2 * cos((1:T) / T * 3 * pi)
effect[1:40] = 0          # zero on positions 1..40 (truncated)
y[i, ] = X[i, pos] * effect + rnorm(T)
```

Same recipe as `mvsusieR/inst/prototypes/diagnose_mvsusie_issue_61.R`.

## D2. Methods compared

```
                COSINE EFFECT ON POSITIONS 41..128
                ════════════════════════════════════
                            ┌──────────────┐
                            │   y (n x T)  │
                            └──────┬───────┘
                                   │
              ┌────────┬───────────┼─────────────┬──────────┐
              ▼        ▼           ▼             ▼          ▼
          Per-CpG    Top-PC      mvSuSiE       mvSuSiE    fSuSiE
           SuSiE     SuSiE       (default)     (true U)
           (T fits)  (5 fits)    L=1, est V    L=1, V=I    L=1
              │        │           │             │          │
              │        │           │             │          │
              ▼        ▼           ▼             ▼          ▼
         Fragments  Mostly      All-zero    Correct CS,  Correct CS,
         with low   fail (PCs   PIPs        smooth eff   smooth eff
         PIPs at   not aligned                           curve via
         most pos  to signal)                            mf_post_smooth
```

For each method we report:
- Credible sets (`fit$cs` or equivalent).
- PIP at the true causal SNP (position 700, or its index after
  column filtering).
- For methods 4 and 5: an effect-curve overlay plot. `coef()`
  on the mvsusie fit, `effect_curve()` (via mfsusie_plot's
  internals) on the fsusie fit. Truth in black, recovered
  curve in CS color.

## D3. Plot strategy

Per-CpG susie: a single ggplot2-or-base-R bar per SNP, height
= `max_t pip_t[snp]`. No effect curve to plot (per-CpG fits
return scalar effects per position).

Top-PC susie: 5 PIP plots side by side or in a 2x3 grid.

mvsusie effect-curve plot: inline `plot()`. The mvsusie fit
`fit_mv` carries `coef(fit_mv)` returning `(p+1) x T`. The
lead SNP's effect curve is `coef(fit_mv)[lead + 1, ]`. Overlay
truth as `plot(effect, type = "l")` then `lines(coef_curve)`.

fsusie effect-curve plot: `mfsusie_plot(fit_s)` after
`mf_post_smooth(fit, method = "TI")`.

## D4. Issue 61 configurations

From `mvsusieR/inst/prototypes/diagnose_mvsusie_issue_61.R`:

```r
# Failing config (default canonical prior + estimated residual variance)
fit_default <- mvsusie(X, y, L = 1)

# Working config (true outer-product prior + fixed identity residual var)
U_true <- outer(effect, effect)
prior_true <- create_mixture_prior(
  mixture_prior = list(matrices = list(U_true), weights = 1))
V_true <- diag(T)
fit_true_prior <- mvsusie(X, y, L = 1,
  prior_variance = prior_true,
  residual_variance = V_true,
  estimate_residual_variance = FALSE,
  estimate_prior_variance = TRUE,
  max_iter = 100)
```

The vignette runs both. The "default" call demonstrates the
all-zero-PIP failure documented in mvsusieR issue 61. The
"true prior" call recovers the correct CS, and we use it to
demonstrate the cross-condition effect curve.

## D5. Excalidraw diagram

`design/methods-comparison-flow.excalidraw` lays out the five
methods as a comparison fan emerging from the shared `(X, y,
truth)` tuple, with success/failure annotations.

## D6. Out of scope

- Any change to `mfsusie_plot()` API.
- Any new R code in mfsusieR.
- Any tests (no numerical test contract; the vignette IS the
  test via the build/render gate).
