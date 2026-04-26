# Add `small_sample_correction` for Johnson-t marginal Bayes factor

## Why

The default Wakefield Normal marginal Bayes factor in the SER
step under-propagates residual-variance uncertainty into the
per-variable Bayes factor when the sample size is small. This
inflates per-variable PIPs at null variants. The legacy
`mvf.susie.alpha::multfsusie(..., cor_small = TRUE)` switch
addresses this by replacing the Wakefield Normal marginal with
a Johnson 2005 scaled Student-t marginal whose degrees of
freedom track the residual degrees of freedom of the
per-outcome OLS fit (`df = n - 1`).

The Stage 4a audit (see `inst/notes/sessions/`) showed that
plumbing through susieR's NIG path is structurally impossible:
mfsusieR overrides every S3 dispatch point that susieR uses to
reach NIG, so `params$use_NIG = TRUE` would be a no-op. A
wavelet-mixture-prior NIG marginal BF is research work, not a
port. The legacy Johnson-t is implementable, has a known
formula, and is what the original mvf.susie.alpha and fsusieR
both expose for the same use case.

## What changes

### 1. Public argument

```
mfsusie(X, Y, ..., small_sample_correction = FALSE)
```

Boolean flag, default `FALSE`. Forwarded through `fsusie()` via
`...`. When `TRUE`, the SER step replaces the Normal marginal
mixture BF with a Student-t marginal mixture BF using
`df = data$n - 1`.

### 2. Kernel

A new R-only kernel `mixture_log_bf_per_scale_johnson` in
`R/posterior_mixture.R` mirrors the structure of the C++
`mixture_log_bf_per_scale_cpp` Normal kernel but uses
`LaplacesDemon::dstp` for each mixture component. Posterior
moments (`mixture_posterior_per_scale`) are unchanged: the
correction acts on variable selection probabilities only.

### 3. Loglik branch

`loglik.mf_individual` reads `params$small_sample_correction`
and `params$small_sample_df` and dispatches to
`mixture_log_bf_per_scale_johnson` instead of
`mixture_log_bf_per_scale` when the flag is set.

### 4. Dependency

`LaplacesDemon` is added to `DESCRIPTION` Imports. It is
conda-forge `noarch`, so no platform-specific gating is
needed; `pixi.toml` adds `r-laplacesdemon = "*"` to the core
dependencies.

### 5. Tests

`tests/testthat/test_small_sample_correction.R` (new):

- Default path (`small_sample_correction = FALSE`) is
  bit-identical to a call without the argument.
- Johnson-t path runs to convergence on small `n`.
- Structural test: at `n = 40` with a known signal variant,
  Johnson-t recovers the signal (`PIP > 0.9`) and the
  aggregate PIP at null variants is not larger than under the
  Wakefield kernel.
- Argument validation: rejects non-logical input.
- Fidelity test: per-variable LBFs from
  `mixture_log_bf_per_scale_johnson` summed across scales
  match `fsusieR::log_BF` with `df = n - 1` at
  `tolerance <= 1e-12`.

### 6. Refactor-exceptions ledger

Entry mapping `mvf.susie.alpha::multfsusie(..., cor_small =
TRUE)` to `mfsusie(..., small_sample_correction = TRUE)`. The
upstream argument name was opaque; the mfsusieR name describes
the purpose.

## Out of scope

- Porting susieR's NIG path. The audit established that
  plumb-through is impossible (every NIG hook is shadowed by
  an `mf_individual` override) and that a wavelet-mixture-prior
  NIG marginal BF is a separate research contribution.
- Changing posterior moments under the Johnson-t kernel. The
  upstream `cor_small` switch only affects the per-variable
  BF; the conditional posterior given inclusion is still
  Normal. This matches mvf.susie.alpha and fsusieR.

## Impact

- New: `R/posterior_mixture.R::mixture_log_bf_per_scale_johnson`,
  `tests/testthat/test_small_sample_correction.R`,
  refactor-exceptions entry.
- Changed: `R/mfsusie.R` signature; `R/individual_data_methods.R`
  loglik dispatch; `R/mfsusieR-package.R` import;
  `DESCRIPTION` Imports; `pixi.toml`.
- Specs: new
  `inst/openspec/specs/mf-small-sample-correction/spec.md`.
- DESCRIPTION: adds `LaplacesDemon` to Imports.
