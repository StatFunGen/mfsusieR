# Replace hardcoded prior-init `(0.8, 0.2/(K-1), …)` with `null_prior_weight`-driven Dirichlet init

## Why

The data-driven prior init in `R/prior_scale_mixture.R::init_scale_mixture_prior_default()`
(line 101) overwrites `ash`'s fitted mixture weights with the
hardcoded vector `c(0.8, 0.2/(K-1), …, 0.2/(K-1))`. Two issues:

1. The number `0.8` is a magic constant that does not reflect any
   parameter on the public API. mfsusieR already exposes
   `null_prior_weight` (default 2) for the user-supplied-grid path
   via `distribute_mixture_weights()`. The two paths (user-grid
   and ash-driven) currently use different init formulas for the
   same conceptual choice.
2. `null_prior_weight` is not even plumbed into
   `init_scale_mixture_prior_default()`, so there is no way for a
   user to express their prior belief on the null weight when the
   ash path runs (the default).

The fix is to use the existing `null_prior_weight` parameter in
both paths, with the same formula
`pi_null = null_prior_weight / (K + 1)`. One parameter, one
behavior.

## What changes

### 1. Plumb `null_prior_weight` into the ash-driven init path

`init_scale_mixture_prior_default()` gains a `null_prior_weight`
argument (default 2, matching the user-grid path). The hardcoded
`c(0.8, 0.2/(K-1), …)` is replaced with
`c(pi_null, rep((1 - pi_null)/(K-1), K-1))` where
`pi_null = null_prior_weight / (K + 1)`.

`mf_prior_scale_mixture()` passes `null_prior_weight` through to
both `distribute_mixture_weights()` (already does) and
`init_scale_mixture_prior_default()` (new).

### 2. Patch the C2 reference tests with a swap hack

`tests/testthat/test_*.R` files that compare against
`fsusieR::susiF` (which hardcodes the `(0.8, …)` init) currently
rely on bit-fidelity at the prior init level. After the change,
mfsusieR's default init is different. The tests are patched to
overwrite `fit$prior$G_prior[[m]][[s]]$fitted_g$pi` to
`c(0.8, 0.2/(K-1), …)` before the upstream comparison so
bit-identity holds at `tol <= 1e-12`. Same swap-hack pattern
already used for the `X_eff` swap in `test_post_smooth_TI.R`.

### 3. Ledger entry

`inst/notes/refactor-exceptions.md`: new entry documenting the
deliberate divergence from upstream's hardcoded `(0.8, …)` init.
Cites the design decision: `null_prior_weight` controls init pi
in both code paths, no magic constant.

## Impact

- **Behavior**: First IBSS iteration's lBF / alpha differ from
  upstream by `O(null_prior_weight / (K+1) - 0.8)` mass on
  `pi_null`. Subsequent iterations refine pi via mixsqp and
  converge to the same fixed point regardless of init (concave
  objective). Late-iteration alphas / PIPs / CSes are essentially
  unchanged within `tol = 1e-6` on the standard fixtures.
- **API**: `null_prior_weight` is already a public parameter of
  `mfsusie()` (default 2). No new user-facing argument.
- **Tests**: C2 tests patched with swap hack (test-only). Other
  tests unchanged.
- **Documentation**: Brief roxygen comment in
  `init_scale_mixture_prior_default()` explaining why `ash`'s
  fitted pi is discarded (sd-grid is what we want; pi we set
  ourselves based on `null_prior_weight`).
