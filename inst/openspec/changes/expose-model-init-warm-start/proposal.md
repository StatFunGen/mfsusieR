# Expose `model_init` for warm-started fits

## Why

A long-running `mfsusie()` fit on hundreds of regions may need
to be checkpointed between calls (out-of-budget walltime, server
preemption) and resumed without restarting from scratch. Both
upstream functional fine-mapping packages support this via a
`multfsusie.obj = m_prev` argument on the legacy entry point.
mfsusieR delegates to susieR's IBSS workhorse, which already
supports warm-start via `model_init`, but the
`ibss_initialize.mf_individual` S3 override explicitly drops
the user-supplied `params$model_init` (the source comment
records "Model_init / warm-start is deferred"). The
`L_greedy` outer loop on the susieR side passes a per-round
`model_init` from one round to the next; this is silently a
no-op on the mfsusieR S3 path, blocking the cross-round warm-
start that susieR designed in.

The current `vignettes/mfsusie_long_running_fits.Rmd` carries
a "planned extension" caveat for this. This change implements
the extension.

## What changes

### 1. Public argument

Add `model_init = NULL` to `mfsusie()` and forward through
`fsusie()`. Roxygen documents the argument as a previously
returned `mfsusie` fit object whose
`alpha`, `mu`, `mu2`, `KL`, and `sigma2` state is used to
seed the next IBSS run.

### 2. `ibss_initialize.mf_individual` lifts the no-op

The override in `R/ibss_methods.R` currently sets the per-
effect state to zero ignoring `params$model_init`. The new
behavior: when `params$model_init` is non-`NULL`, copy
`alpha`, `mu`, `mu2`, `KL`, and `sigma2` from the supplied
fit into the working model object before the IBSS loop
starts. The L of the supplied fit must equal the L requested
in the call; otherwise the function errors with a clear
message.

The fix mirrors `ibss_initialize.default` in susieR.

### 3. `L_greedy` cross-round warm-start now honored

With (2) in place, the susieR `L_greedy` workhorse path that
threads `model_init` from round `k` to round `k + 1` flows
through the mfsusieR S3 override correctly. No mfsusieR-side
plumbing change at the workhorse level is required.

### 4. Tests

`tests/testthat/test_warm_start.R` (new):

- Convergence-shortcut: cold fit on a small fixture takes
  `niter_cold` iterations; resuming the same data from the
  cold result via `model_init = m_cold` takes `niter_warm
  <= 1` iteration. The differential is the test invariant.
- Bit-identity vs susieR scalar M=1, T_1=1 case with
  `s_init = m_prev`.
- Reference test against
  `mvf.susie.alpha::multfsusie(..., multfsusie.obj = m_prev)`
  on a multi-outcome fixture; tolerance set per the existing
  `mvf.susie.alpha` C3 contract (Pattern A allowed for the
  documented sigma2 bug).
- Error path: `model_init` with `L != requested L` errors
  cleanly.

### 5. Vignette

`vignettes/mfsusie_long_running_fits.Rmd` drops the
"planned extension" caveat. Adds a worked example: a 4-
iteration cap, save the partial fit, resume with
`model_init = partial_fit`, observe convergence in 1 - 2
extra iterations.

### 6. Refactor-exceptions ledger

Entry for `mvf.susie.alpha/R/multfsusie.R::multfsusie.obj`:
"replaced-by-`R/ibss_methods.R::ibss_initialize.mf_individual`
honoring `params$model_init`. The legacy `multfsusie.obj`
plus `max_step` workaround is retired in favor of the
single `model_init` argument and the single `max_iter`
budget."

## Impact

- New: tests/testthat/test_warm_start.R, refactor-exceptions
  entry, public-API docs.
- Changed: `R/mfsusie.R`, `R/fsusie.R` signatures;
  `R/ibss_methods.R` override; vignette.
- C2/C3 contracts: warm-start agreement vs upstream is added
  as a new sub-contract on the `multfsusie.obj` use case.
- DESCRIPTION: no changes.
- Specs: new
  `inst/openspec/specs/mf-warm-start/spec.md`.

## Out of scope

- A separate `max_step` per-call iteration cap. Our single
  `max_iter` continues to mean "iterations from current
  state to convergence"; legacy `max_step` was a workaround
  for the absent `model_init`. Once `model_init` lands the
  legacy split is unnecessary.
- Mid-call resume from disk; the user is responsible for
  saving the fit object between calls.
