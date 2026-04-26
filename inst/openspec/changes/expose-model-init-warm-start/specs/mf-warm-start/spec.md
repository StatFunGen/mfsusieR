# mf-warm-start capability

## ADDED Requirements

### Requirement: `mfsusie()` accepts a `model_init` argument

`mfsusie()` and `fsusie()` SHALL accept `model_init = NULL`.
When non-`NULL`, `model_init` SHALL be a previously returned
fit object whose `alpha`, `mu`, `mu2`, `KL`, and `sigma2`
state seeds the next IBSS run.

#### Scenario: signature

```
mfsusie(X, Y, pos = NULL, L = ..., model_init = NULL, ...)
```

The function SHALL forward `params$model_init <- model_init`
to the IBSS workhorse.

### Requirement: `ibss_initialize.mf_individual` honors `params$model_init`

The S3 override in `R/ibss_methods.R` SHALL copy `alpha`,
`mu`, `mu2`, `KL`, and `sigma2` from `params$model_init` into
the working model when `params$model_init` is non-`NULL`. The
copy SHALL mirror `susieR::ibss_initialize.default`.

#### Scenario: L compatibility check

When `nrow(params$model_init$alpha) != params$L`, the override
SHALL error with a message identifying both values.

### Requirement: warm-start converges fast on the same data

A warm-started fit on the same data as a previously converged fit SHALL converge in at most one additional IBSS iteration.

A fit returned by `mfsusie(X, Y, ..., L)` and immediately resumed via `mfsusie(X, Y, ..., L, model_init = fit)` on the same `(X, Y)` SHALL report `niter <= 1`.

#### Scenario: convergence-shortcut on a small fixture

Cold fit takes `niter_cold >= 4` iterations on a fixed-seed
fixture. The warm fit started from the cold result reports
`niter_warm <= 1`. The numerical state at convergence agrees
to `tolerance <= 1e-12` between the cold continuation and the
warm result.

### Requirement: cross-package warm-start contract

`mfsusie(..., model_init = m_prev)` SHALL match
`mvf.susie.alpha::multfsusie(..., multfsusie.obj = m_prev)`
on the multi-outcome warm-start case at the established
C3 contract tolerance.

#### Scenario: reference comparison

A multi-outcome fixture is fit twice: once via mfsusieR with
`model_init = m_prev`, once via the upstream legacy entry
point with `multfsusie.obj = m_prev`. The numeric outputs
agree at the same tolerance the existing C3 contract uses for
non-warm-start runs (Pattern A allowed for the documented
sigma2 deviation).
