# Tasks

## 1. Code change

- [ ] 1.1 Add `null_prior_weight` argument to
  `init_scale_mixture_prior_default()` in
  `R/prior_scale_mixture.R`. Default 2.
- [ ] 1.2 Replace the hardcoded `c(0.8, 0.2/(K-1), …)` line at
  `R/prior_scale_mixture.R:101` with
  `c(pi_null, rep((1 - pi_null) / (K - 1), K - 1))` where
  `pi_null = null_prior_weight / (K + 1)`.
- [ ] 1.3 Plumb `null_prior_weight` through
  `mf_prior_scale_mixture()` to `init_scale_mixture_prior_default()`.

## 2. Tests

- [ ] 2.1 Identify which reference tests in `tests/testthat/`
  break under the new init.
- [ ] 2.2 For each broken test, add a swap hack that overwrites
  `fit$prior$G_prior[[m]][[s]]$fitted_g$pi` to the upstream
  `c(0.8, 0.2/(K-1), …)` vector before the comparison.
  Inline comment cites the test-hack pattern.

## 3. Documentation

- [ ] 3.1 Concise roxygen note in
  `init_scale_mixture_prior_default()` explaining the
  Dirichlet-init formula tied to `null_prior_weight`.
- [ ] 3.2 Ledger entry in `inst/notes/refactor-exceptions.md`.

## 4. Verification

- [ ] 4.1 Run `devtools::test()`; expect 0 failures.
- [ ] 4.2 Quick smoke fit on a fixture; verify `fit$niter` and
  `fit$pip` are within tolerance of pre-change values.
- [ ] 4.3 Commit.
