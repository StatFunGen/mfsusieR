# Tasks

## 1. Baseline profvis

- [ ] 1.1 Write `bench/profiling/profile-mfsusie-mstep.R` that
  fits `mfsusie()` on the `practical_dataset` and wraps the call
  in `profvis::profvis()`.
- [ ] 1.2 Run; save `bench/profiling/flamegraphs/<hash>-baseline.html`.
- [ ] 1.3 Identify top three hot paths; record in
  `inst/notes/sessions/<date>-perf-mixsqp-baseline.md`.

## 2. Lift loop invariants

- [ ] 2.1 Add `model$em_cache` slot, repopulated by
  `update_variance_components.mf_individual()` (`sigma2_per_pos`,
  `shat2_m`, `sdmat[m, s]`).
- [ ] 2.2 Refactor `optimize_prior_variance.mf_individual()` to
  read from `model$em_cache` instead of recomputing.
- [ ] 2.3 Verify `devtools::test()` passes at the existing
  tolerances.
- [ ] 2.4 Re-profile; save `bench/profiling/flamegraphs/<hash>-cache.html`.
- [ ] 2.5 Commit.

## 3. Vectorise NaN-imputation

- [ ] 3.1 Replace the `apply(log_L, 2, …)` block in
  `mf_em_likelihood_per_scale` with a pure matrix operation.
- [ ] 3.2 Verify tests pass.
- [ ] 3.3 Re-profile; save flamegraph.
- [ ] 3.4 Commit.

## 4. Adaptive variant subsetting

- [ ] 4.1 Add `mixsqp_alpha_eps` parameter to `mfsusie()`
  (default `1e-6`).
- [ ] 4.2 Inside `optimize_prior_variance.mf_individual()`,
  before building the L-mat, restrict to
  `keep_idx <- which(model$alpha[l, ] > mixsqp_alpha_eps)`.
- [ ] 4.3 Adjust `zeta_l`, `bhat_m`, `shat_m` slices accordingly.
- [ ] 4.4 Re-profile; verify tests pass at `tol = 1e-10`.
- [ ] 4.5 Commit.

## 5. C++ port if needed

- [ ] 5.1 Inspect post-step-4 flamegraph. Skip 5 if
  `mf_em_likelihood_per_scale` is no longer in the top three
  hot paths.
- [ ] 5.2 If kept: cpp11armadillo port with R reference; fidelity
  test at `tol = 1e-10`.
- [ ] 5.3 Commit.

## 6. Final report

- [ ] 6.1 Write summary in `inst/notes/sessions/<date>-perf-mixsqp.md`:
  what was done, baseline → final wall-clock, flamegraph
  before/after, C++ decision.
