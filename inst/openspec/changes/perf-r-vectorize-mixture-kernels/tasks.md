# Tasks

## 1. New helper: vectorized log_dens computation

- [ ] 1.1 Add `mf_em_log_dens_per_scale(bhat_slice, shat_slice,
  sdmat, log_sdmat)` to `R/em_helpers.R`. Returns
  `list(log_dens = (N x K), log_dens_null = length-N)`. Pure
  vectorized R using broadcast arithmetic and the cached
  `sdmat` / `log_sdmat`.
- [ ] 1.2 Unit test in `tests/testthat/test_variance_and_elbo.R`:
  output matches the pre-change cpp `log_dens` at `tol = 1e-12`
  on a fixed-seed fixture.

## 2. Stash log_dens cache on ser_stats

- [ ] 2.1 Modify `compute_ser_statistics.mf_individual` to
  compute the log_dens cache per (m, scope group) using the
  helper from (1) and stash on
  `ser_stats$log_dens_per_m[[m]][[s]]`.
- [ ] 2.2 Verify the cache is read by downstream callers in
  the same SER step; never persisted across effects.

## 3. Replace mf_em_likelihood_per_scale's cpp call with R

- [ ] 3.1 Modify `mf_em_likelihood_per_scale` in
  `R/em_helpers.R`. When `log_dens_cache` is supplied, use
  `L <- exp(log_dens_cache)`. Vectorize the NA imputation
  block (replace the `apply(log_L, 2, ...)` with a
  matrix-level operation using `matrixStats::colMedians`).
- [ ] 3.2 Add `matrixStats` to DESCRIPTION Imports if not
  already (verify; it is currently in Imports per
  `mf_get_ER2_per_position`).
- [ ] 3.3 Test parity: same fit produces identical mixsqp pi
  output at `tol = 1e-12`.

## 4. Replace mixture_log_bf_per_scale_cpp with R aggregator

- [ ] 4.1 In `loglik.mf_individual`, use vectorized R
  log-sum-exp over the cached `log_dens` to compute the
  per-outcome lbf vector. Use `matrixStats::rowMaxs` for the
  stable softmax.
- [ ] 4.2 Verify identical to pre-change at `tol = 1e-12`
  (including the small-sample Johnson-t branch — that one
  still calls `mixture_log_bf_per_scale_johnson` which is
  pure R already).

## 5. Replace mixture_posterior_per_scale_cpp with R aggregator

- [ ] 5.1 In `calculate_posterior_moments.mf_individual`,
  build posterior moments via the closed form using cached
  `log_dens`, `sd_grid`, `shat`, `V`. Vectorized R.
- [ ] 5.2 Test parity: mu and mu2 identical to pre-change at
  `tol = 1e-12`.

## 6. Delete cpp kernels

- [ ] 6.1 Remove the three kernels from
  `src/posterior_mixture.cpp` (or delete the file entirely if
  no other kernels remain).
- [ ] 6.2 Run `cpp11::cpp_register()` to regenerate
  `src/cpp11.cpp` and `R/cpp11.R`.
- [ ] 6.3 Rename the R reference oracles in
  `R/posterior_mixture.R` to drop the `_R` suffix; promote
  to production. Remove redirector wrappers
  `mixture_log_bf_per_scale` and `mixture_posterior_per_scale`
  if they only call the cpp kernel — replace with the new
  vectorized R impls.

## 7. Benchmark + commit

- [ ] 7.1 Re-run the n=84, p=3500, M=6, T=128 0/1/2-genotype
  fixture; verify ≥ 10% improvement on top of
  `perf-skip-elbo-on-pip`.
- [ ] 7.2 Re-run the apples-to-apples mvf comparison
  (n=100, p ∈ {500, 1000, 2000}, M=2). Verify mfsusieR
  strictly faster per iter than mvf at every p.
- [ ] 7.3 Update `NEWS.md` with the cpp-removal note.
- [ ] 7.4 Commit. Single focused commit
  `perf(mixture): replace cpp kernels with vectorized-R + log_dens cache`.

## 8. Documentation

- [ ] 8.1 Update `inst/notes/sessions/<date>-cpp-removal.md`
  documenting the rationale (mvf/fsusieR achieve same with R;
  3-kernel duplication; cache strategy).
