# Tasks

## 1. Implement gate

- [ ] 1.1 Modify `get_objective.mfsusie` in
  `R/individual_data_methods.R` to short-circuit and return
  `NA_real_` when `params$convergence_method == "pip"`.
- [ ] 1.2 Update roxygen for `mfsusie()`'s `convergence_method`
  arg to note that `fit$elbo` is `NA` on the PIP path.
- [ ] 1.3 Update roxygen for `get_objective.mfsusie` (the
  `@noRd` internal docstring) to reflect the gate.

## 2. Tests

- [ ] 2.1 New test in `tests/testthat/test_variance_and_elbo.R`:
  `mfsusie(..., convergence_method = "pip")` produces a fit
  whose `$elbo` is all NA except possibly the initial NA.
- [ ] 2.2 New test: `mfsusie(..., convergence_method = "elbo")`
  produces a finite, monotone-non-decreasing `$elbo` (preserves
  existing behavior).
- [ ] 2.3 New test: PIPs and alpha bit-identical to pre-change
  fit at `tol = 1e-12` on a fixed-seed fixture (parity test
  comparing `convergence_method = "pip"` to a
  `convergence_method = "elbo"` fit on the same data; alpha
  and PIP must agree at convergence regardless of which gate
  guarded the loop).

## 3. Benchmark

- [ ] 3.1 Re-run the n=84, p=3500, M=6, T=128 0/1/2-genotype
  fixture with default config; record per-iter wall time.
  Verify ≥ 30% improvement.
- [ ] 3.2 Re-run the apples-to-apples mvf comparison (n=100,
  p ∈ {500, 1000, 2000}, M=2, T=128, matched
  `max_SNP_EM = p`, `greedy = FALSE`, `backfit = FALSE`).
  Verify mfsusieR is strictly faster per iter at every p.

## 4. Documentation + commit

- [ ] 4.1 Update `NEWS.md` with the perf change and the
  `fit$elbo` semantics under PIP convergence.
- [ ] 4.2 Commit. Single focused commit
  `perf(elbo): skip refresh_lbf_kl on PIP convergence path`.
