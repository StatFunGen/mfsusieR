# Speed up the per-effect mixture-prior M-step

## Why

`optimize_prior_variance.mf_individual()` is called once per
(effect `l`, IBSS iteration). Inside, for each (modality `m`,
scale `s`) it builds a likelihood matrix
(`mf_em_likelihood_per_scale`) and runs one mixsqp solve. Each
build is dominated by quantities that do not depend on `l`:

- `shat2_m = outer(1/xtx_diag, sigma2_per_pos)` (depends only on
  data-fixed `xtx_diag` and per-iter `sigma2`).
- `sdmat = sqrt(outer(svec^2, sd_grid^2, +))` (depends on
  `shat` and the init-fixed sd-grid).

Both are recomputed inside every effect call, ~`L`-fold redundant
within one IBSS iteration. The L-mat builder also has an `apply()`
NaN-imputation pass that loops in R rather than vectorising.

A separate concern is M-step input scale: every effect's M-step
solves mixsqp on all `p` SNPs even though SuSiE's per-effect alpha
posterior is concentrated on a handful. SNPs with `alpha_j ≈ 0`
contribute nothing to the M-step gradient and could be safely
dropped from the L-mat input.

## What changes

Stage the work into measurable, independently-revertible steps,
each with a profvis snapshot before and after:

1. **Profile baseline**. profvis on a representative `mfsusie()`
   fit (`practical_dataset` or equivalent). Save flamegraph under
   `bench/profiling/flamegraphs/` with the commit hash. Identify
   the top three hot paths.
2. **Lift loop invariants** out of the per-effect M-step. Cache
   `sigma2_per_pos`, `shat2_m`, and `sdmat[m, s]` in
   `model$em_cache`, repopulated by
   `update_variance_components.mf_individual()` (the per-iter
   step that touches `sigma2`). The M-step reads from the cache
   instead of recomputing. Numerical identity at `tol = 1e-12`.
3. **Vectorise NaN-imputation** in `mf_em_likelihood_per_scale`.
   Replace the R-level `apply(log_L, 2, function(col) …)` with a
   pure matrix operation. Numerical identity at `tol = 1e-12`.
4. **Adaptive variant subsetting**. Inside
   `optimize_prior_variance.mf_individual()`, drop SNPs with
   `alpha[l, j] < eps` (default `eps = 1e-6`) from the L-mat
   input. Truncation error bounded analytically by
   `sum_{j outside} alpha_j · max_k(L_jk)`, well under
   floating-point precision for typical sparse posteriors.
   Tests at `tol = 1e-10`.
5. **C++ port if profile says so**. Only if (1)-(4) leave
   `mf_em_likelihood_per_scale` as a hot path on the post-fix
   profile. Pure-R reference comparison at `tol = 1e-10` per
   CLAUDE.md Phase 7 rule.

Each step is a separate commit on this OpenSpec change.

## Acceptance criteria

- End-to-end fit time of `mfsusie()` on the
  `practical_dataset` reduces by `≥ 30%` after step 4 (target
  set conservatively; flamegraph drives the actual achievable
  number).
- All steps preserve `tol ≤ 1e-10` numerical agreement with the
  pre-change output on the standard test fixtures.
- Reference tests at `tol = 1e-12` (C2, C3) continue to pass
  modulo any port-source-bug fix already documented.

## Impact

- `optimize_prior_variance.mf_individual()` rewrites internal but
  the public API is unchanged.
- `model` object gains an `em_cache` slot (internal, not
  user-facing).
- `update_variance_components.mf_individual()` repopulates the
  cache after sigma2 changes; one extra computation per IBSS
  iter, recovers L-fold redundancy avoidance.
- `mf_em_likelihood_per_scale` gains a vectorised NaN-imputation
  branch.
- Adaptive subsetting is gated by an `eps` parameter (default
  `1e-6`); user can disable by setting `eps = 0`.

## Out of scope

- Recovering ash's empirical pi at init (separate change).
- Reorganising the per-(m, s) loop structure (current order is
  manuscript-aligned).
- mixsqp control parameters (left at upstream defaults).
