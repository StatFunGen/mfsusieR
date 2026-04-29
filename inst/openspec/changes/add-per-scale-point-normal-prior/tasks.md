# Tasks

## Implementation discipline

Unit tests in this plan are *authored* as part of the spec.
Keep them all. **Do not run the full unit-test suite after
trivial mechanical steps** (variable rename, sed-style
replacements, single-line default flips, removing a stale
comment). Running adds minutes to the feedback loop and catches
nothing on those cases. Run the suite only after non-trivial
code lands: a new function, a dispatch change, a new test, or
anything that touches numerical paths.

## 1. Dependency

- [ ] 1.1 Add `ebnm (>= 1.1.0)` to `DESCRIPTION` `Imports`.
  ebnm enters through `susieR`'s closure already
  (`mixsqp`, `ashr`); the marginal weight is `ebnm` itself.
- [ ] 1.2 Add `importFrom(ebnm, ebnm_point_normal)` and
  `importFrom(ebnm, ebnm_point_laplace)` to `NAMESPACE` (or
  let roxygen2 emit them via `@importFrom` tags on the M-step
  helper).
- [ ] 1.3 Confirm `R CMD check` is clean on a fresh tree.

## 2. Public API: `prior_variance_scope` enum extension

- [ ] 2.1 In `R/mfsusie.R`, extend the formal:
  `prior_variance_scope = c("per_outcome", "per_scale",
  "per_scale_normal", "per_scale_laplace")`.
- [ ] 2.2 Update the roxygen `@param prior_variance_scope` to
  describe all four values, the per-scale point-* motivation
  (per-scale data scarcity + intrinsic regularization of the
  2-parameter parametric form), and that the new values are
  fit via `ebnm` (Normal slab vs Laplace slab).
- [ ] 2.3 Forward the value through `mf_prior_scale_mixture()`
  in `R/prior_scale_mixture.R`.

## 3. Internal `prior_class` infrastructure

- [ ] 3.1 In `R/prior_scale_mixture.R`, map the new public scopes
  to `prior_class = "mixture_point_normal_per_scale"` and
  `prior_class = "mixture_point_laplace_per_scale"`. Set
  `class(G_prior) <- prior_class` on the init helper's return
  (single-class tag, mirroring the foreseeable-set cardinality;
  no parent-class hierarchy).

## 4. Init helper for ebnm-backed per-scale prior

- [ ] 4.1 New function in `R/prior_scale_mixture.R`:
  `init_ebnm_prior_per_scale(Y_m, X, scale_index, prior_class,
  ...)`. The helper computes
  `bs <- compute_marginal_bhat_shat(X, Y_m)` once, then per
  scale `s` picks
  `lead_s <- which.max(rowMeans(bs$Bhat[, idx_s, drop = FALSE]^2))`
  and calls the ebnm function selected by `prior_class` on
  `(bs$Bhat[lead_s, idx_s], bs$Shat[lead_s, idx_s])`.
- [ ] 4.2 Store `fitted_g <- fit$fitted_g` and `idx <- idx_s` on
  `G_prior[[s]]`. Record `lead_init_s = lead_s` on the same
  entry (used by tests to assert the lead picker selected a
  signal-bearing variable).
- [ ] 4.3 Set `class(G_prior) <- prior_class` (single-class tag).
- [ ] 4.4 Wire `mf_prior_scale_mixture()` to call this helper
  when `prior_variance_scope ∈ {"per_scale_normal",
  "per_scale_laplace"}`.
- [ ] 4.5 Unit test: shape of returned G_prior (length M, each a
  list of length S_m, each entry a `fitted_g` record with `idx`
  set correctly; class on the outer list is the single
  `prior_class` tag).
- [ ] 4.6 Unit test: marginal-data lead picker selects a
  signal-bearing variable on the synthetic sparse fixture
  (`set.seed(2L); n = 200; p = 30; T_m = 64; signal at variables
  c(7, 18) and times c(20, 44)`). The recorded
  `G_prior[[1L]][[s]]$lead_init_s` should be in `c(7L, 18L)` for
  every scale `s` whose `idx_s` overlaps the signal positions
  `c(20L, 44L)`, and ebnm's fitted `sigma` should be at least 5x
  the per-scale empirical SD of pure-noise variables.

## 5. M-step dispatch, `.opv_<class>` helpers mirroring `.opv_mixsqp`

- [x] 5.1 New helper in `R/individual_data_methods.R`:
  `.opv_ebnm_point(data, params, model, ser_stats, keep_idx,
  zeta_keep, ebnm_fn)`. Signature matches `.opv_mixsqp`'s
  `(data, params, model, ser_stats, keep_idx, zeta_keep)` plus
  an `ebnm_fn` argument.
- [x] 5.2 Body: for each `m`, slice
  `bhat_m <- ser_stats$betahat[[m]][keep_idx, , drop = FALSE]`
  and `shat_m <- sqrt(ser_stats$shat2[[m]][keep_idx, , drop = FALSE])`.
  For each `s` in `seq_along(G_m)`, call
  `ebnm_fn(x = as.vector(bhat_m[, idx, drop = FALSE]),
  s = as.vector(shat_m[, idx, drop = FALSE]),
  g_init = G_m[[s]]$fitted_g,
  fix_g = !isTRUE(params$estimate_prior_variance))`. Write
  `fit$fitted_g` into `G_m[[s]]$fitted_g` and
  `fit$fitted_g$pi` into `model$pi_V[[m]][s, ]`. Multi-variable
  by design (mirrors mixsqp's data shape minus the per-row
  alpha weighting); avoids the lead-only `pi_0` collapse.
- [x] 5.3 Add two thin shims `.opv_ebnm_point_normal` and
  `.opv_ebnm_point_laplace` that delegate to `.opv_ebnm_point`
  with the matching `ebnm_fn`.
- [x] 5.4 Extend the dispatch arm in
  `optimize_prior_variance.mf_individual()` with two new
  `else if` branches keyed on
  `inherits(G_prior[[1L]], "mixture_point_normal_per_scale")`
  and `inherits(G_prior[[1L]], "mixture_point_laplace_per_scale")`.
  Each calls the matching `.opv_ebnm_point_*` helper.
- [x] 5.5 The new helpers do not consume
  `params$mixture_null_weight` (no-op on the parametric form).
  They share the alpha-thinned `(keep_idx, zeta_keep)` already
  computed by the caller using `params$alpha_thin_eps`. No
  deferred-M-step gate: at iter 1 alpha is uniform and ebnm sees
  all p variables, so `pi_0` fits close to 1.
- [x] 5.6 Rename `mixsqp_alpha_eps` to `alpha_thin_eps` in
  `R/mfsusie.R` (formal + roxygen + body forwarding) and
  `R/individual_data_methods.R` (read site). The rename
  reflects that the threshold now drives both M-step solvers.

## 6. `iter_cache` with class-gated slots

- [x] 6.1 Rename `refresh_em_cache.mf_individual()` to
  `refresh_iter_cache.mf_individual()` and the model slot
  `model$em_cache` to `model$iter_cache`. The previous name was
  misleading because the cache is consumed by SER calls (loglik
  + posterior moments) too, not just the M-step.
- [x] 6.2 Drop the `sigma2_per_pos` slot. It was populated by
  the cache build but never read by any consumer; every reader
  calls `mf_sigma2_per_position()` directly.
- [x] 6.3 Add a one-line gate
  `is_mixsqp_prior <- inherits(model$G_prior[[1L]],
  c("mixture_normal", "mixture_normal_per_scale"))` early in the
  function. Always build `iter_cache$shat2`. Build
  `iter_cache$sdmat` and `iter_cache$log_sdmat` only when
  `is_mixsqp_prior` (the K-axis precompute is mixsqp-only;
  ebnm has no K-axis aggregate).
- [x] 6.4 Unit test: with
  `prior_variance_scope ∈ {"per_scale_normal",
  "per_scale_laplace"}`, `model$iter_cache$sdmat` and
  `model$iter_cache$log_sdmat` MUST be NULL after
  `refresh_iter_cache.mf_individual()` runs;
  `model$iter_cache$shat2` MUST be set (one
  `p × T_basis[m]` matrix per outcome `m`).

## 7. Degenerate-case fidelity tests (machine precision)

All `tol = 1e-12` unless noted. Tests live in
`tests/testthat/test_per_scale_normal_degeneracy.R` (kept; the
file already exists).

### 7a. Scalar T=1 vs `susieR::susie`, primary correctness locks

Common fixture: `n = 200`, `p = 50`, two causal SNPs, scalar Y,
`L = 5`, **`max_iter = 200`**, **`tol = 1e-10`**,
**`convergence_method = "elbo"`**,
`estimate_prior_variance = FALSE`,
`estimate_residual_variance = TRUE`,
`L_greedy = NULL`, `null_prior_init = 0`,
`residual_variance_scope = "per_outcome"`.
Reference: same fixture in `susieR::susie()`, matching `max_iter`
and `tol`.

Tolerance protocol:

- **`tol = 1e-12`** on: `alpha`, `pip`, `niter`,
  `lapply(sets$cs, sort)`, `sigma2[[1]]`. These survive ebnm-
  free arithmetic and only depend on the IBSS converging.
- **`tol = 1e-10`** on: `mu` list, `mu2` list, `lbf`, `KL`,
  `tail(elbo, 1)`. These can pick up residual numerical noise
  from the per-iter loglik kernel even at convergence.

- [ ] 7a.1 `prior_variance_scope = "per_outcome"` +
  `prior_variance_grid = sigma2` (length-1), already exists,
  ensure it still passes.
- [ ] 7a.2 `prior_variance_scope = "per_scale"` +
  `prior_variance_grid = sigma2` (length-1).
- [ ] 7a.3 `prior_variance_scope = "per_scale_normal"` +
  `sigma_init = sqrt(sigma2)`. **Primary lock for the new
  ebnm-backed path.** Configured with
  `estimate_prior_variance = FALSE` so ebnm's M-step is
  short-circuited and the prior stays at the fixed Gaussian.

### 7b. Cross-flavor consistency at the degenerate point

Same fixture as 7a. Pairs that should bit-match each other; any
two of the three should agree at `tol = 1e-12`.

- [ ] 7b.1 `per_outcome` ≡ `per_scale` (both with length-1 grid).
- [ ] 7b.2 `per_outcome` ≡ `per_scale_normal`.
- [ ] 7b.3 `per_scale` ≡ `per_scale_normal`.

Transitivity is implied; assert all three for redundancy.

### 7c. Functional T > 1, single-Gaussian-per-scale degenerate

`per_scale` (length-1 grid + null=0 + no EB) and
`per_scale_normal` (null=0 + sigma_init fixed + no EB) both
collapse to a single Gaussian per scale; they should bit-match.

- [ ] 7c.1 T=64 fixture, single causal: `per_scale` ≡
  `per_scale_normal` at `tol = 1e-12` on the full fit shape
  (alpha, mu list, mu2 list, pip, lbf, KL, sigma2, elbo[-1],
  niter, sets$cs, dwt_meta).

### 7d. Vary `L`

The dispatch must be `L`-invariant. Run 7a.3 with
`L ∈ {1, 5, 10}` and verify each matches `susie(L = same)` at
`tol = 1e-12`.

- [ ] 7d.1 L = 1 (single effect): bit-match.
- [ ] 7d.2 L = 5: bit-match (= 7a.3).
- [ ] 7d.3 L = 10: bit-match.

### 7e. ebnm wrapper fidelity (no IBSS)

Lower-level tests that exercise the M-step wrapper directly,
bypassing the IBSS dispatch.

- [ ] 7e.1 Forwarding contract (Normal): the wrapper passes
  `(x, s, g_init, fix_g)` to `ebnm::ebnm_point_normal` exactly
  as constructed; the returned `fitted_g` is unpacked into
  `G_m[[s]]$fitted_g` without modification. Assert by mocking
  `ebnm::ebnm_point_normal` with `local_mocked_bindings()` and
  capturing the call.
- [ ] 7e.2 Forwarding contract (Laplace): same as 7e.1 against
  `ebnm::ebnm_point_laplace`.
- [ ] 7e.3 `g_init` warm-start: the wrapper forwards
  `g_init = G_m[[s]]$fitted_g` to ebnm at every IBSS iter
  beyond the first. Assert by mocking and capturing two
  successive calls; the second call's `g_init` equals the
  first call's returned `fitted_g`.
- [ ] 7e.4 `fix_g` short-circuit: when
  `estimate_prior_variance = FALSE`, the wrapper either passes
  `fix_g = TRUE` to ebnm OR does not call ebnm at all. The
  G_prior slot's `fitted_g` MUST equal the input `fitted_g` at
  `tol = 1e-12` after the M-step.

## 7f. Numerical-property tests (do not require susie reference)

- [ ] 7f.1 **M-step idempotence**: two back-to-back M-step calls
  on the same `(bhat, shat)` input return bit-identical `pi`
  and `sd` (`tol = 1e-14` for Normal, `tol = 1e-12` for
  Laplace; ebnm's optim is deterministic on identical inputs).
- [ ] 7f.2 **Probability constraints**: returned
  `fitted_g$pi` satisfies `length(pi) == 2`, `sum(pi) == 1`
  at `tol = 1e-12`, both components in `[0, 1]`.
- [ ] 7f.3 **Slab non-negativity**: returned
  `fitted_g$sd[2] >= 0` (Normal) and `fitted_g$scale[2] >= 0`
  (Laplace) at `tol = 1e-12`; the spike component (`sd[1]`
  / `scale[1]`) equals 0 exactly.
- [ ] 7f.4 **Signal-recovery sanity**: synthetic data with half
  null half slab `N(0, 1)`, `n = 1000` per scale, `s = 0.3` →
  `pi_0_hat ∈ [0.4, 0.6]`, `sigma_hat ∈ [0.85, 1.15]` (Normal
  path). Same data + Laplace path → `pi_0_hat ∈ [0.4, 0.7]`,
  Laplace `scale` consistent with the marginal moment.
- [ ] 7f.5 **Noise-recovery sanity**: synthetic data all from
  `N(0, shat^2)` (no signal) → `pi_0_hat >= 0.95` on both
  Normal and Laplace paths over `n = 1000` observations.
- [ ] 7f.6 **idx_size = 1 doesn't crash**: ebnm wrapper on a
  single observation returns a valid `fitted_g` (length-2
  `pi`, sum-to-1) on both Normal and Laplace. Estimates may
  be unstable; verify the *shape* contract only.

## 7g. End-to-end shape / contract

- [ ] 7g.1 `mfsusie(prior_variance_scope = "per_scale_normal", ...)`
  on a small fixture returns a fit with the documented shape:
  `class(fit$G_prior[[m]])` equals
  `"mixture_point_normal_per_scale"` for every m; each
  `fit$G_prior[[m]][[s]]$fitted_g` has `length(pi) == 2`,
  `sd == c(0, sigma)` for some `sigma >= 0`; `fit$pi_V[[m]]` is
  `S_m × 2`; `fit$alpha`, `fit$mu`, `fit$mu2`, `fit$pip`,
  `fit$sets$cs` populated.
- [ ] 7g.2 Same test for
  `prior_variance_scope = "per_scale_laplace"`:
  `class(fit$G_prior[[m]])` equals
  `"mixture_point_laplace_per_scale"`; `fitted_g` has the
  laplacemix shape (`pi`, `scale`, `mean` of length 2); fit
  object populated.
- [ ] 7g.3 The fit object inherits `c("mfsusie", "susie")`
  classes (no regression vs the current paths).
- [ ] 7g.4 `predict.mfsusie`, `coef.mfsusie`, `fitted.mfsusie`,
  `summary.mfsusie`, `print.mfsusie` all work on both new-path
  fits with no errors.

## 8. End-to-end power tests

- [ ] 8.1 **Sparse-wavelet recovery (synthetic).** The fixture
  from the proposal's acceptance criteria
  (`set.seed(2L); n = 200; p = 30; T_m = 64; signal at
  variables c(7, 18) and times c(20, 44)`):
  - `mfsusie(per_scale_normal)` recovers PIP on variables 7 and
    18 above 0.9; returns `nCS == 2`.
  - `mfsusie(per_scale_laplace)` recovers the same.
  - `mfsusie(per_outcome)` is the side-by-side reference and
    also recovers both signals.
- [ ] 8.2 On the why_functional simulation (clean Gaussian
  smooth effect, `n = 100`, `p ≈ 1000`, `T = 128`):
  `mfsusie(per_scale_normal)` and `mfsusie(per_scale_laplace)`
  each recover the planted causal in a 95% CS at PIP comparable
  to the `per_outcome` default (within 0.05 absolute).
- [ ] 8.3 On the susie-coloc multi-outcome fixture (planted
  shared causals at SNPs 2 and 7):
  `mfsusie(per_scale_normal)` and `mfsusie(per_scale_laplace)`
  each recover both CSes, diagonal-PP.H4 > 0.5, no spurious
  third CS.
- [ ] 8.4 Multi-outcome correctness: M=3, planted shared
  causal across all outcomes, both new paths recover the
  shared CS with `marginal_prob >= 0.5` for every outcome
  under
  `susie_post_outcome_configuration(by = "outcome",
  method = "susiex")`.
- [ ] 8.5 Null-locus stability: simulation with NO planted
  effect (all noise), both new paths return 0 CSes (no
  spurious findings). The point-* parametric regularization
  should make this clean.
- [ ] 8.6 Sparse-coverage stress: simulation with
  `wavelet_magnitude_cutoff > 0` (some near-zero columns get
  masked), both new paths produce a valid fit (no NaN, no
  crash) with the masked columns honored.

## 8b. Vignette sweep (acceptance gate)

After all unit tests pass, run a side-by-side sweep across every
shipped vignette except `_practical_data_applications.Rmd`
(hidden from the build, uses pinned fits). Verify both new
scope values agree with the default `per_outcome` on
credible-set recovery and run measurably faster.

- [ ] 8b.1 Extend
  `inst/bench/profiling/per_scale_normal_vignette_sweep.R` so
  that for each of the 9 user-facing vignettes
  (`fsusie_colocalization.Rmd`,
  `fsusie_covariates_adjustment.Rmd`,
  `fsusie_dnam_case_study.Rmd`, `fsusie_gtex_case_study.Rmd`,
  `fsusie_intro.Rmd`, `fsusie_why_functional.Rmd`,
  `getting_started.Rmd`, `mfsusie_intro.Rmd`,
  `mfsusie_long_running_fits.Rmd`, `post_processing.Rmd`),
  it reproduces the headline `fsusie()` / `mfsusie()` call from
  that vignette and re-runs it three times: at the default
  scope, at `per_scale_normal`, and at `per_scale_laplace`.
  Capture for each run: `system.time()`,
  `length(fit$sets$cs)`,
  `lapply(fit$sets$cs, function(cs) cs[which.max(fit$pip[cs])])`
  (lead variable per CS), and the full CS membership lists.
- [ ] 8b.2 For each workload, compute (a)
  `|nCS_new - nCS_po|`, (b) per-CS Jaccard between matched CSes
  (match by minimum lead-variable LD distance, or by index when
  CS counts agree), (c) lead-variable equality or "lead in the
  other CS", (d) wallclock ratio
  `t_per_outcome / t_<new_scope>`.
- [ ] 8b.3 Acceptance gate (per new scope, evaluated
  independently):
  - `|nCS_new - nCS_po| <= 1` on every workload.
  - Lead variable equal or LD-tag-shifted for every matched CS.
  - Per-CS Jaccard >= 0.5.
  - Speedup >= 2x on at least 4 of the 9 workloads (median
    speedup recorded as a reference number, not a hard
    threshold beyond the 2x-on-4 rule).
  Output the table to
  `inst/bench/profiling/per_scale_normal_vignette_sweep.tsv`
  with one row per (workload, scope) pair.
- [ ] 8b.4 Commit the sweep script and the output TSV. Note
  the commit hash so the sweep can be re-run on future changes
  for regression tracking.
- [ ] 8b.5 If any (workload, scope) pair fails the gate (CS
  count diverges, lead variable lands outside the matched CS,
  Jaccard < 0.5, or speedup is < 1x on a fixture small enough
  that the gain should be visible), block the merge and
  investigate. Likely causes: the marginal-data lead picker
  hits a noise variable on a low-signal scale; ebnm fits a
  near-null prior at init that the M-step cannot escape.

## 9. Documentation

- [ ] 9.1 Update `@param prior_variance_scope` roxygen in
  `R/mfsusie.R` with the per_scale_normal and per_scale_laplace
  descriptions.
- [ ] 9.2 Update `R/prior_scale_mixture.R` roxygen for the new
  `init_ebnm_prior_per_scale` helper and the dispatch.
- [ ] 9.3 Mention both new values in
  `vignettes/fsusie_intro.Rmd` and `vignettes/mfsusie_intro.Rmd`
  where the per_outcome / per_scale scope is currently
  discussed. Use **lead variant** in vignette prose (genetics
  audience); code and roxygen use **lead variable**.
- [ ] 9.4 Re-document via `devtools::document()`.

## 10. Final regression check

- [ ] 10.1 Run the full unit-test suite. All tests pass at the
  documented tolerances.
- [ ] 10.2 Both
  `mfsusie(prior_variance_scope = "per_scale_normal")` and
  `mfsusie(prior_variance_scope = "per_scale_laplace")` run to
  completion on a small fixture and return a fit with the
  documented shape.
