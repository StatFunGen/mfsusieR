# Plan handoff: ebnm-based per-scale priors + dispatch fix

This is a handoff for the next session. It replaces the
"alpha-weighted MLE" implementation in
`inst/openspec/changes/add-per-scale-point-normal-prior/` with an
ebnm-based design. The OpenSpec change should be amended in place;
do not open a new change unless Gao asks for one.

## 1. Why this exists (in two paragraphs)

The 2026-04-29 session implemented `prior_variance_scope = "per_scale_normal"` with a custom L-BFGS-B MLE
(`mf_em_point_normal`) on alpha-weighted observations. The vignette
sweep at `inst/bench/profiling/per_scale_normal_vignette_sweep.R`
showed the alpha-weighted approach produced `nCS = 0` and uniform
`PIP = 1/p` on every realistic vignette, because diffuse alpha
collapses sigma to ~1e-5 (slab and spike both Diracs at zero).
A pivot to lead-variable-only weighting recovered signal on three
of four bundled-data workloads, but exposed two further bugs.

**Bug 1 (uniform-alpha iter-1 trap):** at IBSS iter 1
`alpha[l, ]` is `1/p` uniform; `which.max(zeta_keep)` returns
variable 1 (deterministic tie-break) and the M-step fits a noise
column. The resulting degenerate prior produces uniform alpha
again -- a fixed point at "no signal." Trace at
`/tmp/trace_psn.R` confirms iters 1-5 stuck at uniform alpha,
lead = 1, on the synthetic sparse-fine-scale fixture (`set.seed(2)`,
`n=200, p=30, T=64`, signal at `var c(7, 18)` and `T c(20, 44)`).

**Bug 2 (sparse-signal sigma underestimation in init):**
`init_point_normal_prior_per_scale` uses
`sqrt(mean(Bhat^2) - mean(Shat^2))`. For a wavelet-domain spike
where 1-2 of `idx_size` positions per scale carry signal, the mean
averages signal in with the bulk of noise -- sigma_init ~ 0.03 vs
the true signal scale ~ 0.4 (factor of 10x underestimate).
Even with bug 1 fixed, the SER under this weak prior cannot
amplify the planted causal enough to concentrate alpha onto it.

## 2. The agreed pivot

Three pieces, addressed together:

1. **ebnm replaces the custom M-step solver.** `mf_em_point_normal()`
   goes away. The M-step calls
   `ebnm::ebnm_point_normal(x, s)` for `per_scale_normal` and
   `ebnm::ebnm_point_laplace(x, s)` for `per_scale_laplace`.
   ebnm handles its own init heuristics, optim convergence,
   and boundary cases (the `qlogis(0) = -Inf` clamp, the
   `n <= 2` freeze, the `factr` choice, the optim convergence
   warnings -- all become ebnm's problem).

2. **Option A's dispatch fix for the cold start.** Pick the
   lead variable per scale at *init* from the marginal data
   (`argmax_j mean_t(Bhat[j, idx_s]^2)`). This is a well-defined
   choice (no ties, no `alpha`-dependent threshold) and gives a
   sigma matched to the strongest-looking variable at that
   scale. After IBSS effect 1 has converged once, the M-step
   uses the conventional `argmax(alpha[l, ])` lead.

3. **Add `per_scale_laplace` as a parallel prior class.**
   The Laplace slab has heavier tails than Normal; for sparse
   wavelet-domain signals (1-2 active positions per scale) the
   tail accommodates outliers without averaging them out.
   Users pick by `prior_variance_scope` enum extension:
   `c("per_outcome", "per_scale", "per_scale_normal",
   "per_scale_laplace")`.

## 3. State of the tree at session start

- Branch: `main` (a previous session committed Steps 1-3 of the
  custom-EM design plus 1543 unit tests passing). Verify with
  `git log --oneline -5`.
- `R/em_helpers.R` has `mf_em_point_normal()` -- to remove.
- `R/individual_data_methods.R` has `.opv_point_normal()` with
  the `2.0/p` deferred-M-step threshold and the lead-variable
  dispatch -- to be replaced.
- `R/prior_scale_mixture.R` has
  `init_point_normal_prior_per_scale()` with the per-variable-mean
  init heuristic (mid-experiment state, not a defended choice) --
  to be replaced.
- `tests/testthat/test_per_scale_normal_degeneracy.R` has 213
  passing tests across sections 7e/7f/7g/7b-d/8. Most still
  apply; the M-step solver tests (7e.1-7e.4) need rewriting
  for the ebnm wrapper.
- `inst/bench/profiling/per_scale_normal_vignette_sweep.R` has
  the 4-workload sweep harness. Reusable.
- OpenSpec change at
  `inst/openspec/changes/add-per-scale-point-normal-prior/` --
  amend in place (proposal.md, design.md, tasks.md, spec.md).

## 4. Implementation steps

### Step 0 -- amend the OpenSpec change (no code yet)

- Edit `proposal.md`: replace the "alpha-weighted MLE" rationale
  with the ebnm + lead-from-marginal-data design. Note that
  point-laplace is added alongside point-normal.
- Edit `design.md`: replace section "Alpha-weighted MLE (no
  lead-variant approximation)" with "ebnm-driven M-step +
  cold-start lead from marginal data."
- Edit `tasks.md`: rewrite section 5 (init helper + custom EM)
  and section 6 (M-step dispatch). Add section 5c for
  `per_scale_laplace`.
- Edit `specs/scale-mixture-prior/spec.md`: ADD requirement that
  the cold-start lead is picked from marginal data; ADD
  requirement that ebnm handles the M-step.
- Run `cd inst && openspec validate` to confirm.

### Step 1 -- add ebnm to DESCRIPTION

- `Imports: ebnm` (or `Suggests:` if we want it optional with a
  graceful fallback to a stub solver -- discuss with Gao before
  picking).
- `R CMD check` on a clean tree.

### Step 2 -- replace `init_point_normal_prior_per_scale` (Option A)

New helper signature:

```
init_marginal_lead_per_scale(Y_m, X, groups, prior_class)
```

For each scale s with index set `idx_s`:
1. Compute `bs <- compute_marginal_bhat_shat(X, Y_m)`.
2. `lead_s <- which.max(rowMeans(bs$Bhat[, idx_s, drop = FALSE]^2))`.
3. Run `ebnm::ebnm_point_normal(bs$Bhat[lead_s, idx_s],
    bs$Shat[lead_s, idx_s])` (or `ebnm_point_laplace` if
    `prior_class == "mixture_point_laplace_per_scale"`).
4. Store `fitted_g$pi`, `fitted_g$sd`, `fitted_g$mean` (and
    `fitted_g$scale` for laplace) in the `G_prior[[s]]$fitted_g`
    slot. Record `lead_init_s = lead_s` on the entry for tests.
5. Set `class(G_prior) <- prior_class`.

This eliminates the moment-estimator init entirely; ebnm fits
the prior at init from a real signal-tracking observation set.

### Step 3 -- replace `.opv_point_normal` with ebnm dispatch

Two prior classes share the same dispatch shape, so route through
one helper:

```
.opv_ebnm <- function(data, params, model, ser_stats,
                      keep_idx, zeta_keep, ebnm_fn) {
  lead <- keep_idx[which.max(zeta_keep)]
  for (m in seq_len(data$M)) {
    bhat_m <- ser_stats$betahat[[m]]; shat_m <- sqrt(ser_stats$shat2[[m]])
    G_m <- model$G_prior[[m]]
    for (s in seq_along(G_m)) {
      idx <- G_m[[s]]$idx
      fit <- ebnm_fn(x = bhat_m[lead, idx], s = shat_m[lead, idx])
      model$G_prior[[m]][[s]]$fitted_g <- fit$fitted_g
    }
  }
  model
}
```

Then dispatch from `optimize_prior_variance.mf_individual` on
`inherits(model$G_prior[[1L]], "mixture_point_normal_per_scale")`
or `"mixture_point_laplace_per_scale"`, with `ebnm_fn` set to
`ebnm::ebnm_point_normal` or `ebnm::ebnm_point_laplace`
respectively. Drop the `2.0/p` deferred-M-step guard --
the cold-start init from marginal data already handles the
chicken-and-egg, so the M-step at iter 2+ runs unconditionally.

### Step 4 -- delete obsolete code

- `mf_em_point_normal()` in `R/em_helpers.R` -- delete.
- The `qlogis(min(max(pi_0_init, 1e-3), 1 - 1e-3))` clamp -- gone
  with the function.
- Any `sigma_init` / `pi_0_init` parameters threaded through
  the dispatch.
- The `params$mixsqp_alpha_eps` thinning at the M-step: ebnm
  does its own data-quality filtering. Keep `mixsqp_alpha_eps`
  on the mixsqp branch.

### Step 5 -- update tests

- 7e.1-7e.4 (M-step solver fidelity): rewrite as wrapper tests
  that confirm we forward `(x, s)` to ebnm correctly and unpack
  `fitted_g` correctly. The "idempotence at 1e-14" claim no
  longer applies (ebnm is the authority).
- 7f.1-7f.8 (numerical properties on IBSS): mostly still hold;
  re-run and adjust tolerances if needed.
- 7g.1-7g.3 (end-to-end shape): add `per_scale_laplace` cases
  alongside `per_scale_normal`.
- 7b/c/d (cross-flavor degenerate equivalence): ebnm with a
  fixed prior (length-1 grid + no M-step) collapses to the same
  Gaussian as `per_outcome` and `per_scale`. Bit-match contract
  should still hold; verify.
- 8.1-8.5 (end-to-end power): rerun. The signal-recovery
  contract on `scenario_minimal` should now pass for both
  `per_scale_normal` and `per_scale_laplace`.

### Step 6 -- vignette sweep acceptance gate

Run `Rscript inst/bench/profiling/per_scale_normal_vignette_sweep.R`
twice -- once with `per_scale_normal`, once with
`per_scale_laplace`. Acceptance criteria from
`tasks.md` 8b.3:

- `|nCS_psn - nCS_po| <= 1` on every workload.
- Lead variable equality on matched CSes.
- Mean Jaccard >= 0.5 per workload.
- Speedup >= 2x on >= 4 of the 9 workloads (the script ships 4;
  add the remaining 5 vignette workloads if Gao approves).

### Step 7 -- documentation

- `R/mfsusie.R` `@param prior_variance_scope`: list four values.
- `R/prior_scale_mixture.R` `@param prior_variance_scope`: same.
- `vignettes/fsusie_intro.Rmd` and `vignettes/mfsusie_intro.Rmd`:
  add `per_scale_laplace` to the prior-choice prose.

  Note: vignette prose uses **lead variant**, not "lead variable"
  (per `~/.claude/projects/-home-gw-GIT-mfsusieR/memory/feedback_lead_variable_terminology.md`).
  Code, comments, roxygen use **lead variable**.

## 5. Stress-test fixture for the sparse-signal failure mode

The synthetic fixture that originally exposed bug 2 (kept here so
the next session can reproduce):

```
set.seed(2L)
n <- 200L; p <- 30L; T_m <- 64L
X <- matrix(rnorm(n * p), n)
shape <- numeric(T_m); shape[c(20L, 44L)] <- 1
beta  <- numeric(p);    beta[c(7L, 18L)]  <- c(1.0, -1.0)
Y <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
     matrix(rnorm(n * T_m, sd = 0.5), n)
```

`per_outcome` recovers PIP `[7] = 1.000`, `[18] = 1.000` in 2 CSes.
`per_scale_normal` (current alpha-weighted then lead-variable
attempts) gets PIP `[7] ≈ 0.41-0.74`, `[18] ≈ 0.49-0.85`,
returning 0 CSes. Add this as a unit test in section 8.

## 6. Risks Gao should weigh before approving

- **New runtime dependency.** ebnm pulls in `RcppArmadillo`, `mixsqp`,
  `ashr`. mfsusieR already depends on `susieR` which depends on
  `mixsqp`, so the marginal weight is mostly `ebnm` itself + `ashr`.
- **ebnm version pinning.** The point-laplace path is mature but
  the API is occasionally tweaked. Pin a known-good version in
  the C2/C3-style fixture machinery.
- **Bug 1 mitigation depends on Step 2 being correct.** If
  `argmax(rowMeans(Bhat^2))` happens to pick a noise variable on
  some pathological fixture, we are back to bug 1. Add a unit
  test that asserts the marginal-data lead picker selects a
  signal-bearing variable on `scenario_minimal` and on the
  synthetic sparse fixture.

## 7. How to start the new session

1. Open a fresh Claude Code session in `~/GIT/mfsusieR`. Use
   `/clear` if you want to keep the current shell, or close and
   reopen the terminal.
2. The first prompt to Claude:

   > Read `inst/notes/sessions/2026-04-29-ebnm-per-scale-pivot.md`
   > and execute the plan starting at Step 0. Stop after each
   > step's tests pass and report the diff before moving on.

3. Claude will read the OpenSpec change, this note, the relevant
   `R/` files, and start with the proposal/design/tasks
   amendments. Verify the plan-doc summary matches your
   intent before letting it touch code.
4. Each step in the plan ends with a verifiable artefact (passing
   tests, openspec validate, sweep TSV). Use those as natural
   checkpoints.

## 8. Memory entries the new session will already see

- `feedback_propose_before_commit.md` -- numerical fixes follow
  try-then-summarize; magic thresholds need discussion before
  committing.
- `feedback_lead_variable_terminology.md` -- code/roxygen use
  "lead variable"; vignettes use "lead variant".
- `refactor_skill_principles.md`, `review_loop_methodology.md`,
  `test_fidelity_and_manuscript_xref.md` -- standard mfsusieR
  conventions.
