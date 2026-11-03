# Fix objective dispatch dead code and make the reported ELBO coherent

## Why

The 2026-10-23 review (Claude + Codex, two rounds; record in
`optimization/bug1_discussion.md`) found that `get_objective.mfsusie`
never runs:

- susieR's workhorse calls `get_objective(data, params, model)`
  (`susieR/R/susie_workhorse.R:87`), and the generic dispatches on
  `data`, whose class is `c("mf_individual", "individual")`. mfsusieR
  registers the override on the model class `"mfsusie"`
  (`R/zzz.R:86-108`).
- So neither the `pip` short-circuit (commit dc5515e, change
  `perf-skip-elbo-on-pip`) nor `refresh_lbf_kl.mf_individual` (commit
  13b33f2) has ever executed. Under the default `pip` convergence,
  `fit$elbo` is computed every iteration by susieR's
  `get_objective.default` (`Eloglik - sum(KL)`). A repository search
  finds no other caller of either function.

Three facts shape the fix:

1. **The default path is the intended one.** Commit a7ba30a removed a
   trivial `get_objective.mf_individual` on purpose, and
   `tests/testthat/test_variance_and_elbo.R:168-198` asserts the
   default path.
2. **The refresh would be wrong if enabled.** It recomputes `KL_l` with a
   new residual while freezing `q_l`. The identity
   `KL = E_q log p(R|b) - log Z(R, g)` holds only when `q_l` is the exact
   posterior for `(R, g)`. Codex's counterexample: true KL 0.0966,
   refresh -0.903. The refresh also reads the shared `G_prior`
   scratchpad instead of `fitted_g_per_effect[[l]]`, and it mutates only
   a local copy of the model.
3. **Skipping the ELBO under `pip` is not worth it.**
   `optimization/scripts/elbo_cost.R` measures `get_objective` at 0.3%
   to 0.8% of `mfsusie()` time today, and 0.5% to 1.6% after the planned
   kernel/M-step speedups (n = 1000 to 5000, p = 1000 to 2000). Two
   vignettes, `print`/`summary`, and `test_mfsusie_smoke.R` use
   `fit$elbo` under `pip`. The premise of `perf-skip-elbo-on-pip` (the
   refresh doubles per-iteration work) is false, because the refresh
   never ran.

A real coherence gap remains. `post_loglik_prior_hook.mf_individual`
ends the inner EM with an M-step and no matching E-step
(`R/individual_data_methods.R:704-709`). So `KL[l]` is evaluated
against the prior in force before that M-step, while `G_prior`,
`fitted_g_per_effect`, `pi_V` and `V[l]` hold the new prior. The
reported value is `F(q, h)`, a valid lower bound under the per-effect
priors `h` used at each effect's last E-step, but not the free energy of
the returned state `(q, g)`. Codex notes that this bookkeeping alone does
not break monotonicity.

## What changes

**Part A: remove dead code. Fitting behavior unchanged.**
- Delete `get_objective.mfsusie`, `refresh_lbf_kl.mf_individual`, and the
  `update_alpha` argument of `loglik.mf_individual`.
- Remove `"get_objective"` from `mfsusie_generics` in `R/zzz.R`.
- Fix the stale comments (`R/zzz.R:90-95`, `R/mfsusie.R:487`).
- `fit$elbo` keeps its current values under both convergence methods.

**Part B: make the reported objective describe the returned state.**
Option B2 (recommended; see `design.md`): end the inner EM with an
E-step, so the last action per effect is
`loglik -> calculate_posterior_moments -> compute_kl` against the final
prior. `alpha`, `mu`, `mu2`, `KL`, `lbf` and `fitted` are then all
consistent with the returned prior. This changes fitted values
numerically, so it is gated by the acceptance criteria below, including
the FDR null simulations.

**Supersedes** `perf-skip-elbo-on-pip`. Proposed disposition: withdraw
that change (move it to `archive/` with a "withdrawn: premise invalid"
note). Do not implement it.

## Acceptance criteria

Part A:
- `identical()` holds on complete fit objects before and after, on
  fixed-seed fixtures for all four `prior_variance_scope` values, under
  both `convergence_method` values, with and without `L_greedy`.
- Full test suite: no new failures relative to the pre-change tree.

Part B:
- With `convergence_method = "elbo"`, `fit$elbo` is non-decreasing
  (`diff(elbo)[-1] >= -1e-8 * abs(elbo)`) on the smoke fixtures for
  every scope.
- For every effect, the stored `KL[l]`, recomputed from the stored
  `alpha`, `mu`, `mu2` and returned prior via the SER identity, matches
  at `<= 1e-8`.
- `fitted[[m]]` equals `X %*% sum_l alpha_l * mu_l[[m]]` at `<= 1e-8`
  after `trim_null_effects` (covers BUG-13).
- Runtime overhead at most one extra loglik + moments + KL pass per
  effect per iteration, reported on the
  `optimization/scripts/e2e_compare.R` fixtures.
- FDR null simulations (the existing `per_scale_normal` and
  `per_outcome` harness): calibration no worse than the pre-change
  baseline. Any change in power or FDR is reported to the user before
  this becomes the default.

## Impact

- Part A: internal only; no public API change.
- Part B: changes posterior output numerically (alpha, PIP, CS can
  move). It needs a `NEWS.md` entry and a decision on whether it ships
  as default or behind an option first.

## Out of scope

- Skipping ELBO computation under `pip` (dropped; see Why, point 3).
- The `log(pi + sqrt(eps))` prior-weight offset (same expression in
  susieR, `susieR/R/susie_utils.R:1197`).
- Changing susieR's workhorse or generic signatures.
