# Tasks

## 0. Decisions (user)

- [ ] 0.1 Approve this proposal.
- [ ] 0.2 Part B: B2 as default, or behind an option first (design.md,
  open question 1).
- [ ] 0.3 Withdraw `perf-skip-elbo-on-pip` (move to `archive/` with a
  "withdrawn: premise invalid" note).

## 1. Part A: remove dead objective code (fitting behavior unchanged)

- [ ] 1.1 Save baseline fits: fixed-seed fixtures, four
  `prior_variance_scope` values x `convergence_method` in
  {`"pip"`, `"elbo"`} x `L_greedy` in {default, `NULL`}.
- [ ] 1.2 Delete `get_objective.mfsusie` and
  `refresh_lbf_kl.mf_individual` (`R/individual_data_methods.R`).
- [ ] 1.3 Remove the `update_alpha` argument and branch from
  `loglik.mf_individual`.
- [ ] 1.4 Remove `"get_objective"` from `mfsusie_generics`
  (`R/zzz.R`); update the comments there and at `R/mfsusie.R:487`.
- [ ] 1.5 Test: `getS3method("get_objective", "mfsusie", optional =
  TRUE, envir = asNamespace("susieR"))` is NULL, and the fit's ELBO
  equals `get_objective.default` at the final state.
- [ ] 1.6 Verify `identical()` against 1.1 on every fixture; full
  test suite shows no new failures.
- [ ] 1.7 `devtools::document()`; commit
  `refactor(objective): remove unreachable get_objective.mfsusie and refresh`.

## 2. Part B: final E-step after the last inner M-step (B2)

- [ ] 2.1 In `post_loglik_prior_hook.mf_individual`, after the M-step
  on the cap exit, run `loglik -> calculate_posterior_moments ->
  compute_kl` against the final prior (also when
  `max_inner_em_steps = 0`).
- [ ] 2.2 Test: stored `KL[l]` matches the SER identity recomputed from
  the stored state and returned prior (`<= 1e-8`).
- [ ] 2.3 Test: `fitted[[m]] == X %*% sum_l alpha_l mu_l[[m]]` after
  `trim_null_effects` (`<= 1e-8`), on a fixture that hits the cap.
- [ ] 2.4 Test: `convergence_method = "elbo"` ELBO non-decreasing on
  the smoke fixtures, every scope.
- [ ] 2.5 Runtime overhead on the `optimization/scripts/e2e_compare.R`
  fixtures.
- [ ] 2.6 FDR null simulations (existing harness), before vs after;
  report to the user before the default changes.
- [ ] 2.7 `NEWS.md` entry; commit `fix(inner_em): end inner EM with an E-step`.

## 3. Close

- [ ] 3.1 Update `inst/notes/refactor-exceptions.md` if B2 deviates
  from `mvf.susie.alpha`'s inner-EM order.
- [ ] 3.2 Archive this change.
