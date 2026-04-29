# Tasks

## 1. Default change

- [ ] 1.1 In `R/mfsusie.R` formals, change
  `mixsqp_alpha_eps = 1e-6` to `mixsqp_alpha_eps = 5e-5`.
- [ ] 1.2 Update the `@param mixsqp_alpha_eps` roxygen to state
  the new default and re-state the truncation-error bound
  (`relative_error ≈ eps × O(10)`, ⇒ ~5e-4 with the new
  default — still 200× below the IBSS `tol = 1e-4` default).

## 2. Numerical-fidelity check vs old default

- [ ] 2.1 Bench script
  `inst/bench/profiling/alpha-eps-fidelity.R` that runs `mfsusie()`
  on three fixtures (`small`, `medium`, `large` p) under
  both `mixsqp_alpha_eps = 1e-6` and `mixsqp_alpha_eps = 5e-5`,
  records `alpha`, `pip`, `lbf`, `niter`, `sets$cs`, and
  outputs:
  - `max|alpha_new - alpha_old|` ≤ `1e-3`.
  - `max|pip_new - pip_old|` ≤ `1e-3`.
  - `niter` and `sets$cs` membership identical.
- [ ] 2.2 Run the script; commit the output to
  `inst/notes/sessions/<date>-alpha-eps-bump.md`.

## 3. Speedup check

- [ ] 3.1 Bench script
  `inst/bench/profiling/alpha-eps-speed.R` on `p = 5000`,
  `T = 128`, `M = 1` fixture; report total `mfsusie()` runtime
  and per-iter M-step time at both eps values.
- [ ] 3.2 Acceptance: ≥ 30% reduction in M-step time at the new
  default. If not, abandon the change.

## 4. Test suite

- [ ] 4.1 Run the full unit-test suite at the new default.
  All existing tests pass at their tolerances. (Tests that
  pin alpha values to high precision should not exist; if
  they do, document the impact and decide case-by-case.)
- [ ] 4.2 Verify `mixsqp_alpha_eps = 0` still recovers the
  full-set behavior (regression test for the escape hatch).
