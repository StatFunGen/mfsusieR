# Test fidelity and manuscript cross-reference policy

Binding policy for Phase 3+ unit tests in mfsusieR. Companion to
`inst/notes/refactor-discipline.md` (which covers the porting and
Pattern A / Pattern B rules).

This file is the canonical source. Memory copies in
`~/.claude/projects/.../memory/` are short pointers at this file
so the policy is reproducible across sessions and contributors.

## Rule 1: Test fidelity

Every ported function, ELBO formula, posterior moment, log-Bayes
factor, residual variance update, prior variance update, PIP
formula, and CS construction step in mfsusieR SHALL have a unit
test that compares against the same quantity from the relevant
port source on a fixed seed. The three apple-to-apple equivalence
contracts (design.md D11a / D11b / D11c) drive the comparison
target:

- C1 (`M = 1`, `T_1 = 1`, scalar prior, etc.) -> `susieR::susie`
  at tolerance `<= 1e-10`.
- C2 (`M = 1`, `T_1 > 1`) -> `fsusieR::susiF` at tolerance
  `<= 1e-8`.
- C3 (`M >= 1`, legacy variance mode) ->
  `mvf.susie.alpha::multfsusie` at tolerance `<= 1e-8`.

"All details" means: same `alpha`, same `mu`, same `mu2`, same
`lbf`, same `lbf_variable`, same `KL` per effect, same `sigma2`,
same `elbo` trajectory, same `pip`, same `cs` membership at the
appropriate tolerance. Aggregate "fits roughly the same" checks
are not enough.

### Why

The port has to be numerically faithful before any redesign. If a
later FDR fix changes a number, we need to know whether the
change came from the fix or from a porting bug.

### How to apply

- Phase 4 tests follow the canonical shape in CLAUDE.md but
  tighten the `expect_equal` calls: every numeric output that
  exists in both packages gets a paired comparison, not just
  `pip` and `length(cs)`.
- Tolerance is binary per the refactor discipline
  (`inst/notes/refactor-discipline.md` section 3): `<= 1e-8` for
  apple-to-apple comparisons, `<= 1e-10` for deterministic
  analytic intermediates and the susieR-degeneracy test. Smoke
  tests (no numerical comparison) for apple-to-orange. Graduated
  tolerances (`1e-6`, `1e-2`, `5e-2`) are FORBIDDEN.
- The Pattern B `1e-12` tolerance for fsusieR-bug-fix unit tests
  is a documented exception; see
  `inst/notes/refactor-discipline.md` section 4.
- When mfsusieR deviates from the port source intentionally
  (e.g., the PIP-after-CS-filter fix, Pattern A in the
  refactor-discipline policy), the test asserts the deviation
  explicitly with a comment citing the OpenSpec change that
  authorized it.
- If a comparison fails and the discrepancy is small but nonzero,
  do NOT reframe it as "acceptable". Flag it to Gao per CLAUDE.md
  hard rule #6.

## Rule 2: Math cross-reference

Every formula in the codebase (docstrings, design docs) SHALL be
tagged with the corresponding manuscript reference. The
manuscript lives at `~/GIT/MultifSuSiE_Manuscript`. Key files:

- `methods/derivation.tex` - posterior derivations, KL, ELBO algebra
- `methods/algorithms.tex` - IBSS pseudocode, EM updates, convergence
- `methods/online_method.tex` - main methods section
- `results/1_methods_overview.tex` - figure 1 / overview equations

The recent commit history shows ongoing math revisions; treat the
current `main` branch of the manuscript as the working source of
truth.

### Why

The original code (`mvf.susie.alpha`, `fsusieR/susiF`) and the
manuscript's math are not in 1:1 correspondence; some
discrepancies are the FDR bug. Tagging every formula makes
mismatches visible at code-review time rather than at
simulation time.

### How to apply

- In roxygen and design.md, cite manuscript locations as
  `manuscript:methods/derivation.tex#sec:posterior` (line numbers
  are brittle; section labels or equation tags like
  `eq:elbo_mvf` are preferred).
- `@manuscript_ref methods/derivation.tex eq:foo` is required in
  main package roxygen. `@references_original
  mvf.susie.alpha/R/...` and bare `mvf.susie.alpha` / `fsusieR` /
  `susiF` / `multfsusie` strings are FORBIDDEN in main package
  roxygen per design.md D12 and
  `inst/notes/refactor-discipline.md` section 7. Original-code
  provenance lives in test file headers,
  `inst/notes/refactor-exceptions.md`, and dev notes only.
- Phase 1 paradigm notes do NOT need to be retro-fitted with
  manuscript refs. New writing from Phase 2 onward does.
- If the code's derivation differs from the manuscript, that is
  a finding, not a free implementation choice. Open an OpenSpec
  proposal with the diff written out before writing code.
- Do not edit the manuscript repo. It is read-only reference,
  like the other paradigm packages.
