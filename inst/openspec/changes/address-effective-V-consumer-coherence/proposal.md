# Address effective-V consumer coherence

## Why

Commit 605f395 ("Make V more informative; prepare code audit")
changed `post_loglik_prior_hook.mf_individual` to populate
`model$V[l]` with the per-effect effective slab variance
(mean over (m, s) of `sum_k pi[l, m, s, k] * var_k`) instead of
the placeholder value 1, and dropped V from the cleanup strip
list so it survives onto the final fit.

The 2026-04-30 cross-package audit (Agents A, B, C) flagged four
downstream consequences that need consistent handling:

- **A-1 (susieR / mfsusieR)**: `trim_null_effects.mf_individual`
  is still a documented no-op based on the obsolete "V[l] = 1
  always" assumption, but `susie_get_pip` now drops effects with
  `V[l] < prior_tol`. The two pruners are out of sync: `pip` no
  longer counts dropped effects but `alpha`/`mu`/`mu2`/`lbf`/`KL`
  for those same effects are kept on the fit.
- **B-2 (mfsusieR vs fsusieR)**: mfsusieR's PIP now V-gates
  collapsed effects automatically; fsusieR has no analog. Need
  to document the gate in the `mfsusie()` return-value docstring
  so users do not silently see different PIPs from the same data.
- **C-4 (mfsusieR vs mvf.susie.alpha)**: weak-signal regimes can
  push the effective V below `prior_tol = 1e-9`, silently
  dropping a real-but-weak effect from the PIP. Need a regression
  test on a controlled null fixture and an explicit decision on
  the `prior_tol` default.
- **C-5 (per-effect SER step)**: `model$V[l]` dual-encodes two
  semantically different things at different points in the IBSS
  step (1 during BF computation, effective-V between iters). The
  pre-hook silently discards `V_init`; this is correct but
  non-obvious. Needs a docstring note.

These are all consequences of the same design choice; addressing
them in one change keeps the V-semantics narrative coherent.

## What changes

### 1. `trim_null_effects.mf_individual` — coherent with the V filter

Either (a) delete the override so it falls through to
`trim_null_effects.default` (which zeroes alpha/mu/mu2/lbf/KL for
`V[l] < prior_tol`, matching the PIP filter), or (b) keep the
override but make it perform the same zeroing on the per-effect
list-of-list `mu`/`mu2`. Pick (a) if the default works on
mfsusieR's mu shape; (b) otherwise.

### 2. Document the V-filter in the public API

Add a paragraph to the `mfsusie()` return-value roxygen for
`pip`, naming the V-based filter and the `prior_tol` argument.
Also note in NEWS.md.

### 3. Regression test on null-only fixture

Fit a synthetic null-only fixture (no causal signal) at default
settings and assert `model$pip` matches what
`susie_get_pip(model, prior_tol = 0)` returns when the V filter
is disabled. The test exists to flag (and accept or revise) any
silent-drop behavior on weak-signal cases.

### 4. Pre-hook docstring note

Add a one-line comment to `pre_loglik_prior_hook.mf_individual`
explaining that `V_init` is intentionally discarded — V=1 is the
load-bearing semantic for `loglik.mf_individual`'s mixture sd
grid; the post-hook's effective V lives only on the model object
between IBSS iters.

## Impact

- **Severity**: medium. Coherence fix; not blocking release but
  affects user-visible PIP semantics.
- **Source-code changes**: `R/individual_data_methods.R`
  (trim_null_effects override and pre-hook docstring),
  `R/mfsusie.R` (return-value roxygen).
- **Tests**: one new regression test on the null-only fixture.
- **Documentation**: NEWS.md entry for the V-filter exposure.

## Out of scope

- Performance tuning of the V-filter threshold itself.
- Changes to the effective-V computation (already locked in
  605f395).
- Phase 5 FDR comparison of mfsusieR's V-filtered PIP vs the
  pre-605f395 unfiltered PIP — that is a Phase 5 investigation,
  not part of this change.

## References

- Audit Agent A report:
  `inst/notes/sessions/2026-04-30-cross-package-audit-a-susier-posthooks.md`
  (Finding A-1)
- Audit Agent B report:
  `inst/notes/sessions/2026-04-30-cross-package-audit-b-fsusier-posthooks.md`
  (Finding B-2)
- Audit Agent C report:
  `inst/notes/sessions/2026-04-30-cross-package-audit-c-mvf-posthooks.md`
  (Findings C-4, C-5)
- Audit summary:
  `inst/notes/cross-package-audit-summary-posthooks.md`
- Source change:
  `inst/openspec/changes/audit-cross-package-post-hooks/`
