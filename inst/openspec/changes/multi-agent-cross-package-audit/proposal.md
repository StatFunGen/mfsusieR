# Multi-agent cross-package correctness and duplication audit

## Why

After the round-3 changes land (preprocessing, model_init,
NIG, practical-data vignette, scalewise SD fix, smash, and
the coverage broadening), three independent fresh-context
opus-model agents read every public-API function in mfsusieR
side-by-side with the corresponding function in `susieR`,
`fsusieR`, and `mvf.susie.alpha`. Each agent flags two
categories:

1. Numerical divergence not documented in
   `inst/notes/refactor-exceptions.md`.
2. Code that duplicates an upstream helper that mfsusieR
   should be delegating to instead (most relevant for
   susieR; less for the port-source packages).

The audit catches drift introduced by the round-3 changes
(cumulative algorithmic deviations not flagged at PR time)
and identifies any helpers we re-implemented when we should
have delegated.

## What changes

### 1. Audit step

Three opus agents are launched in parallel:

- Agent A: `mfsusieR` <-> `susieR`. Focus on
  delegation: every numerical helper in mfsusieR that has
  a susieR equivalent should be delegated, not duplicated.
  Flag duplicates plus any divergence.
- Agent B: `mfsusieR` <-> `fsusieR`. Focus on the
  single-modality functional case (M=1, T_1>1). Every
  numerical output should match upstream at the
  established tolerance, modulo documented exceptions.
- Agent C: `mfsusieR` <-> `mvf.susie.alpha`. Focus on
  the multi-modality functional case. Every numerical
  output should match upstream at the C3 tolerance modulo
  documented exceptions.

Each agent returns a structured report listing findings
with file path and line range citations. The reports are
saved at `inst/notes/sessions/<date>-cross-package-audit-{a,b,c}.md`.

### 2. Triage

Findings are triaged into one of three buckets:

- `accept-as-documented`: deviation is in
  refactor-exceptions, no action.
- `fix-now`: clear bug, opens a targeted change in the
  same session.
- `track-for-later`: deviation is real but bounded; opens
  a follow-up OpenSpec change post-audit.

### 3. Audit summary

A single `inst/notes/cross-package-audit-summary.md`
collates the three reports plus the triage decisions and
any follow-up OpenSpec change names.

## Impact

- No source-code changes from this proposal directly.
- Findings drive zero or more follow-up changes.
- Documentation: three session notes plus a summary.

## Out of scope

- The audit runs after rounds 3-6 land. Earlier audits are
  out of scope for this proposal.
- Implementing fixes for findings is out of scope here;
  fixes are tracked in their own OpenSpec changes.
