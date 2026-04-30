# Tasks

Severity: HIGH. Block feature work in the SER scaffolding,
prior storage, and convergence path until Phase B is
complete.

## 1. Phase A. Three Opus agents in parallel

Each agent runs from a fresh context. Each agent prompt
must include:

1. The "Known divergences (do NOT re-flag)" section of
   `proposal.md` verbatim.
2. The full text of `inst/notes/refactor-exceptions.md`.
3. The recent commit list from proposal.md (so the agent
   knows what changed since 2026-04-26).
4. The agent-specific scope (which upstream package, which
   pathway).

- [ ] 1.1 Spawn Agent A (mfsusieR vs susieR). Scope =
      Agent A section of proposal.md. Report saved to
      `inst/notes/sessions/<date>-cross-package-audit-a-susier-posthooks.md`.
- [ ] 1.2 Spawn Agent B (mfsusieR vs fsusieR). Scope =
      Agent B section of proposal.md. Report at
      `<date>-cross-package-audit-b-fsusier-posthooks.md`.
- [ ] 1.3 Spawn Agent C (mfsusieR vs mvf.susie.alpha).
      Scope = Agent C section. Report at
      `<date>-cross-package-audit-c-mvf-posthooks.md`.

## 2. Phase B. Triage

- [ ] 2.1 Read all three reports. For each finding,
      classify as `accept-known`, `add-to-divergence-ledger`,
      `fix-now`, or `track-for-later`.
- [ ] 2.2 For `add-to-divergence-ledger`: add the entry to
      `~/Documents/mfsusier-mvf-divergences/divergences.tex`
      and a stanza to `inst/notes/refactor-exceptions.md`
      in the SAME session.
- [ ] 2.3 For `fix-now`: open a targeted OpenSpec change.
      Link from this proposal.
- [ ] 2.4 For `track-for-later`: open an OpenSpec change
      with a deferred task list.

## 3. Phase C. Summary

- [ ] 3.1 Write `inst/notes/cross-package-audit-summary-posthooks.md`
      collating the three reports, the triage classification
      per finding, and the names of any follow-up OpenSpec
      changes.
- [ ] 3.2 Append to `inst/notes/refactor-exceptions.md` a
      "Cross-package audit (post hooks, <date>)" subsection
      pointing to the summary.

## 4. Process gate

- [ ] 4.1 After this audit archives, every future PR that
      touches `R/individual_data_methods.R`,
      `R/model_class.R`, `R/ibss_methods.R`, or
      `R/em_helpers.R` requires a cross-package audit-
      style review before merge. Document the gate in
      `inst/notes/review-loop-methodology.md`.
