# Tasks

Depends on: rounds 3a (`add-preprocessing-low-count-and-quantile-norm`,
`expose-model-init-warm-start`, `expose-nig-residual-variance-method`),
3b (`add-practical-data-applications-vignette`,
`fix-scalewise-sd-and-add-smash-method`,
`broaden-reference-tests-coverage`) all archived.

## 1. Launch three opus-model agents in parallel

- [ ] 1.1 Agent A (mfsusieR <-> susieR). Prompt: read every
      public-API function in mfsusieR side-by-side with the
      corresponding susieR function. Flag (a) numerical
      divergence not documented in refactor-exceptions and
      (b) helpers duplicated when we should be delegating.
      Save report to
      `inst/notes/sessions/<date>-cross-package-audit-a-susier.md`.
- [ ] 1.2 Agent B (mfsusieR <-> fsusieR). Same task, focused
      on the single-modality functional case.
      Report at `<date>-cross-package-audit-b-fsusier.md`.
- [ ] 1.3 Agent C (mfsusieR <-> mvf.susie.alpha). Same task,
      multi-modality case.
      Report at `<date>-cross-package-audit-c-mvf.md`.

## 2. Triage findings

- [ ] 2.1 Read all three reports. For each finding, classify
      as `accept-as-documented`, `fix-now`, or
      `track-for-later`.
- [ ] 2.2 For `fix-now` findings, open a targeted OpenSpec
      change. For `track-for-later`, open an OpenSpec
      change with task list deferred to a future round.

## 3. Summary

- [ ] 3.1 Write `inst/notes/cross-package-audit-summary.md`
      listing the three reports, the triage classification
      per finding, and the names of any follow-up OpenSpec
      changes opened.

## 4. Archive

- [ ] 4.1 `openspec archive multi-agent-cross-package-audit`
      after all `fix-now` findings are addressed.
