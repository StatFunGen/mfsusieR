# mf-cross-package-audit capability

## ADDED Requirements

### Requirement: three independent audits per package family

A multi-agent cross-package audit SHALL run after every round of changes that touches numerical paths.

Three opus-model agents SHALL run in parallel:

- Agent A: mfsusieR side-by-side with susieR.
- Agent B: mfsusieR side-by-side with fsusieR.
- Agent C: mfsusieR side-by-side with mvf.susie.alpha.

Each agent SHALL flag:

- numerical divergence not already documented in
  `inst/notes/refactor-exceptions.md`, and
- code that duplicates an upstream helper that mfsusieR
  should be delegating to instead.

#### Scenario: per-package report saved

Each agent's findings SHALL be saved to a session note at
`inst/notes/sessions/<date>-cross-package-audit-{a,b,c}-<pkg>.md`
with file path and line-range citations for every finding.

### Requirement: findings are triaged into three buckets

Each finding SHALL be classified as one of:

- `accept-as-documented`: the deviation is recorded in
  refactor-exceptions and the audit confirms the rationale.
- `fix-now`: the deviation is a bug; a targeted OpenSpec
  change is opened in the same session.
- `track-for-later`: the deviation is bounded and not blocking;
  a follow-up OpenSpec change is opened with a deferred
  task list.

#### Scenario: triage summary present

A single `inst/notes/cross-package-audit-summary.md` SHALL
list the three reports, the triage classification per
finding, and the names of any follow-up OpenSpec changes
opened.
