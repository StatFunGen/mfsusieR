# Cross-package audit capability spec delta

## ADDED Requirements

### Requirement: post-hook cross-package audit cycle

mfsusieR SHALL run a fresh-context multi-agent audit
comparing the current state to susieR, fsusieR, and
mvf.susie.alpha after any substantial refactor of the SER
scaffolding, prior storage, M-step dispatch, or
convergence path. The audit MUST produce one report per
upstream package plus a triage summary, before the next
release.

#### Scenario: change touches the SER hot path

- **WHEN** a PR modifies `R/individual_data_methods.R`,
  `R/model_class.R`, `R/ibss_methods.R`, or
  `R/em_helpers.R`
- **THEN** the PR description SHALL link an open
  cross-package audit OpenSpec change OR flag the PR as
  `audit-deferred` with a written justification

#### Scenario: agent receives the known-divergence ledger

- **WHEN** an audit agent is launched
- **THEN** its prompt SHALL include verbatim
  (i) the ledger at
  `~/Documents/mfsusier-mvf-divergences/divergences.tex`
  (current entries), and
  (ii) the full `inst/notes/refactor-exceptions.md`,
  AND the agent SHALL NOT flag any finding that matches
  an existing ledger entry.
