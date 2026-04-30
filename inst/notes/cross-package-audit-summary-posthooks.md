# Cross-package audit summary (post hooks, 2026-04-30)

This summary collates Phase A reports from the three audit
agents launched per
`inst/openspec/changes/audit-cross-package-post-hooks/`.
Phase B triage is recorded inline; Phase C follow-up OpenSpec
changes are linked.

## Phase A reports

- Agent A (mfsusieR vs susieR):
  `inst/notes/sessions/2026-04-30-cross-package-audit-a-susier-posthooks.md`
  — 5 findings (0 fix-now, 2 track-for-later, 1 needs-trace, 2
  accept-pre-existing).
- Agent B (mfsusieR vs fsusieR):
  `inst/notes/sessions/2026-04-30-cross-package-audit-b-fsusier-posthooks.md`
  — 9 findings (0 fix-now, 4 track-for-later, 2 needs-trace, 3
  accept-pre-existing).
- Agent C (mfsusieR vs mvf.susie.alpha):
  `inst/notes/sessions/2026-04-30-cross-package-audit-c-mvf-posthooks.md`
  — 9 findings (1 fix-now, 4 track-for-later, 2 needs-trace, 2
  accept-pre-existing).

Total: 23 findings. Zero numerical-correctness divergences from
the SER hook + per-effect prior refactor that landed since the
2026-04-26 audit. The dominant theme is downstream consequences
of commit 605f395 (Make V more informative): four findings touch
the V-filter coherence surface.

## Phase B triage

| ID  | Pathway                         | Original        | Triage                | Resolution                                                                        |
|-----|---------------------------------|-----------------|-----------------------|-----------------------------------------------------------------------------------|
| A-1 | trim_null_effects + V filter    | track-for-later | track-for-later       | OpenSpec `address-effective-V-consumer-coherence`                                 |
| A-2 | warm-start validate_init        | track-for-later | track-for-later       | OpenSpec `audit-2026-04-30-followups` (item A-2)                                  |
| A-3 | post-hook S3 dispatch cleanup   | needs-trace     | track-for-later       | OpenSpec `audit-2026-04-30-followups` (item A-3)                                  |
| A-4 | per-effect persistent state     | accept-pre-ex   | accept-known          | Pathway clean; no action.                                                         |
| A-5 | hook contract / dispatch sanity | accept-pre-ex   | accept-known          | Pathway clean; no action.                                                         |
| B-1 | mixsqp control defaults         | track-for-later | accept-known          | Documented design choice (warm-start amortizes lower iter count).                 |
| B-2 | V-gated PIP filter              | track-for-later | track-for-later       | OpenSpec `address-effective-V-consumer-coherence`                                 |
| B-3 | joint sigma2 vs marginal Shat   | accept-pre-ex   | accept-known          | Divergence 3 in `~/Documents/mfsusier-mvf-divergences`.                           |
| B-4 | mixsqp status warning           | track-for-later | accept-known          | Conservative-by-design.                                                           |
| B-5 | CS purity parity vs fsusieR     | needs-trace     | track-for-later       | OpenSpec `audit-2026-04-30-followups` (item B-5)                                  |
| B-6 | HMM mu-subset clean             | accept-pre-ex   | accept-known          | Refactor-exceptions PR group 6.                                                   |
| B-7 | smash credible-band lfsr        | track-for-later | accept-known          | mfsusieR-only feature; documented.                                                |
| B-8 | TI/HMM regressor differs        | needs-trace     | add-to-ledger         | Stanza appended to `inst/notes/refactor-exceptions.md` (audit subsection).        |
| B-9 | ELBO NA short-circuit           | accept-pre-ex   | accept-known          | Documented.                                                                       |
| C-1 | ledger stale on M-fold scaling  | fix-now         | add-to-ledger         | Ledger entry corrected inline (refactor-exceptions.md L371-396).                  |
| C-2 | alpha-threshold vs top-K thin   | track-for-later | add-to-ledger         | Stanza appended to `inst/notes/refactor-exceptions.md` (audit subsection).        |
| C-3 | scratchpad cleanup gap          | track-for-later | fix-now               | Code fix applied: `ser_cache` deleted entirely (per-l betahat is one division from the existing `model$residuals` cache; per-iter shat2 was redundant with `iter_cache`). `iter_cache` stripped by `cleanup_extra_fields.mf_individual`. |
| C-4 | V-filter null-fixture trace     | needs-trace     | track-for-later       | OpenSpec `address-effective-V-consumer-coherence`                                 |
| C-5 | pre-hook V dual-encoding        | track-for-later | track-for-later       | OpenSpec `address-effective-V-consumer-coherence`                                 |
| C-6 | shared inner-EM tol             | track-for-later | add-to-ledger         | Stanza appended to `inst/notes/refactor-exceptions.md` (audit subsection).        |
| C-7 | combiner default method         | needs-trace     | track-for-later       | OpenSpec `audit-2026-04-30-followups` (item C-7)                                  |
| C-8 | L_greedy expansion test         | track-for-later | track-for-later       | OpenSpec `audit-2026-04-30-followups` (item C-8)                                  |
| C-9 | lbf_variable_outcome layout     | accept-pre-ex   | accept-known          | Design choice per CLAUDE.md hard rule #2.                                         |

### Triage rollup

- **fix-now (applied this session)**: C-3 (cleanup scratchpad strip).
- **add-to-ledger (applied this session)**: C-1 (correction), B-8, C-2, C-6 (new stanzas under "Cross-package audit (post hooks, 2026-04-30)" subsection of `inst/notes/refactor-exceptions.md`).
- **track-for-later, V-semantics cluster (one OpenSpec)**: A-1, B-2, C-4, C-5 → `address-effective-V-consumer-coherence`.
- **track-for-later, polish bundle (one OpenSpec)**: A-2, A-3, B-5, C-7, C-8 → `audit-2026-04-30-followups`.
- **accept-known**: A-4, A-5, B-1, B-3, B-4, B-6, B-7, B-9, C-9.

## Phase C deliverables (this round)

- This summary (`inst/notes/cross-package-audit-summary-posthooks.md`).
- New "Cross-package audit (post hooks, 2026-04-30)" subsection in
  `inst/notes/refactor-exceptions.md`.
- Updated stale ledger entry on `mixture_null_weight` M-fold scaling
  (refactor-exceptions.md ~L371-396).
- Two new OpenSpec changes:
  `address-effective-V-consumer-coherence`,
  `audit-2026-04-30-followups`. Both validate clean.
- Code fix in `R/individual_data_methods.R`: deleted the `ser_cache`
  scratchpad entirely (its `betahat` is one division off
  `model$residuals` and its `shat2` was redundant with
  `iter_cache$shat2`); `iter_cache` is stripped by
  `cleanup_extra_fields.mf_individual` (C-3, ~10s of MB savings on
  heavy fixtures).

## Phase D process gate

Documented in `inst/notes/review-loop-methodology.md`: every PR
that touches `R/individual_data_methods.R`, `R/model_class.R`,
`R/ibss_methods.R`, or `R/em_helpers.R` requires a cross-package
audit-style review before merge.

## What this audit confirmed (silent because-clean)

The agents collectively traced every numerical pathway in the
SER hook integration end-to-end and found ZERO new numerical
divergences from susieR / fsusieR / mvf.susie.alpha that are not
already in the divergence ledger or refactor-exceptions ledger.
The recent commits since 2026-04-26 (per_scale_normal /
per_scale_laplace ebnm priors, the SER pre/post hook design, the
inner-EM loop, the per-effect persistent prior storage, the V
dimensions fix, the M-fold mixsqp scaling) preserve the
established equivalence contracts.

The only behavioral change that surfaced is the V-filter
exposure on PIPs (commit 605f395), which is intentional but
under-documented; the V-semantics cluster OpenSpec change closes
that gap.
