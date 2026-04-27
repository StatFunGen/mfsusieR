# mfsusieR writing style policy

This file binds the writing voice for every artifact that ships
in the package: `R/` roxygen, `data-raw/` description fields,
`vignettes/*.Rmd`, plot titles, figure captions, NEWS entries,
and any prose visible to a downstream user.

The voice is the Wang group `statgen-writing-style` (loaded as
the `anthropic-skills:statgen-writing-style` skill). It is
**mechanical, direct, quantitative**. The reader is a
methods-literate collaborator who reads R help pages and
manuscripts, not an outsider who needs motivation or teaching.

Internal session notes (`inst/notes/sessions/*`), OpenSpec
planning artifacts (`inst/openspec/changes/*/proposal.md`,
`design.md`, `tasks.md`), `inst/CLAUDE.md`, and developer test
files (`tests/testthat/`) MAY use a more relaxed register
because they are not user-facing. Spec files
(`inst/openspec/specs/*/spec.md`) MUST follow this policy
because they form the public contract.

## Hard ban list

These words and phrases SHALL NOT appear in any user-facing
artifact, regardless of context:

- "didactic", "didactic example"
- "fixture" or "fixtures" when referring to a shipped
  example dataset (testthat fixtures inside `tests/` are fine)
- "teaching aid", "teaching tool", "teaching example"
- "designed to surface ...", "designed to demonstrate ..."
- "shaped after", "fsusie()-shaped", "mfsusie()-shaped"
- "tuned to ..."
- "dramatic", "dramatic contrast"
- "clean and reproducible", "clean recovery"
- "visible at a glance", "easy to see"
- "story" used as a synonym for "section" or "narrative"
- "showcase", "showcases"
- "out of the box"
- "knockoff", "cheap knockoff"
- Tier-1 phrases from the statgen-writing-style skill:
  "striking", "strikingly", "pervasive", "widespread"
  (as motivational framing), "intriguing", "fascinating",
  "paradigm shift", "cutting-edge", "compelling"
- Tier-2 phrases require earned emphasis: "notably",
  "importantly", "remarkably", "novel"

### Internal-jargon ban

The following internal references SHALL NOT appear in user-facing
artifacts (vignettes, `R/` roxygen, `R/` inline comments, NEWS,
README, plot titles, spec files). They are valid in
`inst/notes/`, `inst/openspec/`, `tests/testthat/` headers,
and session notes only.

- "upstream" (and "upstream-equivalent", "matches upstream",
  "upstream defaults", etc.) when referring to port sources.
  Describe the behaviour directly; do not contrast against an
  unnamed external reference.
- "port source", "port-source", "port source's"
- Names of port-source packages and their internal symbols:
  `mvf.susie.alpha`, `multfsusie`, `EBmvFR`, `cor_small`,
  `fsusieR::univariate_TI_regression`,
  `fsusieR::univariate_HMM_regression`,
  `fsusieR::univariate_smash_regression`,
  `fsusieR::log_BF`, `fsusieR::cal_Bhat_Shat`, etc. Mentioning
  `susieR::susie()` or `susieR::compute_marginal_bhat_shat`
  is fine (susieR is a public runtime dependency).
- Contract IDs `C1`, `C2`, `C3` (the apple-to-apple equivalence
  contracts).
- "PR group", "PR group N" labels and any "(PR group ...)"
  parentheticals.
- "legacy mode", "legacy variance mode", "legacy output" when
  referring to the original behaviour pre-port. State the
  current mode by its public name.
- "fidelity contract", "degeneracy contract", "binary tolerance
  philosophy" — these are doctrine terms.
- "Pattern A", "Pattern B", "kept-upstream-bug",
  "port-source-bug fix" — refactor-discipline doctrine terms.

When user-facing prose needs to explain why a particular
formula or default was chosen, describe the mathematics or the
empirical behaviour directly. The reasoning lives in the
manuscript or in `inst/notes/` for internal readers.

## Acceptable replacements

| Banned | Use instead |
|---|---|
| fixture | "example dataset" / "dataset" / describe what it is |
| didactic example | "small simulated example", or describe the size and what is in it |
| shaped after Y | "An `fsusie()` example with ..." or name the model directly |
| designed to surface X | describe what is in the data; let the reader infer the use |
| visible at a glance | drop the phrase; the reader will look at the plot themselves |
| story (the prior-comparison story) | "section" / "panel" / "subsection" |
| teaching aid | drop; describe the mechanical role |

## Voice and structure

- Open each paragraph with the finding or claim, not framing.
- Past tense for completed analyses (*we ran*, *fsusie identified
  three credible sets*); present tense for method descriptions
  (*fsusie returns the per-effect mean and credible band*).
- First person plural for collaborative prose; passive only when
  the agent is irrelevant.
- Quantitative claims state numbers (`n = 100`, `T = 32`,
  `PIP = 0.97`), not adjectives (*small*, *high*).
- No em dashes. Use commas, semicolons, parentheses.
- No colons in headings. Use a descriptive phrase.

## Audit before commit

Before staging any change to user-facing prose, search the
modified files for the hard ban list above. If a banned phrase
appears, rewrite. If unsure, invoke the `statgen-writing-style`
skill and feed it the paragraph for review.

The audit applies to:

- `R/data.R`, `R/*.R` roxygen blocks
- `data-raw/*.R` `description` fields and comment headers
- All `vignettes/*.Rmd`
- `inst/openspec/specs/*/spec.md` (the public contract)
- `NEWS.md`, `README.md`
- Plot titles, axis labels, legend entries built into `R/`

Internal artifacts (sessions, planning proposals, CLAUDE.md)
are exempt.
