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
