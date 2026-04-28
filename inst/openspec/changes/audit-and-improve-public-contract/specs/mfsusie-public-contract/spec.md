# mfsusieR public contract

## ADDED Requirements

### Requirement: fit object SHALL expose a documented and trimmed field set

`mfsusie()` SHALL return a list whose every field has a
one-line description in the `@return` roxygen of `mfsusie()`.
Fields that duplicate information available from a base
`susie()` field, or that are derivable from other fields
without measurable cost, SHALL be removed. Fields that are
mfsusieR-specific extensions SHALL be flagged as such in the
roxygen.

#### Scenario: documented field set

- **WHEN** a user runs `mfsusie(X, Y, ...)` and reads
  `?mfsusie`
- **THEN** every field in the returned object SHALL appear in
  the `@return` documentation with a one-line description
- **AND** any field not present in `susieR::susie()`'s output
  SHALL be flagged as an mfsusieR extension.

### Requirement: trivial S3 overrides SHALL be removed

Each S3 method registered for `mf_individual` or `mfsusie` SHALL be retained only when its body materially diverges from
`susieR`'s default for that generic. Methods whose body
reproduces or trivially wraps the default SHALL be removed and
the generic SHALL fall through to `susieR`'s default.

#### Scenario: a removed override defers to susieR's default

- **WHEN** mfsusieR removes its trivial override of a generic
  `g`
- **THEN** `g(data, params, model, ...)` for an
  `mf_individual` data object SHALL dispatch to
  `g.default()` from susieR
- **AND** the resulting `mfsusie()` fit on the standard test
  fixtures SHALL be byte-identical to the pre-removal fit.

### Requirement: `verbose = TRUE` SHALL print per-iteration diagnostics

When `mfsusie(verbose = TRUE)` is called, the IBSS loop SHALL
print one line per iteration containing at minimum the
iteration number, ELBO, and maximum per-effect alpha change.

#### Scenario: per-iter verbose output

- **WHEN** `mfsusie(X, Y, verbose = TRUE)` is called
- **THEN** for each IBSS iteration `k`, exactly one line SHALL
  be printed of the form
  `iter k | elbo=<value> | max_alpha_diff=<value> [...]`
- **AND** the iteration number SHALL be printed before the
  iteration body executes.

### Requirement: `convergence_metric` SHALL select between PIP-diff and ELBO

`mfsusie()` SHALL accept a `convergence_metric` argument with
values `"pip_diff"` (default) and `"elbo"`. The IBSS loop SHALL
stop when the chosen metric falls below `tol`.

#### Scenario: ELBO-based convergence on a standard fixture

- **WHEN** the user runs `mfsusie(X, Y, convergence_metric =
  "elbo", tol = 1e-3)`
- **THEN** the fit SHALL converge in a number of iterations
  that differs from the PIP-diff-based fit by at most 2

### Requirement: HMM smoother SHALL return a credible band

`mf_post_smooth(method = "HMM")` SHALL populate
`fit$smoothed$HMM$credible_bands[[m]][[l]]` with a per-
position lower/upper bound at the requested coverage level.
The band formula SHALL be derived from the HMM posterior
mixture and SHALL be documented in the function roxygen.

#### Scenario: HMM band populated

- **WHEN** the user runs `mf_post_smooth(fit, method = "HMM")`
  on a fit with at least one credible set
- **THEN** for each `(m, l)` with a non-null effect curve,
  `fit$smoothed$HMM$credible_bands[[m]][[l]]` SHALL be a
  `T_basis[m] x 2` matrix
- **AND** the lower-bound column SHALL not exceed the
  upper-bound column at any position.
