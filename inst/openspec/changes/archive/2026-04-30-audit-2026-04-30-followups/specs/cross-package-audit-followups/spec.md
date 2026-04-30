# Cross-package audit follow-ups spec delta

## ADDED Requirements

### Requirement: warm-start fit validator

mfsusieR SHALL validate `params$model_init` shape and content
before the IBSS loop runs. The validator SHALL check that
`alpha`, `mu`, `mu2`, `pi_V`, `fitted_g_per_effect`, `sigma2`,
`V` contain no NA or Inf and conform to the list-of-list shapes
expected by the `mf_individual` IBSS path.

#### Scenario: warm-start fit with NA in alpha

- **WHEN** a user calls `mfsusie(..., model_init = bad_fit)`
  with `bad_fit$alpha` containing an NA value
- **THEN** the validator SHALL error with a message naming the
  offending field, BEFORE the IBSS loop runs

### Requirement: extension seam for cross-outcome combiners

The `combine_outcome_lbfs` generic SHALL provide a `.default`
method that errors with a message naming the registered
combiner subclasses, so future combiner subclasses surface a
clear "method not implemented" failure rather than UseMethod's
opaque dispatch error.

#### Scenario: user supplies a combiner with no registered method

- **WHEN** an internal call dispatches `combine_outcome_lbfs` on
  an object whose class chain has no matching method
- **THEN** the `.default` method SHALL error with a message
  listing the methods available on the generic
