# mf-small-sample-correction capability

## ADDED Requirements

### Requirement: `small_sample_correction` is exposed on `mfsusie()` and `fsusie()`

`mfsusie()` and `fsusie()` SHALL accept `small_sample_correction` as a public logical argument with default `FALSE`. When `TRUE`, the SER step SHALL replace the per-variable Wakefield Normal marginal Bayes factor with a Johnson 2005 scaled Student-t marginal Bayes factor with `df = data$n - 1` degrees of freedom.

#### Scenario: signature

```
mfsusie(X, Y, pos = NULL, L = ...,
        small_sample_correction = FALSE,
        ...)
```

`fsusie()` forwards the argument via `...`.

### Requirement: default path is bit-identical to a call without the argument

`mfsusie(..., small_sample_correction = FALSE)` SHALL produce the same fit as a call that omits the argument. Both calls use the Wakefield Normal kernel.

The numeric outputs `alpha` and `pip` SHALL agree at tolerance `0`.

#### Scenario: default no-op

For a small fixed seed dataset, the explicit-`FALSE` and implicit-default fits are equal at tolerance `0`.

### Requirement: Johnson-t kernel matches the upstream marginal at machine precision

The mfsusieR Johnson-t kernel `mixture_log_bf_per_scale_johnson` summed across scales SHALL match `fsusieR::log_BF` with `df = n - 1` at machine precision. This is the apple-to-apple parity contract for the Johnson-t branch (the upstream kernel is the same scaled-t mixture density). The two implementations sum the same terms in slightly different order; floating-point non-associativity produces sub-ULP-relative differences that vanish as the summation order is fixed.

#### Scenario: per-variable LBF parity

For a fixed seed dataset with `n = 50`, `J = 15`, `T_m = 64`, mixture grid `sd = c(0, 0.5, 1.5)`, weights `pi = c(0.6, 0.2, 0.2)`, and `df = n - 1`, the per-variable LBFs from the two implementations agree at `tolerance <= 1e-14` (max absolute deviation ~3.55e-15, max relative deviation ~2.08e-16, i.e. one ULP).

### Requirement: Johnson-t reduces null-variant PIP relative to Wakefield on small `n`

`mfsusie(..., small_sample_correction = TRUE)` SHALL not increase the aggregate PIP at known-null variants relative to the Wakefield default on small-sample data. The signal variant SHALL still be recovered (`PIP > 0.9`).

#### Scenario: structural correction at small `n`

A simulated dataset with `n = 40`, `T_m = 32`, `p = 30`, and one signal variant is fit twice. The signal-variant PIP under Johnson-t exceeds `0.9`. The aggregate PIP across non-signal variants under Johnson-t is no larger than under the Wakefield default (within numerical tolerance `1e-6`).

### Requirement: argument validation rejects non-logical input

The argument validator SHALL error when `small_sample_correction` is not a length-1 logical (`TRUE` or `FALSE`).

#### Scenario: invalid input

A call with `small_sample_correction = "yes"` errors with a message identifying the required type.
