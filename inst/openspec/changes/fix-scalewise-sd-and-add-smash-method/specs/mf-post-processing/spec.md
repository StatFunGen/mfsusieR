# mf-post-processing capability (deltas)

## MODIFIED Requirements

### Requirement: `mf_post_smooth` supports four methods with `"TI"` as default

`mf_post_smooth(fit, method = ...)` SHALL accept
`method in c("TI", "scalewise", "HMM", "smash")` with
default `"TI"`. The `"smash"` value SHALL gate on
`requireNamespace("smashr", quietly = TRUE)` and error
clearly when smashr is not installed.

#### Scenario: default routes to TI

`mf_post_smooth(fit)` SHALL produce the same result as
`mf_post_smooth(fit, method = "TI")`.

#### Scenario: smash availability check

When `smashr` is not installed and the user calls
`mf_post_smooth(fit, method = "smash")`, the function SHALL
error with a message that names smashr and points to the
install instructions.

### Requirement: scalewise pointwise SD uses the linear-combination variance formula

The scalewise pointwise standard deviation SHALL use the linear-combination variance formula `sd_pos[t] = sqrt( sum_k W[t, k]^2 * var_w[k] )` where `W` is the inverse-DWT operator on the dyadic grid.

The previous formula `inverse-DWT(sqrt(var_w))` SHALL NOT be used; the "Parseval" justification SHALL be removed from the post-processing vignette.

#### Scenario: closed-form match on a small example

For a dyadic example dataset of length 8 with a fixed input variance vector `var_w`, the returned `sd_pos` SHALL match the exact `sqrt(W^2 %*% var_w)` reference computation at tolerance `<= 1e-14` (machine precision).

### Requirement: smash method matches upstream bit-identically

`mf_post_smooth(fit, method = "smash")` SHALL match
`fsusieR::univariate_smash_regression` output at tolerance
`<= 1e-12` on `effect_curve` and `cred_band` for the M=1
case across the standard fixture sweep.

#### Scenario: bit-identity sweep

For `n = 60`, `T = 64`, fixed seed, the per-position
`effect_estimate` and `cred_band` returned by the in-package
`univariate_smash_regression` SHALL match the upstream
functional fine-mapping post-processor at tolerance
`<= 1e-12`.

### Requirement: precondition check is method-aware

`mf_post_smooth` SHALL require `fit$residuals` and
`fit$lead_X` only for `method in c("TI", "HMM", "smash")`.
For `method = "scalewise"` neither field is required.

#### Scenario: scalewise on a fit without residuals

`mf_post_smooth(fit_lacking_residuals, method = "scalewise")`
SHALL succeed. Calling the same fit with `method = "TI"`
SHALL error with a message naming the missing slot.
