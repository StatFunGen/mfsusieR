# mf-residual-variance-method capability

## ADDED Requirements

### Requirement: `estimate_residual_method` is exposed on `mfsusie()` and `fsusie()`

`mfsusie()` and `fsusie()` SHALL accept
`estimate_residual_method` as a public argument with values
`c("MoM", "MLE", "NIG")` and default `"MoM"` to match
`susieR::susie()`. The argument SHALL be forwarded to the
IBSS workhorse via `params$use_NIG`.

#### Scenario: signature

```
mfsusie(X, Y, pos = NULL, L = ...,
        estimate_residual_method = c("MoM", "MLE", "NIG"),
        ...)
```

`fsusie()` forwards the same argument.

### Requirement: NIG matches `susieR::susie` bit-identically on the scalar case

`mfsusie(..., estimate_residual_method = "NIG")` SHALL match `susieR::susie(..., estimate_residual_method = "NIG")` bit-identically on the scalar (M=1, T_1=1) case.

The numeric outputs `sigma2`, `alpha`, `mu`, `mu2`, `lbf`, and `pip` SHALL agree at tolerance `<= 1e-12`.

#### Scenario: scalar bit-identity sweep

For `n in {50, 80}`, `p in {200, 500}`, fixed seed, the bit-
identity holds.

### Requirement: NIG corrects small-sample inflation on multi-outcome data

`mfsusie(..., estimate_residual_method = "NIG")` SHALL produce a tighter `sigma^2` and lower null-variant PIPs than the MoM default on small-sample multi-outcome data.

For multi-outcome data with `n = 80` and a known set of null variants, the per-outcome `sigma2` under NIG SHALL be `>=` the corresponding MoM `sigma2`, and the maximum PIP at known-null variants under NIG SHALL be `<=` the maximum PIP at the same variants under MoM.

#### Scenario: small-sample correction sanity

A simulated fixture with `n = 80`, `p = 100`, `M = 3` is fit
twice. The asserted directions of the differences hold; this
is a structural, not numerical, test.
