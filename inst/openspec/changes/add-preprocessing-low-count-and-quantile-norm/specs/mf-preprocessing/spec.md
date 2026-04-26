# mf-preprocessing capability

## ADDED Requirements

### Requirement: `low_count_filter` masks sparse wavelet columns

`mfsusie()`, `fsusie()`, and `mf_adjust_for_covariates()` SHALL
accept a `low_count_filter` argument with default `0`. The
argument SHALL be a non-negative scalar threshold; columns of
the wavelet-domain matrix `Y_wd = cbind(W$D, W$C)` with
`median(|column|) <= low_count_filter` SHALL be flagged. The
flagged column index set SHALL be threaded through both the
prior initialization step and the per-iteration marginal
regression step. In the marginal regression step the masked
columns SHALL have `Bhat` set to `0` and `Shat` set to `1`
after the call to `susieR::compute_marginal_bhat_shat()`.

#### Scenario: zero-median column case at default threshold

When called with `low_count_filter = 0` on data containing at
least one wavelet column whose absolute-value median is exactly
`0`, `mfsusie()` SHALL produce the same fit as the upstream
functional fine-mapping pathway at the same default. Tolerance
`<= 1e-12` on `alpha`, `mu`, `mu2`, `lbf`, `pip`, and exact
match on `sets$cs`.

#### Scenario: positive threshold sweep

For `low_count_filter in {0.1, 0.5}` the fit SHALL match the
upstream result at tolerance `<= 1e-12` on the same numeric
fields.

### Requirement: `quantile_norm` applies column-wise rank-INT to wavelet coefficients

`mfsusie()`, `fsusie()`, and `mf_adjust_for_covariates()` SHALL
accept a `quantile_norm` argument with default `FALSE`. When
`TRUE`, the column-wise rank-based normal quantile transform
`qqnorm(rank(., ties.method = "random"))$x` SHALL be applied to
`Y_wd = cbind(W$D, W$C)` after the low-count flagging and
before the IBSS loop. The transform SHALL use `set.seed(1)` for
tie reproducibility, matching the upstream convention.

#### Scenario: bit-identity at `quantile_norm = TRUE`

`mfsusie(..., quantile_norm = TRUE)` SHALL match the upstream
functional fine-mapping pathway with the analogous flag set,
at tolerance `<= 1e-12` on every numeric output.

### Requirement: marginal regression delegates to susieR

The IBSS loop SHALL NOT re-implement marginal regression in
mfsusieR. The wavelet-domain `Bhat` / `Shat` SHALL be obtained
by calling `susieR::compute_marginal_bhat_shat()` and applying
the low-count mask post-compute via a thin wrapper. The wrapper
SHALL be the only mfsusieR function authorized to write into
the masked positions.

#### Scenario: wrapper signature

`mf_compute_bhat_shat_with_low_count_mask(X, Y_wd, lowc_idx)`
SHALL return a list with `Bhat` and `Shat` matrices identical
to `susieR::compute_marginal_bhat_shat(X, Y_wd)` except at
columns indexed by `lowc_idx`, where `Bhat` is `0` and `Shat`
is `1`.

### Requirement: `mf_adjust_for_covariates` honors both preprocessing options

`mf_adjust_for_covariates()` SHALL accept and act on `low_count_filter` and `quantile_norm` arguments with the same semantics as `mfsusie()`.

When `quantile_norm = TRUE` the function SHALL return `Y_adjusted` on the original (un-INT'd) position scale so that `fsusie(Y_adjusted, X)` operates on the same units as the input `Y`.

#### Scenario: round-trip on `Y_adjusted` units

For `quantile_norm = TRUE`, the absolute deviation between
`Y` and `Y_adjusted` shall be small (controlled by the
adjustment magnitude) and on the same scale as `Y`, with no
column of `Y_adjusted` having a standardized variance close
to `1` purely as a consequence of the transform.

#### Scenario: previously rejected flags now succeed

`mf_adjust_for_covariates(Y, Z, low_count_filter = 0.1)` and
`mf_adjust_for_covariates(Y, Z, quantile_norm = TRUE)` SHALL
return successfully (no v1 reject), and their output SHALL
match the corresponding upstream EBmvFR-style result at
tolerance `<= 1e-12` on `fitted_func`.
