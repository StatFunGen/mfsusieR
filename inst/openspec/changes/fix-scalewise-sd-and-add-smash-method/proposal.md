# Fix `scalewise` pointwise SD and add `smash` smoother

## Why

Two issues with the post-fit smoother surface:

1. `mf_post_smooth(method = "scalewise")` uses
   `inverse-DWT(sqrt(var_w))` for the position-space pointwise
   standard deviation. The inverse DWT is linear, so this
   formula is the absolute value of a linear combination of
   per-scale standard deviations, NOT the position-space
   posterior SD. The correct formula for an orthonormal
   wavelet basis `W` is
   `sd_pos[t] = sqrt( sum_k W[t, k]^2 * var_w[k] )`. The
   current formula over- or under-states the per-position SD
   wherever the wavelet basis mixes coefficients across
   positions. The vignette
   `vignettes/post_processing.Rmd` justifies the wrong
   formula by appealing to Parseval's theorem, which applies
   to total energy, not per-position decomposition.
2. `fsusieR::univariate_smash_regression` (selected by
   `post_processing = "smash"`) is one of the three
   established post-smoothers in the upstream functional
   fine-mapping pathway. mfsusieR currently only exposes
   `"scalewise"`, `"TI"`, and `"HMM"`. Adding `"smash"`
   closes this gap.

This change also retires the precondition check in
`mf_post_smooth` that requires `fit$residuals` and
`fit$lead_X` for every method; `scalewise` does not read
either field.

## What changes

### 1. Fix scalewise pointwise SD

Replace the per-position SD formula in `.post_smooth_scalewise`:

```
old: sd_pos <- abs(mf_invert_dwt(matrix(sqrt(var_w), 1L), ...))
new: sd_pos <- sqrt(.invert_variance_curve(var_w, ...))
```

`.invert_variance_curve` is a new internal helper that
applies the squared-filter-coefficient inverse-DWT pattern,
mirroring the `wd_variance` / `av_basis_variance` helpers
already in `R/utils_wavelet.R`. The result is the correct
position-space variance vector.

Update `vignettes/post_processing.Rmd`: drop the
"Parseval"-based justification and replace with the linear-
combination derivation.

### 2. Add `method = "smash"`

`mf_post_smooth(method = "smash")`. Requires the smashr
package; gated on `requireNamespace("smashr",
quietly = TRUE)`. Body:

- For each effect `l` and outcome `m`, compute the per-
  position regression estimate of the lead variable against
  the in-sample residual response.
- Pass to `smashr::smash.gaus(est, sigma = sd_grid)`.
- Pointwise band reuses the same sigma grid.

The implementation matches `fsusieR::univariate_smash_regression`.
Bit-identity test at tolerance `<= 1e-12`.

### 3. Default method change (pre-1.0)

`mf_post_smooth(method = c("TI", "scalewise", "HMM",
"smash"))` defaults to `"TI"` (was `"scalewise"`). TI is the
default in the upstream functional fine-mapping case studies
because it carries pointwise credible bands without extra
dependencies beyond what mfsusieR already imports.

The pre-1.0 default change is mentioned in `NEWS.md`.

### 4. Method-aware precondition check

`mf_post_smooth` precondition check is method-aware. The
new contract:

- `"scalewise"`: requires only `fit$mu`, `fit$mu2`, the
  scale index, and the wavelet metadata.
- `"TI"`, `"HMM"`, `"smash"`: require `fit$residuals` and
  `fit$lead_X`.

A user calling
`mf_post_smooth(scalewise_only_fit, method = "scalewise")`
on a fit that lacks the residuals slot succeeds. Calling it
with `method = "TI"` errors with a clear message.

### 5. Tests

- `tests/testthat/test_post_smooth_scalewise.R`: tighten
  the existing structural test with a closed-form check
  on the corrected pointwise SD against an exact reference
  computation on a small fixture.
- `tests/testthat/test_post_smooth_smash.R` (new): bit-
  identity vs `fsusieR::univariate_smash_regression` at
  tolerance `<= 1e-12`. Gated on
  `requireNamespace("smashr")`.

### 6. Docs

- Roxygen for `mf_post_smooth` enumerates all four methods
  and notes the dependency requirements.
- `vignettes/post_processing.Rmd` rewritten to present
  TI as the recommended default (cred bands available), the
  three alternatives (scalewise: fast no-extra-deps;
  HMM: lfsr; smash: wavelet-shrinkage reference), and the
  corrected scalewise SD derivation.

## Impact

- Changed: `R/mfsusie_methods.R` (default + method-aware
  precondition + scalewise SD); `R/utils_wavelet.R` (new
  variance-inverse helper); `vignettes/post_processing.Rmd`
  (rewrite + corrected derivation).
- New: `R/post_smooth_smash.R`,
  `tests/testthat/test_post_smooth_smash.R`.
- Pre-1.0 default change for `mf_post_smooth`. Documented
  in `NEWS.md`.
- DESCRIPTION: smashr stays in Suggests; the smash method
  gates on `requireNamespace`.
- pixi.toml: r-smashr stays under
  `[target.linux-64.dependencies]`. osx-arm64 vignette
  builds skip the smash chunk.

## Out of scope

- Reimplementing `smash.gaus`. We delegate to smashr.
