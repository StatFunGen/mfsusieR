# Design — `mf_adjust_for_covariates()`

## D1. Algorithm map

The wavelet-EB path reproduces `fsusieR::EBmvFR(Y, X = Z, adjust = TRUE)`
at upstream defaults by assembling existing primitives:

```
INPUT
  Y_pos (n x T_pos)     functional response
  Z     (n x K)         covariates to adjust away
  X     (n x p, opt)    genotype matrix to FWL-residualize
                        (extension; not in upstream)

PIPELINE
  1.  Y_centered = col_scale(Y_pos, scale = FALSE)
  2.  Z_centered_scaled = col_scale(Z, scale = TRUE)
  3.  W = mf_dwt_matrix(Y_centered, filter, family)
      Y_wavelet = cbind(W$D, W$C)            # n x T_padded
  4.  prior  = init_scale_mixture_prior_default(
                 Y_wavelet, Z_centered_scaled,
                 prior_class = "mixture_normal_per_scale",
                 indx_lst, ...)
  5.  obj.fitted_wc  = matrix(0, K, T_padded)
      obj.fitted_wc2 = matrix(0, K, T_padded)
      obj.MLE_wc     = matrix(0, K, T_padded)
      obj.MLE_wc2    = matrix(0, K, T_padded)
      obj.sigma2     = mean(apply(Y_centered, 2, var))
      obj.G_prior    = prior$G_prior
      obj.pi_hist    = list(initial pi from G_prior)

  6.  iter = 1
      while (iter < max_iter && check > tol):
        # 6a. per-covariate Gauss-Seidel
        for j in 1..K:
          partial_resid = Y_wavelet - Z[, -j] %*% obj.fitted_wc[-j, ]
          (Bhat, Shat) = compute_marginal_bhat_shat(
                            X = matrix(Z[, j], ncol = 1),
                            Y = partial_resid)
          Shat = pmax(Shat, 1e-32)
          obj.MLE_wc[j, ]  = Bhat
          obj.MLE_wc2[j, ] = Shat^2
          obj.fitted_wc[j, ]  = ashr::postmean(obj.G_prior, Bhat, Shat,
                                               indx_lst)
          obj.fitted_wc2[j, ] = ashr::postsd  (obj.G_prior, Bhat, Shat,
                                               indx_lst)^2
        # 6b. residual variance (KEPT-UPSTREAM-BUG formula;
        #     see refactor-exceptions Pattern B entry)
        ER2  = sum((Y_wavelet - Z %*% obj.fitted_wc)^2)
             - sum(obj.fitted_wc^2)
             + sum(obj.fitted_wc2)
        obj.sigma2 = ER2 / (n * T_padded)
        # 6c. update prior weights
        L = .compute_per_scale_likelihood_matrix(
              obj.MLE_wc, sqrt(obj.MLE_wc2),
              obj.G_prior, indx_lst, is_ebmvfr = TRUE)
        new_pi = .m_step_per_scale(L, ..., is_ebmvfr = TRUE)
        obj.G_prior = update_prior(obj.G_prior, new_pi)
        obj.pi_hist[[length + 1]] = new_pi
        # 6d. convergence
        check = sum(abs(new_pi - prev_pi)) / log(length(new_pi))
        iter += 1

  7.  fitted_func = matrix(0, K, T_pos)
      for j in 1..K:
        fitted_func[j, ] = mf_invert_dwt(obj.fitted_wc[j, ],
                                         filter, family,
                                         T_pos)
        # rescale by csd_X[j] (because Z was scaled)
        fitted_func[j, ] = fitted_func[j, ] / csd_X[j]
  8.  Y_adjusted = Y_pos - Z %*% fitted_func
      X_adjusted = (if X supplied) X - Z %*% (Z'Z)^{-1} Z' X
  9.  return list(Y_adjusted, X_adjusted, fitted_func, sigma2,
                  niter, converged, method = "wavelet_eb")
```

## D2. Why the buggy `get_ER2.EBmvFR` formula is kept

**See `inst/notes/refactor-exceptions.md` for the formal entry.**

The fsusieR formula `sum(resid^2) - sum(postF^2) + sum(postF2)`
omits `predictor_weights = colSums(Z^2)`. The mathematically
correct formula (per the susieR-style decomposition we landed
in `mf_get_ER2_per_position`) would include the `pw` factor.

Three reasons we keep the upstream formula:

1. **Apple-to-apple contract.** The unit-test contract for
   `mf_adjust_for_covariates(method = "wavelet_eb")` is bit-
   identity to `fsusieR::EBmvFR(adjust = TRUE)` at
   `tolerance <= 1e-8`. Using the corrected formula would
   produce a different `obj$sigma2` and break the contract.
2. **The bug is benign for `Y_adjusted`.** `sigma2` is dead-
   end in EBmvFR's outer loop: `cal_Bhat_Shat` ignores its
   `resid_var` argument and computes `Shat` directly from
   `colVars(Y - X*Bhat) / (n - 1)`. The buggy `sigma2` only
   affects the stored field in the returned object, not any
   subsequent posterior. `Y_adjusted` is identical under
   either formula.
3. **Scope discipline.** Per CLAUDE.md hard rule #1 the full
   `EBmvFR` model is out of scope. We are not re-deriving
   EBmvFR's variance estimator; we are reproducing what
   upstream does for the `adjust = TRUE` use case.

The corrected `mf_get_ER2_per_position` (used in our SuSiE
path) is NOT touched. The two paths use different ER2
functions by design.

## D3. The two `is_ebmvfr` flag plumbing points

The opus-model verification flagged that fsusieR's
`L_mixsq` and `m_step` both take an `is.EBmvFR` flag that
changes:

- `L_mixsq` with `is.EBmvFR = TRUE` does NOT prepend a null-
  component penalty row to the likelihood matrix.
- `m_step` with `is.EBmvFR = TRUE` does NOT prepend the
  `null_prior_weight * length(...)` weight to the EM weight
  vector.

If we reused our existing SuSiE-flavored M-step verbatim, we
would get a different L matrix and different weights and the
final mixture-weights vector would diverge.

**Refactor:** factor out the two functions in
`R/em_helpers.R` with an `is_ebmvfr` flag. The SuSiE call site
passes `is_ebmvfr = FALSE` (current behavior). The new
EBmvFR-adjust call site passes `is_ebmvfr = TRUE`. Both call
sites exercise the same code, which catches refactor breakage
on either side.

## D4. The two unsupported preprocessing branches

Upstream `EBmvFR` has two preprocessing options that are NOT
supported in v1 of `mf_adjust_for_covariates`:

| Upstream flag | Default | v1 behavior | Tracked |
|---|---|---|---|
| `thresh_lowcount` | `0` | error: "not supported in v1" | follow-up |
| `quantile_trans`  | `FALSE` | error: "not supported in v1" | follow-up |

Both are off by default upstream, so the bit-identity contract
holds at defaults. A user who passes either flag gets a clear
error rather than silent divergence.

## D5. Why the OLS branch exists at all

The closed-form `Y_adjusted = (I - Z(Z'Z)^{-1}Z') Y` is
mathematically the no-shrinkage limit of the wavelet-EB path.
It is faster (no IBSS loop, no DWT, no ash) and is appropriate
when the covariates are known to be small in number, well-
conditioned, and the user does not want EB shrinkage.

It is exposed because:

- Some users find the EB-shrunk output suspicious and want a
  closed-form sanity check.
- The OLS path needs no bug-preservation contract: it is
  defined by its math.
- The unit test for the OLS path is closed-form sanity
  (orthogonality, idempotency), not bit-match to upstream.

## D6. The `X` (genotype) FWL extension

Upstream `EBmvFR(adjust = TRUE)` returns only `Y_adjusted`. It
does NOT residualize the genotype matrix. This is a legitimate
choice when Z and X are orthogonal (e.g., Z = principal
components computed on a different cohort), but is incorrect
when Z and X are correlated.

We expose `X` as an optional argument; when supplied, we apply
the FWL correction:

```
H_Z = Z (Z'Z)^{-1} Z'
X_adjusted = (I - H_Z) X
```

This is the SAME projection used in our existing OLS hat-
matrix code in `fsusie_covariates_and_coloc.Rmd`. We retain it
because it is the statistically correct treatment when Z and
X share variance.

For bit-identity to upstream the user calls
`mf_adjust_for_covariates(Y, Z, X = NULL)` (no X argument);
the result on `Y_adjusted` matches `fsusieR::EBmvFR(adjust =
TRUE)$Y_adjusted` exactly.

## D7. Vignette structure

`vignettes/fsusie_covariates_adjustment.Rmd` follows the
shape of `fsusieR::Adjusting_covariate.Rmd`:

```
1. Overview
2. Simulated example
   - Build Y with confounding by 3 random covariate functions
3. Naive fit (no adjustment)
   - mfsusie_plot of mis-attributed signals
4. Adjusted fit
   - mf_adjust_for_covariates(Y, Z) -> Y_adjusted
   - fsusie(Y_adjusted, X) -> recovered fit
5. Compare
   - PIPs of adjusted vs naive
6. The OLS option
   - Brief note on method = "ols" for the no-shrinkage case
```

No port-source mentions per the standing rule.

## D8. The `fsusie_colocalization.Rmd` rename

Mechanical: take `vignettes/fsusie_covariates_and_coloc.Rmd`,
delete the covariates section (lines ~36-69 inclusive), keep
the colocalization section, rename the file. Update the
intro paragraph to drop "(i) regressing out covariates" from
the scope-statement so the vignette description is honest.

## D9. Test plan

`tests/testthat/test_adjust_covariates.R`:

```
test_that("wavelet_eb matches fsusieR::EBmvFR(adjust=TRUE) bit-identically", {
  skip_if_not_installed("fsusieR")
  for (n in c(100L, 200L)) {
    for (T_m in c(64L, 128L)) {
      for (K in c(2L, 3L, 5L)) {
        sim <- build_sim(seed = ..., n = n, T_m = T_m, K = K)
        ours <- mf_adjust_for_covariates(sim$Y, sim$Z,
                                         method = "wavelet_eb")
        ref  <- fsusieR:::EBmvFR(sim$Y, X = sim$Z, adjust = TRUE,
                                  verbose = FALSE)
        expect_equal(ours$Y_adjusted, ref$Y_adjusted,
                     tolerance = 1e-12)
        expect_equal(ours$fitted_func, ref$fitted_func,
                     tolerance = 1e-12)
        expect_equal(ours$sigma2, ref$sigma2,
                     tolerance = 1e-12)
        expect_equal(ours$niter, ref$niter)
      }
    }
  }
})

test_that("ols path is closed-form sanity-correct", {
  ...orthogonality and idempotency assertions...
})

test_that("error paths fire cleanly", {
  expect_error(mf_adjust_for_covariates(Y, Z,
                                        thresh_lowcount = 0.1),
               "not supported in v1")
  expect_error(mf_adjust_for_covariates(Y, Z,
                                        quantile_trans = TRUE),
               "not supported in v1")
})
```

`tests/testthat/test_em_helpers.R`:

```
test_that(".compute_per_scale_likelihood_matrix matches both flag values", {
  ...is_ebmvfr = TRUE: no penalty row...
  ...is_ebmvfr = FALSE: penalty row prepended...
})

test_that(".m_step_per_scale matches both flag values", {
  ...is_ebmvfr = TRUE: no nullweight prepend...
  ...is_ebmvfr = FALSE: nullweight prepend...
})
```

If the SuSiE-track tests still pass after the helper refactor,
the `is_ebmvfr = FALSE` plumbing is correct by construction
(end-to-end IBSS test).

## D10. Out of scope (see proposal §Out of scope)

- `thresh_lowcount > 0` and `quantile_trans = TRUE`
- ELBO branch
- Full EB-multivariate-functional-regression model (the
  original out-of-scope EBmvFR ruling stands for SuSiE-
  alternative use)

## D11. Diagram

See `design/covariate-adjustment-flow.excalidraw` for the
data-flow diagram. Render via the official Excalidraw MCP in
Claude Desktop, or open at https://excalidraw.com.
