# mfsusieR vs mvf.susie.alpha — three numerical-formula divergences

Date: 2026-04-29
Scope: consolidate every formula-level divergence between mfsusieR
and `mvf.susie.alpha` found across prior audits and this session's
head-to-head bench. Drives Phase 5 FDR investigation and informs
the per-iter perf gap discussion.

## Summary table

| # | Where in mvf | Where in mfsusieR | Type | Status |
|---|---|---|---|---|
| 1 | `R/operation_on_multfsusie_obj.R:1075-1090` `get_ER2.multfsusie` | `R/individual_data_methods.R::mf_get_ER2_per_position` | **mvf bug** (Pattern A fix) | mvf formula misses the `predictor_weights` factor in the bias correction and aggregates effects via sum-then-square. mfsusieR has the correct decomposition. Documented in `refactor-exceptions.md` PR group 6 (lines 262-287). Pattern A regression guard: `tests/testthat/test_mvf_alpha_degeneracy.R:114-129`. |
| 2 | `R/ELBO_mutlfsusie.R:269-282` `get_objective.multfsusie` | `R/individual_data_methods.R::get_objective.mfsusie` | **mvf bug** (sign error on `loglik_SFR_post`); the added `o` entropy term is a partial workaround that does not fully correct it | mvf has `KL_l = -loglik_SFR(l) - loglik_SFR_post(l)` where the second term should enter with `+` per the susieR derivation. The `o = sum_l[log(p) + H(alpha_l)]` term is bounded by `L*log(p)` and cannot offset the sign error (which is of order `L*n*T`). Full derivation in `inst/notes/cross-package-audit-derivations.md` section 2. mfsusieR follows the correct susieR decomposition: `Eloglik - sum(KL)`. |
| 3 | `fsusieR::cal_Bhat_Shat` (called by mvf via `cal_Bhat_Shat_multfsusie`) | `R/individual_data_methods.R::compute_ser_statistics.mf_individual` + `mf_per_outcome_bhat_shat` | **Different SE estimator**, not strictly a bug | mvf uses a per-(variable, outcome) marginal residual variance: `Shat[j,m]^2 = colVars(Y[,m] - X[,j] * Bhat[j,m]) / (n-1)`. Each variable carries its own `sigma_jm^2` from a separate univariate regression. mfsusieR uses the SuSiE backbone form: a single `sigma2[m]` per outcome (or `sigma2[m,s]` per scale), then `Shat[j,m]^2 = sigma2[m] / colSums(X^2)[j]`. New finding 2026-04-29; ledger entry pending. |

The two bugs are documented and fixed. The third is a model-formulation
choice that interacts with bug (1) in ways the prior audit did not
quantify; this session uncovered the interaction.

## How they interact (this session's new finding)

### Bench result on identical data (n=84, p=500, M=6, T=128, L=10, true variables {11, 200, 450})

| | mfsusieR | mvf |
|---|---|---|
| niter | 2 | 2 |
| PIP at variable 11 | 1.000 | 1.000 |
| PIP at variable 200 | 1.000 | 1.000 |
| PIP at variable 450 | **0.016** | **1.000** |
| sigma2 (per-outcome mean) | 0.72-0.86 | 0.54 |

True noise variance is 0.36 (data generated with `sd_noise = 0.6`).
Both estimates exceed truth because the unmodelled third signal
inflates residual variance, but mvf's estimate is closer.

### Mechanism

Direction of bug (1) on `sigma2`:

- mvf's ER2 bias correction is `-sum(postF^2) + sum(postF2_sd2)`.
  The corresponding correct term is
  `sum_l((alpha_l . pw)^T mu2_l - ||X(alpha_l mu_l)||^2)`.
- The bug is missing the `predictor_weights = colSums(X^2) ≈ n-1`
  factor on each side. mvf's bias correction is roughly the correct
  one divided by `n`, so it is much smaller than the correct value.
- The correct bias correction is positive (it captures the
  variational uncertainty in `beta`). With the bug, mvf's ER2 is
  approximately the raw residual sum of squares with ~no bias
  correction, while ours adds the full positive correction.
- Therefore mvf's `sigma2 = ER2 / (n*T)` is **smaller** than ours.

Direction of estimator (3) on per-variable BFs:

- For null variables, `colVars(Y - X[,j] * Bhat[,j]) ≈ Var(Y) ≈ sigma2`,
  so the two SE formulas converge.
- For true-signal variables, regressing `Y` on `X[,j]` soaks up the
  signal correlated with `X[,j]`, so
  `colVars(Y - X[,j] * Bhat[,j]) < sigma2`. mvf's per-(j) `Shat^2`
  is therefore **smaller at signal variables** than the SuSiE-style
  joint `Shat^2 = sigma2 / colSums(X^2)`.
- Smaller Shat means larger Bayes factor at signal variables, more
  aggressive detection.

Both effects push mvf toward higher detection power on signal
variables, which is what we observed on this dataset (mvf detects
variable 450; we miss it).

The bug (1) and the estimator (3) compound: bug (1) gives mvf a
smaller `sigma2` than ours; estimator (3) further reduces the SE
at signal variables relative to the joint SuSiE form. From a
point-estimate perspective both push the BF up.

### What this implies

- The "bug helps detection" reading is consistent with Anjing's
  earlier observation that mvf has more PIP mass at non-causal
  variants under null simulations (FDR miscalibration).
  Pattern: an under-estimated `sigma2` inflates BFs uniformly,
  which inflates PIPs at both signal and null variants. The
  apparent power gain on this seed is genuine on the signal
  variable but is the same mechanism that causes the FDR
  miscalibration on null seeds.
- The PIP convergence path (`pip_stall_window = 5`) declares done
  when no PIP changes for 5 iters. Our IBSS at iter 2 has not yet
  given variable 450 enough alpha mass to move; PIP convergence
  fires; we stop. With `convergence_method = "elbo"` and a tighter
  tol we may pick up variable 450 over more iters as `sigma2`
  converges. Worth testing before any conclusion about detection
  power.

## What is verified (read source)

### (1) ER2 — exact lines

mvf, `R/operation_on_multfsusie_obj.R:1046-1096` `get_ER2.multfsusie`:
the bias correction is computed as `sum(postF2_sd2) - sum(postF^2)`
where `postF = sum_l alpha_l * mu_l`. Two errors:

- `sum(postF^2)` aggregates effects via sum-then-square. The correct
  per-effect form is `sum_l ||X * (alpha_l * mu_l)||^2` (square-then-sum).
- The sums of squared coefficients are not weighted by
  `predictor_weights = colSums(X^2)`, which carry an `O(n)` factor
  for scaled X.

mfsusieR, `R/individual_data_methods.R::mf_get_ER2_per_position`
(line 417 in current source):
```r
res_m <- (data$D[[m]] - model$fitted[[m]])[idx_m, , drop = FALSE]
rss_t <- colSums(res_m * res_m)
bias_t <- numeric(T_pad)
for (l in seq_len(model$L)) {
  alpha_l  <- model$alpha[l, ]
  mu_lm    <- model$mu[[l]][[m]]
  mu2_l_m  <- model$mu2[[l]][[m]]
  XEb_lm   <- (X %*% (alpha_l * mu_lm))[idx_m, , drop = FALSE]
  bias_t   <- bias_t + colSums((alpha_l * pw) * mu2_l_m) -
              colSums(XEb_lm * XEb_lm)
}
rss_t + bias_t
```
This is the susieR `get_ER2.individual` decomposition generalized to
per-(modality, position).

### (2) ELBO `o` term — exact lines

mvf, `R/ELBO_mutlfsusie.R:269-282`: returns `Eloglik - sum(KL) + o`,
where `o` aggregates per-effect categorical entropy.

`mvf::cal_KL_l` returns `KL_l = -loglik_SFR(l) - loglik_SFR_post(l)`.
Per the susieR derivation:
```
KL_l_correct  = -lbf_l - L_null + E_q[log p(Y|b_l)]
KL_l_mvf      = -lbf_l - L_null - E_q[log p(Y|b_l)]
KL_l_mvf - KL_l_correct = -2 * E_q[log p(Y|b_l)]
```
`E_q[log p(Y|b_l)]` is the SER posterior expected log-likelihood,
typically a large negative number (`O(L * n * T)`). The bug inflates
the sum of `KL_l` by `O(L * n * T)`. The `o` correction is bounded
by `L * log p`, much smaller. So `o` does not offset the bug.

mfsusieR `R/individual_data_methods.R::get_objective.mfsusie`
(line 745 in current source):
```r
get_objective.mfsusie <- function(data, params, model) {
  if (identical(params$convergence_method, "pip")) return(NA_real_)
  model <- refresh_lbf_kl.mf_individual(data, params, model)
  Eloglik(data, model) - sum(model$KL, na.rm = TRUE)
}
```
follows the correct susieR decomposition (no `o` term, correct sign
on the SER posterior expected log-likelihood inside `compute_kl`).

### (3) Bhat/Shat estimator — exact lines

mvf calls `fsusieR::cal_Bhat_Shat` via `cal_Bhat_Shat_multfsusie`. The
fsusieR helper (pure R):
```r
Bhat <- crossprod(X, Y) / d                       # J x T
Shat <- do.call(cbind, lapply(seq_len(J),
    function(j) Rfast::colVars(Y[,j] - sweep(X, 2, Bhat[,j], "*"))))
Shat <- sqrt(pmax(Shat, 1e-64)) / sqrt(n - 1)
```
Each `(j, m)` cell of `Shat` is the SE of the marginal regression
of `Y[,m]` on `X[,j]`.

mfsusieR `R/individual_data_methods.R::compute_ser_statistics.mf_individual`
+ `mf_per_outcome_bhat_shat`:
```r
bhat_m  <- XtR_m / pw                       # XtR_m = X^T R_m
shat2_m <- outer(1 / pw, sigma2_per_pos)    # cached on em_cache
```
where `sigma2_per_pos` is broadcast from `model$sigma2[[m]]` (one
scalar per outcome, or one per scale), updated each iter via
`update_variance_components`. This is the SuSiE backbone form:
single `sigma2` per outcome shared across variables, with
`Shat^2[j] = sigma2 / sum(x_j^2)`.

The susieR fork also exposes `compute_marginal_bhat_shat` at
`~/GIT/susieR/R/univariate_regression.R:207`, which packages both
estimators in one call:
```r
if (!is.null(sigma2)) {
  Shat <- matrix(sqrt(sigma2 / predictor_weights), nrow = J, ncol = T_y)
} else {
  Shat <- vapply(seq_len(T_y),
    function(t) Rfast::colVars(Y[,t] - sweep(X, 2, Bhat[,t], "*")),
    numeric(J))
  Shat <- sqrt(pmax(Shat, 1e-64)) / sqrt(n - 1)
}
```
Calling `compute_marginal_bhat_shat(X, R_m, sigma2 = sigma2_per_pos)`
reproduces our current path; calling without `sigma2` reproduces
fsusieR's path.

## Per-iter perf gap context (this session)

Head-to-head at identical config (`prior_variance_scope = "per_outcome"`
on our side, `prior = "mixture_normal"` on mvf):

| p    | mfsusieR (s/iter) | mvf (s/iter) | gap |
|------|-------------------|--------------|-----|
| 200  | 8.14              | 9.69         | 0.84x (we win 16%) |
| 500  | 18.34             | 17.77        | 1.03x |
| 1000 | 37.46             | 35.56        | 1.05x |

Per-iter wall is comparable. mvf wastes ~45% inclusive of total time
on `irlba` (mixsqp's `tsvd` preconditioner, triggered by leaving
`tol.svd` at the default `1e-6`) but compensates by avoiding our
three-pass dnorm structure. We set `tol.svd = 0` in
`R/em_helpers.R:144` because our `(p*T+1) x K` L matrix is
tall-skinny, so no SVD truncation can occur. Per-effect call count
is identical: M mixsqp calls + 3 dnorm matrix passes.

Per-iter compute is therefore not what makes mvf detect signal we
miss. The detection difference is driven by divergences (1) and (3).

## VI consistency

From the variational lower bound:
- The ELBO depends on a single `sigma2` (per outcome / per scale).
- The per-variable BF is derived from the same `sigma2`.
- `(A)` SuSiE / mfsusieR plug the same `sigma2` into `Shat^2 =
  sigma2 / sum(x_j^2)`, so ELBO and BF are coherent.
- `(B)` fsusieR / mvf plug a per-(j, m) `sigma_jm^2` into `Shat`,
  not the ELBO's `sigma2`. The BF and the ELBO are no longer
  derived from one variational posterior; the resulting "ELBO"
  is not a clean lower bound.

**(A) is correct from a variational standpoint.** (B) is a
heuristic that was the de facto standard in fsusieR, inherited by
mvf. The two formulations agree at null variants and diverge at
signal variants in the direction described above.

## Open questions

1. Does our PIP convergence (`pip_stall_window = 5`) stop too early
   on this dataset? Test: rerun with `convergence_method = "elbo"`,
   `tol = 1e-6`, `max_iter = 50` and check whether PIP at variable
   450 climbs.
2. Is our higher `sigma2` causing the under-detection at variable
   450 directly, or is it the joint-vs-marginal SE form? Test: a
   manual fit that uses our IBSS but plugs in fsusieR-style per-(j,
   m) `Shat` (or vice versa) on a known seed.
3. Should `compute_ser_statistics.mf_individual` add an opt-in to
   the per-(j, m) marginal `Shat` for power on real data, with the
   default kept on the variationally-consistent joint `Shat`? This
   would be a Phase 5 / Phase 6 question, not a perf change.

## File pointers

- ER2 derivation (manuscript-grounded): `R/individual_data_methods.R::mf_get_ER2_per_position`
- ER2 ledger: `inst/notes/refactor-exceptions.md` PR group 6
- ELBO `o` derivation: `inst/notes/cross-package-audit-derivations.md` section 2
- ELBO ledger reference: `inst/notes/sessions/2026-04-26-cross-package-audit-c.md` item 4.4
- cal_Bhat_Shat finding: this note (new)
- Bench scripts: `inst/notes/sessions/scripts/bench_vs_mvf.R`,
  `inst/notes/sessions/scripts/compare_fits.R`
- Profile dumps used to verify mvf irlba path: `/tmp/prof_mfsusie.out`,
  `/tmp/prof_mvf.out` (one-shot, not committed)
