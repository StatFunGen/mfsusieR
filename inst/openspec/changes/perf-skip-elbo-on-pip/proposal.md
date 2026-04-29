# Skip ELBO refresh on PIP convergence path

## Why

mvf vs mfsusieR per-iter benchmark, default config (`prior =
"mixture_normal"` ↔ `prior_variance_scope = "per_outcome"`,
PIP convergence) on n=100, 0/1/2 genotype, M=2 functional outcomes:

| p   | mvf s/iter | mfsusieR s/iter | gap         |
|-----|-----------:|----------------:|------------:|
| 100 |      0.377 |           0.431 | mvf 14% faster |
| 200 |      0.549 |           0.744 | mvf 36% faster |
| 500 |      1.909 |           2.687 | mvf 41% faster |

The structural cause, traced from susieR's `susie_workhorse.R:87`:

```r
elbo[iter + 1] <- get_objective(data, params, model)   # called every iter
model <- check_convergence(data, params, model, elbo, iter)
```

`get_objective` is invoked unconditionally per IBSS iter. Our
override `get_objective.mfsusie` calls
`refresh_lbf_kl.mf_individual`, which performs a **full L-effect
SER stat sweep** (`compute_residuals + compute_ser_statistics +
loglik + compute_kl` for every l) so the reported ELBO is a
coherent variational free energy rather than the susieR-default
hybrid `Eloglik(state_t) - sum_l KL_l(state-when-l-was-updated)`.

That refresh structurally mirrors the IBSS sweep itself, so we
double the per-iter work. mvf and fsusieR both default
`cal_obj = FALSE` and skip ELBO computation entirely; their PIP
convergence reads `delta_alpha`, never the ELBO.

On the **PIP convergence path** (mfsusie's default
`convergence_method = "pip"`), the ELBO is never read by
`check_convergence`. The refresh is wasted work.

## What changes

`get_objective.mfsusie(data, params, model)` SHALL return
`NA_real_` when `params$convergence_method == "pip"`, skipping
`refresh_lbf_kl.mf_individual()` and `Eloglik()`. The full
coherent-ELBO path runs only when `convergence_method = "elbo"`.

Concretely:
```r
get_objective.mfsusie <- function(data, params, model) {
  if (identical(params$convergence_method, "pip")) {
    return(NA_real_)
  }
  model <- refresh_lbf_kl.mf_individual(data, params, model)
  Eloglik(data, model) - sum(model$KL, na.rm = TRUE)
}
```

susieR's `check_convergence` accepts `NA` ELBO without error on
the PIP path (it reads `delta_alpha`, not the ELBO array). The
returned `fit$elbo` array will contain `NA` entries on the PIP
path; documented in the roxygen for `convergence_method`.

## Acceptance criteria

* On the n=84, p=3500, M=6, T=128 0/1/2-genotype fixture, default
  config, per-iter wall time decreases **≥ 30%** vs the
  pre-change baseline.
* On the n=100, p∈{500, 1000, 2000}, M=2, T=128 fixture, default
  config, mfsusieR is **strictly faster per iter than
  `mvf.susie.alpha::multfsusie(prior = "mixture_normal")`** with
  matched `max_SNP_EM = p`, `greedy = FALSE`, `backfit = FALSE`,
  matched seed, hash-verified inputs.
* Numerical fits identical at all sizes (PIPs, alpha, mu, sigma2)
  at `tol = 1e-12` to the pre-change implementation.
* When `convergence_method = "elbo"`, ELBO trajectory unchanged
  from current behavior (still coherent variational free energy).
* All 1326 existing tests pass; new tests cover the PIP-path gate
  and the apples-to-apples mvf comparison.

## Impact

* Public API unchanged.
* `fit$elbo` may contain `NA` values on the PIP path. Documented.
* Negligible code change (~5 LOC).
* Removes a known per-iter cost on the default convergence path
  with no algorithmic compromise.

## Out of scope

* Refactoring susieR's workhorse to skip `get_objective` entirely
  (would require an upstream PR to susieR).
* Changing the default `convergence_method`.
* Other `refresh_lbf_kl` users (none currently — only
  `get_objective.mfsusie`).

## Risk + mitigation

* User scripts that compare ELBO trajectories without setting
  `convergence_method = "elbo"` will receive `NA`s. Mitigation:
  document the gate prominently; existing tests don't rely on
  ELBO under PIP convergence.
* `check_convergence`'s NA handling: verified by reading
  susieR's source — `check_convergence` short-circuits to
  `delta_alpha` on the PIP path before ELBO is read.
