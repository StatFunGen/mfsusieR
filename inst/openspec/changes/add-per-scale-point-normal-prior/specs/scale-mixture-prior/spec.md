# Per-scale point-normal and point-laplace priors (delta)

## ADDED Requirements

### Requirement: `prior_variance_scope = "per_scale_normal"` SHALL be a public option

`mfsusie()` SHALL accept `prior_variance_scope = "per_scale_normal"` as a public scope value that selects a per-(outcome, scale) point-normal prior on the per-effect coefficient `b_lm`.

The model is

```
b_lm[t in scale s] ~ pi_0_ms * delta_0  +  (1 - pi_0_ms) * N(0, sigma_ms^2)
```

with two parameters per (outcome, scale) cell estimated by empirical Bayes via `ebnm::ebnm_point_normal()` on the lead-variable observation slice.

#### Scenario: enum dispatch produces the new prior class

- **WHEN** `mfsusie()` is called with `prior_variance_scope = "per_scale_normal"`
- **THEN** `class(fit$G_prior[[m]])` MUST equal `"mixture_point_normal_per_scale"` for every outcome `m`, and each `fit$G_prior[[m]][[s]]$fitted_g` MUST be a 2-component normalmix with `pi` of length 2, `sd = c(0, sigma_ms)`, and `mean = c(0, 0)`.

### Requirement: `prior_variance_scope = "per_scale_laplace"` SHALL be a public option

`mfsusie()` SHALL accept `prior_variance_scope = "per_scale_laplace"` as a public scope value that selects a per-(outcome, scale) point-laplace prior on the per-effect coefficient `b_lm`.

The model is

```
b_lm[t in scale s] ~ pi_0_ms * delta_0  +  (1 - pi_0_ms) * Laplace(0, scale_ms)
```

with two parameters per (outcome, scale) cell estimated by empirical Bayes via `ebnm::ebnm_point_laplace()` on the lead-variable observation slice.

#### Scenario: enum dispatch produces the laplace prior class

- **WHEN** `mfsusie()` is called with `prior_variance_scope = "per_scale_laplace"`
- **THEN** `class(fit$G_prior[[m]])` MUST equal `"mixture_point_laplace_per_scale"` for every outcome `m`, and each `fit$G_prior[[m]][[s]]$fitted_g` MUST be a `laplacemix` record (`pi`, `scale`, `mean` of length 2 with the spike component at 0).

### Requirement: M-step on the ebnm-backed prior classes SHALL call ebnm on the alpha-thinned multi-variable rectangle

The M-step on prior class `"mixture_point_normal_per_scale"` SHALL call `ebnm::ebnm_point_normal()`, and the M-step on `"mixture_point_laplace_per_scale"` SHALL call `ebnm::ebnm_point_laplace()`, per (outcome, scale) on the flattened alpha-thinned rectangle of (Bhat, Shat) with `g_init` set to the previous IBSS iter's `fitted_g`.

The data passed to ebnm is `(x = as.vector(bhat_m[keep_idx, idx_s]), s = as.vector(shat_m[keep_idx, idx_s]))` where `keep_idx` is the alpha-thinned variable indices computed in `optimize_prior_variance.mf_individual()` using `params$alpha_thin_eps` (the same threshold the mixsqp arm uses). The returned `fit$fitted_g` SHALL be written into `G_m[[s]]$fitted_g` without modification, and `fit$fitted_g$pi` SHALL be written into `model$pi_V[[m]][s, ]`. `mixture_null_weight` SHALL be silently ignored on these prior classes.

#### Scenario: ebnm receives the multi-variable rectangle with warm-start g_init

- **WHEN** the M-step runs at (l, m, s) with prior class `"mixture_point_normal_per_scale"` or `"mixture_point_laplace_per_scale"`, and the previous IBSS iter wrote a `fitted_g` into `G_m[[s]]`
- **THEN** the wrapper MUST call `ebnm_fn(x = as.vector(bhat_m[keep_idx, idx_s]), s = as.vector(shat_m[keep_idx, idx_s]), g_init = G_m[[s]]$fitted_g, fix_g = !isTRUE(params$estimate_prior_variance))`, MUST overwrite `G_m[[s]]$fitted_g` with `fit$fitted_g`, and MUST overwrite `model$pi_V[[m]][s, ]` with `fit$fitted_g$pi`.

#### Scenario: M-step is short-circuited when `estimate_prior_variance = FALSE`

- **WHEN** the M-step is invoked under `estimate_prior_variance = FALSE` for either ebnm-backed prior class
- **THEN** the M-step MUST forward `fix_g = TRUE` to ebnm; the post-call `fitted_g` MUST equal the pre-call `fitted_g` at `tol = 1e-12`.

#### Scenario: M-step is idempotent on identical inputs

- **WHEN** the M-step on either ebnm-backed class is called twice in succession with identical `(bhat_m, shat_m, keep_idx, fitted_g)` inputs
- **THEN** the two outputs MUST be bit-identical at `tol = 1e-12` (ebnm's optim is deterministic on identical inputs).

### Requirement: Init step SHALL pick the per-scale lead from marginal data

`init_ebnm_prior_per_scale()` SHALL pick the per-scale lead variable from the marginal data, independent of `alpha`, before fitting the initial prior with ebnm.

For each scale `s` with index set `idx_s`, the helper computes `bs <- compute_marginal_bhat_shat(X, Y_m)` once and sets `lead_s <- which.max(rowMeans(bs$Bhat[, idx_s, drop = FALSE]^2))`. ebnm then fits the prior on `(bs$Bhat[lead_s, idx_s], bs$Shat[lead_s, idx_s])` and the helper records `lead_init_s = lead_s` on `G_prior[[s]]` for downstream tests.

#### Scenario: lead picker selects a signal-bearing variable on the sparse fixture

- **WHEN** `init_ebnm_prior_per_scale()` runs on the synthetic sparse fixture (`set.seed(2L); n = 200; p = 30; T_m = 64; signal at variables c(7, 18) and times c(20, 44)`) for either ebnm-backed prior class
- **THEN** for every scale `s` whose `idx_s` overlaps the signal positions `c(20L, 44L)`, the recorded `G_prior[[1L]][[s]]$lead_init_s` MUST be in `c(7L, 18L)`.

### Requirement: Cache refresh SHALL skip the mixsqp-only blocks on the ebnm paths

`refresh_iter_cache.mf_individual()` SHALL skip the mixsqp-only `sdmat` and `log_sdmat` builds when `model$G_prior[[1L]]` does not inherit from `"mixture_normal"` or `"mixture_normal_per_scale"`, while still building `shat2` so `mf_per_outcome_bhat_shat()` does not fall through to its per-call recomputation.

#### Scenario: ebnm path skips sdmat / log_sdmat but keeps shat2

- **WHEN** `refresh_iter_cache.mf_individual()` runs with `model$G_prior` tagged as `"mixture_point_normal_per_scale"` or `"mixture_point_laplace_per_scale"`
- **THEN** `model$iter_cache$sdmat` and `model$iter_cache$log_sdmat` MUST be NULL, and `model$iter_cache$shat2` MUST be set (one `p × T_basis[m]` matrix per outcome `m`).

### Requirement: `per_scale_normal` SHALL bit-match `susie()` at the degenerate case

`mfsusie(prior_variance_scope = "per_scale_normal")` SHALL produce a fit numerically identical to `susieR::susie()` at `tol = 1e-12` when invoked at the susie-degenerate parameter point: scalar `Y` (`T_m = 1`), zero null-mass init (`null_prior_init = 0`), a fixed slab variance (`sigma_init`), and `estimate_prior_variance = FALSE` (so ebnm's M-step is short-circuited via `fix_g = TRUE` and the prior stays at the fixed Gaussian).

#### Scenario: scalar `Y`, fixed prior, no EB matches susie

- **WHEN** the same fixture (`X`, `y`, seed) is fit with both
  ```r
  susieR::susie(X, y, L = 5, scaled_prior_variance = sigma2,
                estimate_prior_variance = FALSE,
                estimate_residual_variance = TRUE,
                max_iter = 100, tol = 1e-8)
  ```
  and
  ```r
  mfsusie(X, list(matrix(y, ncol = 1)), L = 5,
          prior_variance_scope = "per_scale_normal",
          null_prior_init = 0, sigma_init = sqrt(sigma2),
          estimate_prior_variance = FALSE,
          residual_variance_scope = "per_outcome",
          L_greedy = NULL, max_iter = 100, tol = 1e-8,
          verbose = FALSE)
  ```
- **THEN** the two fits MUST agree on `alpha`, `pip`, `sigma2[[1]]`, `lbf`, `KL`, `tail(elbo, 1)`, `niter`, and `lapply(sets$cs, sort)` at `tol = 1e-12`.
