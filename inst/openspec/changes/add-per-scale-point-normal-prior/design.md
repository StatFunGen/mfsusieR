# Design: per-scale point-normal and point-laplace priors

## Architecture

Four changes; the rest of the IBSS workhorse is untouched.

The shape mirrors how `mixsqp` is wired today: an `init_*` helper
seeds `G_prior[[m]][[s]]$fitted_g` at IBSS iter 0; an `.opv_<class>`
helper runs the per-(m, s) M-step at every IBSS iter, warm-started
from the previous iter's `fitted_g`. ebnm replaces `mixsqp` inside
the per-(m, s) call; the surrounding scaffold is the same.

### 1. `prior_variance_scope` enum extension

Public API:
```r
prior_variance_scope = c("per_outcome", "per_scale",
                         "per_scale_normal", "per_scale_laplace")
```

Internal `prior_class` mapping (set in `mf_prior_scale_mixture`):

```r
prior_class <- switch(prior_variance_scope,
  per_outcome       = "mixture_normal",
  per_scale         = "mixture_normal_per_scale",
  per_scale_normal  = "mixture_point_normal_per_scale",
  per_scale_laplace = "mixture_point_laplace_per_scale"
)
```

`prior_class` is already an attribute of `G_prior`
(`R/prior_scale_mixture.R:124`). The two new values carry the
same shape (list of K-component fitted-g records, one per scale);
K = 2 for both new classes.

### 2. Init step, marginal-data lead per scale + ebnm

New helper `init_ebnm_prior_per_scale(Y_m, X, scale_index,
prior_class, ...)` builds G_prior with one ebnm fit per
(outcome, scale):

```r
init_ebnm_prior_per_scale <- function(Y_m, X, scale_index,
                                       prior_class, ...) {
  ebnm_fn <- switch(prior_class,
    mixture_point_normal_per_scale  = ebnm::ebnm_point_normal,
    mixture_point_laplace_per_scale = ebnm::ebnm_point_laplace
  )
  bs <- compute_marginal_bhat_shat(X, Y_m)
  G_prior <- vector("list", length(scale_index))
  for (s in seq_along(scale_index)) {
    idx <- scale_index[[s]]
    lead_s <- which.max(rowMeans(bs$Bhat[, idx, drop = FALSE]^2))
    fit <- ebnm_fn(x = bs$Bhat[lead_s, idx],
                   s = bs$Shat[lead_s, idx])
    G_prior[[s]] <- list(
      fitted_g    = fit$fitted_g,
      idx         = idx,
      lead_init_s = lead_s          # recorded for tests
    )
  }
  class(G_prior) <- prior_class
  G_prior
}
```

The lead picker is the deterministic Option A: the variable with
the largest mean squared marginal Bhat across scale `s`'s
positions. This is well-defined at IBSS iter 0 (no `alpha`
dependence; no tie-break trap when `alpha` is uniform `1/p`) and
gives ebnm a real signal-tracking observation set rather than the
moment estimate over the full `(p, idx_size)` rectangle.

`mf_prior_scale_mixture` dispatches on `prior_variance_scope`:

```r
prior_obj <- if (prior_variance_scope %in%
                 c("per_scale_normal", "per_scale_laplace")) {
  init_ebnm_prior_per_scale(Y_m, X, scale_index,
                            prior_class = prior_class, ...)
} else {
  init_scale_mixture_prior_default(...)   # current path
}
```

### 3a. M-step body, `ebnm_fn(x, s, g_init, fix_g)` on the multi-variable rectangle

The M-step body for either ebnm-backed class flattens the
alpha-thinned `(|keep_idx| × |idx_s|)` rectangle of (Bhat, Shat)
to two vectors and hands them to ebnm:

```r
fit <- ebnm_fn(x      = as.vector(bhat_m[keep_idx, idx, drop = FALSE]),
               s      = as.vector(shat_m[keep_idx, idx, drop = FALSE]),
               g_init = G_m[[s]]$fitted_g,
               fix_g  = !isTRUE(params$estimate_prior_variance))
G_m[[s]]$fitted_g  <- fit$fitted_g
model$pi_V[[m]][s, ] <- fit$fitted_g$pi
```

The data shape ebnm sees is the same alpha-thinned rectangle
mixsqp consumes, just unweighted (ebnm cannot accept observation
weights). The unweighted-many-variables design
mirrors mixsqp's structural property: ebnm sees the bulk-noise
distribution within `keep_idx` and the parametric MLE settles on
a sensible conservative `pi_0` (most kept variables look null)
plus a `sigma` driven by the few non-null cells. The lead-only
alternative would fit "this single variable is non-null" and
collapse `pi_0` toward 0, producing spurious CSes from the
too-liberal slab.

`g_init = <previous fit>` is the warm-start path: ebnm seeds
optim from the previous IBSS iter's `fitted_g`. This mirrors
mixsqp's `pi_warm_start = pi_prev` (always passed; no cold/warm
branch). `fix_g = TRUE` collapses the M-step to a no-op so the
prior stays at the value the init helper or the previous iter
wrote (used to honor `estimate_prior_variance = FALSE` and to
hold the prior fixed in the susie-degenerate locks).

`alpha_thin_eps` is the unified per-effect alpha threshold,
renamed from `mixsqp_alpha_eps`. The same threshold drives both
arms: drop variables with `alpha[l, j] < alpha_thin_eps` from
the M-step input. On iter 1 alpha is uniform `1/p` and all
variables are kept (giving ebnm the full bulk-noise
distribution); on later iters `keep_idx` shrinks to the
LD-friends of the lead (still a multi-variable set, so the
per-(m, s) prior cannot collapse to single-variable
over-fitting).

ebnm owns:
- optim transport and convergence (the L-BFGS-B / `ashr` internal
  solver).
- init heuristics when `g_init` is `NULL` (cold start at iter 0).
- the small-n freeze (`n <= 2`).
- the `optim` non-convergence fallback (returns `g_init`).
- the `s == 0`, `x == NaN` adversarial-input guards.
- the format of the returned `fitted_g` record.

### 3b. M-step dispatch, `.opv_<class>` helper mirroring `.opv_mixsqp`

The dispatch follows the existing `.opv_<class>` helper-per-class
pattern. The mixsqp branch is `.opv_mixsqp`
(`R/individual_data_methods.R:584-637`); the ebnm-backed branches
are `.opv_ebnm_point_normal` and `.opv_ebnm_point_laplace`. Each
helper takes the same arguments and returns the updated `model`:

```r
.opv_ebnm_point_normal <- function(data, params, model, ser_stats,
                                    keep_idx, zeta_keep) {
  .opv_ebnm_point(data, params, model, ser_stats,
                  keep_idx, zeta_keep,
                  ebnm_fn = ebnm::ebnm_point_normal)
}

.opv_ebnm_point_laplace <- function(data, params, model, ser_stats,
                                     keep_idx, zeta_keep) {
  .opv_ebnm_point(data, params, model, ser_stats,
                  keep_idx, zeta_keep,
                  ebnm_fn = ebnm::ebnm_point_laplace)
}

.opv_ebnm_point <- function(data, params, model, ser_stats,
                             keep_idx, zeta_keep, ebnm_fn) {
  fix_g <- !isTRUE(params$estimate_prior_variance)
  for (m in seq_len(data$M)) {
    bhat_m <- ser_stats$betahat[[m]][keep_idx, , drop = FALSE]
    shat_m <- sqrt(ser_stats$shat2[[m]][keep_idx, , drop = FALSE])
    G_m    <- model$G_prior[[m]]
    for (s in seq_along(G_m)) {
      idx <- G_m[[s]]$idx
      fit <- ebnm_fn(
        x      = as.vector(bhat_m[, idx, drop = FALSE]),
        s      = as.vector(shat_m[, idx, drop = FALSE]),
        g_init = G_m[[s]]$fitted_g,
        fix_g  = fix_g
      )
      model$G_prior[[m]][[s]]$fitted_g <- fit$fitted_g
      model$pi_V[[m]][s, ]              <- fit$fitted_g$pi
    }
  }
  model
}
```

Dispatch is the existing `if (inherits(...)) ... else if (...)`
arm in `optimize_prior_variance.mf_individual()`, extended with
two new arms:

```r
optimize_prior_variance.mf_individual <- function(data, params, model,
                                                   ser_stats, l, ...) {
  # ... existing keep_idx / zeta_keep computation ...
  G1 <- model$G_prior[[1L]]
  if (inherits(G1, c("mixture_normal", "mixture_normal_per_scale"))) {
    model <- .opv_mixsqp(data, params, model, ser_stats,
                         keep_idx, zeta_keep)
  } else if (inherits(G1, "mixture_point_normal_per_scale")) {
    model <- .opv_ebnm_point_normal(data, params, model, ser_stats,
                                     keep_idx, zeta_keep)
  } else if (inherits(G1, "mixture_point_laplace_per_scale")) {
    model <- .opv_ebnm_point_laplace(data, params, model, ser_stats,
                                      keep_idx, zeta_keep)
  } else {
    stop("Unknown prior class on G_prior[[1]]: ",
         paste(class(G1), collapse = ", "))
  }
  list(V = 1, model = model)
}
```

No new S3 generic. No parent-class hierarchy. Class tags stay
single. The dispatch is one `if/else if` arm per class, ~12
lines, the same shape as the existing branch.

With the foreseeable ebnm family being `point_normal` and
`point_laplace`, two helpers + two arms are simpler than an S3
generic + parent class + switch. A third ebnm family (e.g., a
point-t slab) would refactor then.

No deferred-M-step gate is needed. At IBSS iter 1 alpha is
uniform `1/p`, so `keep_idx` covers every variable and ebnm sees
the bulk-noise distribution; pi_0 fits close to 1. Once alpha
concentrates, `keep_idx` shrinks to the LD-friends of the lead
and ebnm refits on that smaller multi-variable rectangle. The
shared `alpha_thin_eps` thinning suffices; no `2/p` arbitrary
constant, no `greedy_lbf_cutoff` gate.

### 3c. `iter_cache` with class-gated slots

`refresh_iter_cache.mf_individual` (renamed from
`refresh_em_cache.mf_individual`; `em_cache` was misleading
because the cache is consumed by SER calls too, not just the
M-step) holds two kinds of per-IBSS-iter precompute:

- **`iter_cache$shat2[[m]]`** (`p × T_basis[m]`): always built.
  Equals `outer(1/xtx_diag, sigma2_per_pos)`. Read by
  `mf_per_outcome_bhat_shat()` on every SER call across both
  arms; saves one `outer()` per (l, m) per IBSS iter (~`L × M`
  outer products per iter).
- **`iter_cache$sdmat[[m]][[s]]`, `iter_cache$log_sdmat[[m]][[s]]`**
  (`(p × |idx_s|) × K` per (m, s)): built only when the prior
  class is mixsqp-backed (`mixture_normal` or
  `mixture_normal_per_scale`). The K-axis precompute is needed
  by `mf_em_likelihood_per_scale` from inside `.opv_mixsqp`;
  ebnm fits 2 scalars per (m, s) by direct optim, so there is
  no K-axis aggregate to amortize.

The class gate is a one-line boolean:

```r
refresh_iter_cache.mf_individual <- function(data, model) {
  is_mixsqp_prior <- inherits(
    model$G_prior[[1L]],
    c("mixture_normal", "mixture_normal_per_scale"))
  # ... build shat2 always ...
  if (!is_mixsqp_prior) {
    return(model)   # skip the mixsqp-only sdmat / log_sdmat builds
  }
  # ... build sdmat / log_sdmat ...
}
```

The previous `sigma2_per_pos[[m]]` slot was dead storage (no
consumer read it; every reader called `mf_sigma2_per_position()`
directly). Dropped.

Naming: `is_mixsqp_prior` reads as "this prior class needs the
K-axis precompute slots." There is no ebnm-specific iter_cache
slot: ebnm has no per-(m, s) precomputable aggregate; its
per-(m, s) state is `g_init`, carried on
`G_prior[[m]][[s]]$fitted_g`.

### 3d. Covariate-adjustment path stays unchanged

`mf_adjust_for_covariates(method = "wavelet_eb")` calls
`mf_em_likelihood_per_scale` and `mf_em_m_step_per_scale` directly
with `is_ebmvfr = TRUE` (`R/adjust_covariates.R:299-318`). It is
not routed through `optimize_prior_variance.mf_individual` and
therefore not through the new dispatch arms. The covariate-adjust
path is not an IBSS effect, does not have an `alpha`, and does
not share the per-scale-point-prior use case. It stays on the
existing mixsqp-with-`is_ebmvfr=TRUE` path.

A regression test confirms `mf_adjust_for_covariates(method =
"wavelet_eb")` continues to work unchanged after this change
lands.

### 4. Loglik / posterior moments

`loglik.mf_individual` and
`calculate_posterior_moments.mf_individual` already iterate over
`G_m` groups and call `mixture_log_bf_per_scale(...)` on each
group's `fitted_g`. The kernel reads `pi` and `sd` and is
K-agnostic, so the Normal path runs unchanged
(`R/individual_data_methods.R:184-198, 284-...`).

The Laplace path needs one small addition: a posterior-moment
helper that handles the Laplace slab. ebnm exposes
`ebnm::ebnm_point_laplace`'s posterior moments via
`ebnm_object$posterior` (mean, sd2), which we map into the
existing `(post_mean[m, j], post_var[m, j])` layout in
`calculate_posterior_moments.mf_individual`. The wiring is a
single `if (inherits(G_m, "mixture_point_laplace_per_scale"))`
branch around the existing normalmix moment formula.

## Lead variable at init

The marginal-data lead picker only acts at the init helper; the
IBSS-loop M-step uses every variable in `keep_idx`, not a
single lead.

**Init (IBSS iter 0).** Lead per scale is picked from the
marginal data:
```r
bs <- compute_marginal_bhat_shat(X, Y_m)
lead_s = which.max(rowMeans(bs$Bhat[, idx_s, drop = FALSE]^2))
```
This is deterministic, alpha-independent, and selects a
variable whose marginal magnitude across scale `s`'s positions
is largest. ebnm fits the prior on
`(bs$Bhat[lead_s, idx_s], bs$Shat[lead_s, idx_s])` and stores
`fitted_g` in `G_prior[[m]][[s]]$fitted_g`. The init prior gives
the IBSS effect-1 SER step a real signal-tracking starting prior
rather than a moment-estimator default; once `alpha` accumulates
mass through the SER step, the M-step refits ebnm on the
multi-variable rectangle.

The recorded `lead_init_s` is kept on the G_prior entry as a
test artifact (asserts the picker selected a signal-bearing
variable on known fixtures); it is not consumed by the
IBSS loop.

Code uses **lead variable**; vignette prose uses **lead
variant**.

## Degeneracy tests at machine precision

Three pairs of fits, all bit-equivalent to `susieR::susie()` at
`tol = 1e-12`. Common parameter pattern: T=1 (scalar Y), single
non-null variance component, no null mass, no EB:

```r
y_scalar <- ...
Y_list   <- list(matrix(y_scalar, ncol = 1))
sigma2   <- 0.2
common <- list(L = 5, residual_variance_scope = "per_outcome",
               estimate_prior_variance = FALSE, L_greedy = NULL,
               max_iter = 100, tol = 1e-8, verbose = FALSE)

fit_s   <- susie(X, y_scalar, L = 5, scaled_prior_variance = sigma2,
                 estimate_prior_variance = FALSE,
                 estimate_residual_variance = TRUE,
                 max_iter = 100, tol = 1e-8)

# (1) per_outcome: already exists
fit_d_o <- mfsusie(X, Y_list, prior_variance_scope = "per_outcome",
                   prior_variance_grid = sigma2,
                   null_prior_init = 0, !!!common)

# (2) per_scale: NEW
fit_d_s <- mfsusie(X, Y_list, prior_variance_scope = "per_scale",
                   prior_variance_grid = sigma2,
                   null_prior_init = 0, !!!common)

# (3) per_scale_normal: NEW, primary correctness lock
fit_d_p <- mfsusie(X, Y_list, prior_variance_scope = "per_scale_normal",
                   sigma_init = sqrt(sigma2),
                   null_prior_init = 0,        # pi_0 = 0
                   estimate_prior_variance = FALSE,
                   !!!common)

# All three must match fit_s at tol = 1e-12 on alpha, pip, mu,
# sigma2, lbf, KL, elbo, niter, sets$cs.
```

The `per_scale_normal` degenerate setup: `pi_0 = 0` collapses the
prior to `N(0, sigma^2)`. With `estimate_prior_variance = FALSE`
the M-step is short-circuited (we forward `fix_g = TRUE` to ebnm,
or skip the call entirely; the contract is no movement on
`fitted_g`). The prior is held at that single Gaussian throughout,
which is exactly classic SuSiE.

The `per_scale_laplace` path does not enter the susie-degeneracy
locks: a Laplace slab is not equivalent to classic SuSiE's single
Gaussian. Its degeneracy contract is internal: `per_scale_laplace`
with `pi_0 = 0` and a fixed scale produces the same fit on two
back-to-back runs (idempotence at `tol = 1e-12`).

## Failure modes / tradeoffs

1. **Coarse-scale data scarcity (`idx_size <= 2`).** ebnm's
   small-n behavior is documented: at `n <= 2` it freezes the
   prior at `g_init`. Downstream G_prior carries the previous-iter
   values for the affected scales. Net effect: ebnm contributes
   nothing at the coarsest 1-2 scales (depending on `T_basis`), and
   those scales ride on the init from marginal data.

2. **Init-time lead-picker pathology.** If the marginal-data
   lead picker selects a noise variable on a degenerate fixture
   (no signal at scale `s`, all `Bhat[, idx_s]` are pure noise),
   ebnm fits a near-null prior at init. The IBSS-loop M-step
   then refits ebnm on the full `(|keep_idx| × |idx_s|)`
   rectangle (not the single lead variable), so a noisy init
   gets corrected as soon as `alpha` accumulates mass. Net
   effect: a bad init is recovered within a few IBSS iters.
   Unit test: assert the marginal-data lead picker selects a
   signal-bearing variable on the synthetic sparse fixture from
   the proposal's acceptance criteria (`set.seed(2L); n = 200;
   p = 30; T_m = 64; signal at variables c(7, 18) and times
   c(20, 44)`).

3. **ebnm version pinning.** The `point_normal` and
   `point_laplace` interfaces, plus the `g_init` / `fix_g`
   warm-start path and the `fitted_g` shape (`normalmix` for
   the Normal slab, `laplacemix` for the Laplace slab), are
   stable across ebnm 1.0.x and 1.1.x. Pin `ebnm (>= 1.0.0)`
   in `DESCRIPTION`. Tighten the pin only if a future ebnm
   release breaks the `g_init` contract observed here.
