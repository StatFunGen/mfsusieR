# S3 dispatch in mfsusieR

## What problem this solves

`susieR` implements the full IBSS algorithm — initialization, the outer
iteration loop, convergence checking, credible-set construction. Reimplementing
all of that in mfsusieR for multi-functional data would be error-prone and
hard to keep in sync. Instead, mfsusieR plugs into susieR's existing loop by
teaching it how to handle a new data class (`mf_individual`) through R's S3
dispatch mechanism.

---

## Primer: how S3 dispatch works in R

R's S3 system is simple. A **generic function** is a function that looks at the
class of its first argument and calls the right **method**:

```r
# Generic defined somewhere:
compute_residuals <- function(data, ...) UseMethod("compute_residuals")

# Method for class "individual" (susieR's scalar case):
compute_residuals.individual <- function(data, ...) { ... }

# Method for class "mf_individual" (mfsusieR's multi-functional case):
compute_residuals.mf_individual <- function(data, ...) { ... }
```

When susieR calls `compute_residuals(data, ...)`, R checks `class(data)`.
If it is `"mf_individual"`, it runs `compute_residuals.mf_individual`. If it
is `"individual"`, it runs `compute_residuals.individual`. The IBSS loop in
susieR does not need a single `if (multi_functional)` branch anywhere.

---

## The three objects that flow through the loop

Every generic in the system takes the same three arguments:

| Object | Class | Mutable? | What it holds |
|--------|-------|----------|---------------|
| `data` | `c("mf_individual", "individual")` | No | DWT-packed Y matrices, X, `na_idx`, `scale_index`, `T_basis`, etc. Set once by `create_mf_individual()`. |
| `params` | plain list | No | All user settings: `L`, `tol`, `prior`, `residual_variance_scope`, etc. Set once by `mfsusie()`. |
| `model` | `c("mfsusie", "susie")` | Yes | Fit state: `alpha`, `mu`, `mu2`, `sigma2`, `G_prior`, `elbo`, etc. Updated every iteration. |

The dispatch key is `class(data)[1]` = `"mf_individual"`. Every time susieR
calls a generic with `data`, R routes it to the `*.mf_individual` method.

---

## The call chain

```
mfsusie(X, Y, pos, ...)
  │
  ├─ create_mf_individual(...)      → data   (class: mf_individual)
  ├─ mf_prior_scale_mixture(...)    → prior
  ├─ assemble params list           → params
  │
  └─ susie_workhorse(data, params)
       │
       ├─ ibss_initialize(data, params)
       │    └─ initialize_susie_model.mf_individual(...)  ← dispatched
       │    └─ initialize_fitted.mf_individual(...)        ← dispatched
       │    └─ ibss_initialize.mf_individual(...)          ← dispatched
       │
       ├─ [outer loop: iter = 1 … max_iter]
       │    ├─ track_ibss_fit(data, …)                    ← dispatched
       │    ├─ ibss_fit(data, params, model)
       │    │    └─ [for l in 1:L]
       │    │         └─ single_effect_update(data, params, model, l)
       │    │              └─ single_effect_regression(data, params, model, l)
       │    │                   ├─ compute_ser_statistics(data,…,l)   ← dispatched
       │    │                   ├─ pre_loglik_prior_hook(data,…,l)    ← dispatched
       │    │                   ├─ loglik(data,…,l)                   ← dispatched
       │    │                   │    └─ optimize_prior_variance(…)    ← dispatched (inside loglik)
       │    │                   ├─ calculate_posterior_moments(…,l)   ← dispatched
       │    │                   ├─ compute_kl(…,l)                    ← dispatched
       │    │                   └─ post_loglik_prior_hook(data,…,l)   ← dispatched
       │    │              └─ update_fitted_values(data,…,l)          ← dispatched
       │    │
       │    ├─ get_objective(data, params, model)          ← dispatched (mfsusie method)
       │    └─ update_model_variance(data, params, model)  ← dispatched
       │         └─ update_variance_components(…)          ← dispatched
       │         └─ update_derived_quantities(…)           ← dispatched
       │
       ├─ trim_null_effects(data, params, model)           ← dispatched
       ├─ ibss_finalize(data, params, model, …)
       │    └─ get_scale_factors.mf_individual(…)         ← dispatched
       │    └─ get_intercept.mf_individual(…)             ← dispatched
       │    └─ get_fitted.mf_individual(…)                ← dispatched
       │    └─ get_variable_names.mf_individual(…)        ← dispatched
       │    └─ get_zscore.mf_individual(…)                ← dispatched
       │    └─ cleanup_extra_fields.mf_individual(…)      ← dispatched
       │    └─ Eloglik.mf_individual(…)                   ← dispatched
       │    └─ get_cs(…)                                  (susieR internal, called bare)
       │
       └─ return model  →  mfsusie() attaches class "mfsusie", dwt_meta, etc.
```

Every `← dispatched` call looks up `class(data)` at runtime and jumps to the
`*.mf_individual` implementation in `R/individual_data_methods.R` or
`R/ibss_methods.R`. susieR itself never branches on data type.

---

## Where the registration happens: `R/zzz.R`

Two problems must be solved:

1. **The generics live in susieR's namespace**, not mfsusieR's. R only finds
   `compute_residuals.mf_individual` if it is registered in the same namespace
   as the generic `compute_residuals`. A plain `@export` in NAMESPACE does not
   do this — it only exports from mfsusieR's namespace.

2. **susieR exports some internals** (e.g., `SER_posterior_e_loglik`,
   `update_variance_components`, `get_cs`) that mfsusieR's methods need to
   call. They cannot be accessed by `susieR::fn` because they are not
   re-exported. They must be fetched directly from `asNamespace("susieR")` at
   load time.

Both are solved in `.onLoad` in `R/zzz.R`:

```r
.onLoad <- function(libname, pkgname) {
  susie_ns <- asNamespace("susieR")
  pkg_ns   <- asNamespace(pkgname)

  # 1. Cache susieR internals as package-level bindings.
  #    Methods can then call them bare (e.g., get_cs(...)) without
  #    knowing they come from susieR.
  for (fn in c("SER_posterior_e_loglik", "update_variance_components",
               "get_cs", "initialize_susie_model", ...)) {
    assign(fn, get(fn, envir = susie_ns), envir = pkg_ns)
  }

  # 2. Register S3 methods into susieR's namespace so that when
  #    susieR calls compute_residuals(data, ...) with data of class
  #    "mf_individual", R dispatches to our method.
  for (g in mf_generics) {
    method_fn <- get(paste0(g, ".mf_individual"), envir = pkg_ns)
    registerS3method(g, "mf_individual", method_fn, envir = susie_ns)
  }

  # 3. Same for the "mfsusie" model class.
  for (g in mfsusie_generics) {
    method_fn <- get(paste0(g, ".mfsusie"), envir = pkg_ns)
    registerS3method(g, "mfsusie", method_fn, envir = susie_ns)
  }
}
```

`registerS3method(generic, class, method, envir = susie_ns)` is the key call.
It inserts the method directly into susieR's S3 method table, so R's dispatch
finds it even though the method lives in mfsusieR.

---

## Two dispatch classes

### Class `mf_individual` (the data object)

Governs every computation that touches the data — residuals, Bhat/Shat,
likelihoods, variance updates, fitted values. All these methods live in:

- `R/individual_data_methods.R` — per-effect SER-step methods
- `R/ibss_methods.R` — per-iteration methods + init/finalize accessors

The full list registered under `mf_individual`:

| Phase | Generic | What the method does |
|-------|---------|----------------------|
| Init | `get_var_y` | Initial σ² estimate from wavelet Y variance |
| Init | `initialize_susie_model` | Builds model struct with `G_prior`, `sigma2`, `fitted_g_per_effect`, etc. |
| Init | `initialize_fitted` | Sets initial fitted values (zeros for wavelet fits) |
| Init | `ibss_initialize` | Attaches iter-0 state: `pi_V`, caches xtx etc. |
| Per-effect | `compute_residuals` | $X^T R_m$ over `na_idx[[m]]` rows only |
| Per-effect | `compute_ser_statistics` | Computes (Bhat, Shat) per outcome, scale group |
| Per-effect | `loglik` | Calls `optimize_prior_variance` inside; returns log-BF vector |
| Per-effect | `optimize_prior_variance` | Runs mixsqp / ebnm per (outcome, scale) group |
| Per-effect | `calculate_posterior_moments` | Posterior mean/var of effect for each outcome |
| Per-effect | `compute_kl` | KL(q ‖ prior) for effect l |
| Per-effect | `neg_loglik` | Negative log-likelihood (used by outer optimizer) |
| Per-effect | `update_fitted_values` | Updates residual fitted values Xr after effect l |
| Per-effect | `SER_posterior_e_loglik` | E[log p(Y \| β_l)] for ELBO computation |
| Per-iter | `update_variance_components` | Updates σ² per modality (per_outcome or per_scale) |
| Per-iter | `update_derived_quantities` | Rebuilds shat2 cache after σ² update |
| Per-iter | `update_model_variance` | Wrapper that calls the two above |
| Per-iter | `track_ibss_fit` | Records per-iter snapshots when `track_fit=TRUE` |
| Per-iter | `trim_null_effects` | Removes effects with max α below threshold |
| Finalize | `get_scale_factors`, `get_intercept`, `get_fitted`, `get_variable_names`, `get_zscore` | Extract components for the returned fit object |
| Finalize | `cleanup_extra_fields` | Drops internal caches (`D`, `xtx_diag_list`, etc.) |
| Finalize | `Eloglik` | Final ELBO log-likelihood term |

### Class `mfsusie` (the model object)

Governs computations that read from the model struct and need to understand
its multi-outcome, multi-scale structure. Registered under `mfsusie`:

| Generic | What the method does |
|---------|----------------------|
| `get_objective` | Refreshes per-effect lbf/KL before returning the ELBO. Without this, the per-iteration ELBO would be a hybrid quantity (Eloglik at iter-final state but KL at the state when each effect was last updated). |
| `format_sigma2_summary` | Formats σ² as a compact string for the verbose progress line. Needed because σ² is a list of vectors (one per modality) rather than a scalar. |
| `format_extra_diag` | Appends null-mass summary to the verbose line. |
| `get_posterior_mean_l` | Returns the per-outcome posterior mean curve for effect l. |
| `get_posterior_mean_sum` | Sums across effects. |
| `get_posterior_moments_l` | Returns (mu, mu2) for effect l. |

---

## Concrete example: `compute_residuals`

When the IBSS loop updates effect l=3, susieR calls:

```r
compute_residuals(data, params, model, l = 3)
```

1. R evaluates `class(data)` → `c("mf_individual", "individual")`
2. R looks up `compute_residuals.mf_individual` in susieR's method table
   (registered there by `.onLoad`).
3. It finds and calls:

```r
# R/individual_data_methods.R:26
compute_residuals.mf_individual <- function(data, params, model, l, ...) {
  R_list <- lapply(seq_len(data$M), function(m) {
    idx_m <- data$na_idx[[m]]   # observed rows for modality m
    Xr_m  <- model$Xr[[m]]     # n x T_basis[m] fitted values
    Xr_l  <- get_Xr_l(data, model, m, l)  # contribution of effect l
    R_m   <- data$D[[m]][idx_m, ] - Xr_m[idx_m, ] + Xr_l[idx_m, ]
    crossprod(data$X_processed[idx_m, ], R_m)  # p x T_basis[m]
  })
  model$XtR <- R_list
  model
}
```

The method computes X^T R separately per modality, using only the `na_idx[[m]]`
rows for each. susieR's generic never sees any of this — it just calls
`compute_residuals(data, ...)` and gets back an updated `model`.

---

## Concrete example: `update_variance_components`

At the end of each IBSS iteration, susieR calls:

```r
update_model_variance(data, params, model)
  └─ update_variance_components(data, params, model)   # dispatched
  └─ update_derived_quantities(data, params, model)    # dispatched
```

`update_variance_components.mf_individual` (`R/individual_data_methods.R:436`)
reads `params$residual_variance_scope` and updates `model$sigma2[[m]]` as
either a scalar (`per_outcome`) or a length-S_m vector (`per_scale`):

```r
update_variance_components.mf_individual <- function(data, params, model, ...) {
  method <- params$residual_variance_scope %||% "per_outcome"
  for (m in seq_len(data$M)) {
    er2_t <- mf_get_ER2_per_position(data, model, m)
    n     <- length(data$na_idx[[m]])
    if (method == "per_outcome") {
      model$sigma2[[m]] <- sum(er2_t) / (n * data$T_basis[m])
    } else {
      indx_m <- data$scale_index[[m]]
      model$sigma2[[m]] <- vapply(indx_m, function(idx)
        sum(er2_t[idx]) / (n * length(idx)), numeric(1))
    }
  }
  refresh_iter_cache.mf_individual(data, model)
}
```

susieR's `update_model_variance` generic calls this transparently, without
knowing anything about wavelet scales or modalities.

---

## Why `get_objective` needs a `mfsusie` method

susieR's default `get_objective` reads `model$elbo` at the current state.
In mfsusieR the per-effect KL and lbf are updated when effect l is visited
(inside `single_effect_regression`), but by the time all L effects have been
updated in one iteration, the lbf[l] and KL[l] values for early effects (l=1,
l=2, …) were computed against an older α state. The final ELBO should reflect
the end-of-iteration state for all effects uniformly.

`get_objective.mfsusie` calls `refresh_lbf_kl.mf_individual` first, which
recomputes lbf[l] and KL[l] for all l against the iter-final α and π_V before
summing the ELBO. Without this, the per-iteration ELBO is a hybrid quantity
that is neither the start-of-iteration nor the end-of-iteration value.

---

## File map

| File | What it contains |
|------|-----------------|
| `R/zzz.R` | `.onLoad`: caches susieR internals, registers all S3 methods via `registerS3method` |
| `R/individual_data_class.R` | `create_mf_individual()`: constructor that gives `data` its `mf_individual` class |
| `R/individual_data_methods.R` | Per-effect SER-step methods (`*.mf_individual`) + per-iter variance methods |
| `R/ibss_methods.R` | Per-iter + init/finalize methods (`*.mf_individual`), `get_objective.mfsusie` |
| `R/model_class.R` | `initialize_susie_model.mf_individual`: builds the model struct with `G_prior`, `sigma2`, `fitted_g_per_effect` |
| `R/mfsusie_methods.R` | User-facing S3: `predict.mfsusie`, `coef.mfsusie`, `summary.mfsusie`, `print.mfsusie` — these are NOT registered into susieR, just exported normally |
| `R/prior_scale_mixture.R` | `optimize_prior_variance.mf_individual`: mixsqp / ebnm EM per (outcome, scale) |

---

## Why this design

**susieR handles everything algorithmic.** The outer loop, convergence, ELBO
bookkeeping, credible-set construction, PIP aggregation — none of that is
reimplemented. mfsusieR only implements the data-specific computations.

**Adding a new data type requires implementing ~20 methods.** Each is a
self-contained function with a clear contract (inputs / outputs defined by the
generic). No global state, no hidden branches.

**The dispatch boundary makes bugs locatable.** If sigma2 is wrong, the bug
is in `update_variance_components.mf_individual`. If BFs are wrong, it is in
`loglik.mf_individual` or `compute_ser_statistics.mf_individual`. The generic
name tells you the lifecycle phase; the class suffix tells you which file.
