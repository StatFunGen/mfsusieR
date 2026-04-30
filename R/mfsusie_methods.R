# User-facing S3 methods on `mfsusie` fit objects.
#
# - predict / coef / fitted: project the posterior coefficients
#   back to the original Y scale via the per-outcome inverse
#   wavelet transform stored on `fit$dwt_meta`.
# - print / summary: terse and detailed reports.
#
# All methods read `fit$dwt_meta` (a stash of the wavelet pipeline
# parameters set by `mfsusie()`) so the user does not need to keep
# the original `mf_individual` data class around.

# ---- helpers --------------------------------------------------

# Per-outcome inverse-wavelet projection of one p x T_basis
# coefficient matrix back to a p x T_m matrix on the original
# position grid. `coef_wavelet` rows are variables (or some other
# index); each row is inverted independently.
mf_invert_per_outcome <- function(coef_wavelet, m, dwt_meta,
                                  intercept = TRUE) {
  T_pad <- dwt_meta$T_basis[m]
  # Match susieR's `coef.susie` convention where the intercept is a
  # separable component: response reconstructions (`predict`,
  # `fitted`) keep it; per-effect (slope) reconstructions (`coef`)
  # set `intercept = FALSE` so the Y baseline `cm_Y` is not added
  # to an effect curve. `csd_Y` is always applied because mu lives
  # in Y-standardized units.
  cm    <- if (intercept) dwt_meta$column_center[[m]] else rep(0, T_pad)
  csd   <- dwt_meta$column_scale[[m]]
  wcm   <- if (intercept) dwt_meta$wavelet_center[[m]] else rep(0, T_pad)
  wcsd  <- dwt_meta$wavelet_scale[[m]]

  # `mf_invert_dwt` returns one curve per row of `coef_wavelet`,
  # on the post-remap padded grid stored at `dwt_meta$pos[[m]]`
  # (length `T_pad`). Callers (predict / coef / fitted) work on
  # that grid throughout, so no per-call trimming is needed.
  mf_invert_dwt(
    D_packed       = coef_wavelet,
    column_center  = cm,
    column_scale   = csd,
    wavelet_center = wcm,
    wavelet_scale  = wcsd,
    filter_number  = dwt_meta$wavelet_filter,
    family         = dwt_meta$wavelet_family
  )
}

# Per-effect coefficient matrix (p x T_basis) for effect l, outcome m,
# on the standardized X scale: alpha_lj * mu_lj_t.
mf_effect_wavelet <- function(fit, l, m) {
  fit$alpha[l, ] * fit$mu[[l]][[m]]
}

# Coefficient sum across all L effects per outcome, on standardized
# X scale, in the wavelet domain.
mf_total_wavelet <- function(fit, m) {
  p     <- ncol(fit$alpha)
  T_pad <- fit$dwt_meta$T_basis[m]
  out   <- matrix(0, nrow = p, ncol = T_pad)
  for (l in seq_len(nrow(fit$alpha))) {
    out <- out + mf_effect_wavelet(fit, l, m)
  }
  out
}

# ---- predict --------------------------------------------------

#' Predicted response curves for new covariates
#'
#' Projects the posterior coefficient curves through new covariate
#' values to produce predicted response curves on the original Y
#' scale, per outcome. The wavelet pipeline that built the fit
#' (column scaling, padding, DWT) is inverted to yield curves on
#' each outcome's original position grid `pos[[m]]`.
#'
#' @param object an `mfsusie` fit returned by `mfsusie()`.
#' @param newx numeric matrix `n_new x p` of new covariates on the
#'   same scale as the training `X`. `NULL` (default) returns the
#'   training fitted values (equivalent to `fitted(object)`).
#' @param ... ignored.
#' @return list of length `M`; each element a numeric matrix
#'   `n_new x T_m` of predicted curves on the original position
#'   grid for that outcome.
#' @details Prediction uses the per-variable alpha-weighted aggregate
#'   coefficient `b[j, t] = sum_l alpha[l, j] * mu[l, j, t] /
#'   csd_X[j]`, the same form as `susieR::predict.susie()`.
#' @export
predict.mfsusie <- function(object, newx = NULL, ...) {
  if (is.null(newx)) return(fitted.mfsusie(object))
  if (!is.matrix(newx)) stop("`newx` must be a numeric matrix.")
  if (ncol(newx) != ncol(object$alpha)) {
    stop(sprintf(
      "`newx` has %d columns but the fit was trained on X with %d columns.",
      ncol(newx), ncol(object$alpha)))
  }

  meta <- object$dwt_meta
  M    <- meta$M
  out  <- vector("list", M)

  newx_std <- sweep(newx, 2, meta$X_center, "-")
  if (any(meta$X_scale != 1)) {
    newx_std <- sweep(newx_std, 2, meta$X_scale, "/")
  }
  for (m in seq_len(M)) {
    coef_wavelet <- mf_total_wavelet(object, m)        # p x T_basis
    pred_wavelet <- newx_std %*% coef_wavelet           # n_new x T_basis
    out[[m]] <- mf_invert_per_outcome(pred_wavelet, m, meta)
  }
  out
}

# ---- coef --------------------------------------------------

#' Per-effect coefficient curves on the original X scale
#'
#' Returns the per-effect, per-outcome coefficient curves on the
#' original (unstandardized) X scale, projected back through the
#' inverse wavelet transform to the original position grid.
#'
#' @param object an `mfsusie` fit.
#' @param smooth_method optional character. When supplied, returns
#'   the post-smoothed effect curves stored at
#'   `object$smoothed[[smooth_method]]` instead of the raw
#'   wavelet-inverse curves. The named smoother must have been
#'   applied via `mf_post_smooth(object, method = smooth_method)`
#'   beforehand. When `NULL` (default), `coef.mfsusie` returns
#'   the raw curves; smoothing is an explicit opt-in.
#' @param ... ignored.
#' @return list of length `M`; each element an `L x T_m` matrix
#'   whose row `l` is effect `l`'s coefficient curve for outcome
#'   `m`. When `smooth_method` is supplied, the returned curves
#'   are on the same shape but smoothed by the named method;
#'   the returned object carries an attribute
#'   `attr(., "smooth_method")` recording the choice.
#' @export
coef.mfsusie <- function(object, smooth_method = NULL, ...) {
  if (!is.null(smooth_method)) {
    if (is.null(object$smoothed) ||
        is.null(object$smoothed[[smooth_method]])) {
      avail <- if (is.null(object$smoothed)) character(0L)
               else names(object$smoothed)
      stop(sprintf(
        "Post-smoothed curves for method = '%s' are not on the fit. ",
        smooth_method),
        if (length(avail) > 0L)
          sprintf("Applied methods: %s. ",
                  paste(sQuote(avail), collapse = ", "))
        else "No methods have been applied. ",
        "Run `mf_post_smooth(fit, method = '", smooth_method, "')` first.")
    }
    sm  <- object$smoothed[[smooth_method]]$effect_curves
    # `effect_curves` is `list[M]` of `list[L]` of length-T_m
    # numeric vectors. Convert to the same shape `coef.mfsusie`
    # returns for the raw path: `list[M]` of `L x T_m` matrices.
    out <- lapply(sm, function(per_outcome) {
      do.call(rbind, per_outcome)
    })
    attr(out, "smooth_method") <- smooth_method
    return(out)
  }

  meta <- object$dwt_meta
  L    <- nrow(object$alpha)
  M    <- meta$M
  X_scale <- meta$X_scale
  out  <- vector("list", M)
  for (m in seq_len(M)) {
    coef_l_wavelet <- matrix(0, nrow = L, ncol = meta$T_basis[m])
    for (l in seq_len(L)) {
      # Per-variable unscale of mu: mu lives in (X-standardized,
      # Y-standardized-then-DWT'd) units; dividing by `X_scale[j]`
      # converts it to (X-raw, Y-standardized-then-DWT'd) units so
      # the alpha-weighted sum is the per-effect coefficient on the
      # raw X scale.
      mu_raw_X <- sweep(object$mu[[l]][[m]], 1L, X_scale, "/")
      coef_l_wavelet[l, ] <- colSums(object$alpha[l, ] * mu_raw_X)
    }
    # Inverse DWT with `column_center = 0` so the Y intercept
    # `cm_Y` is NOT added back to a per-effect estimate. The `csd_Y`
    # scale is still applied (mu is in Y-standardized units).
    out[[m]] <- mf_invert_per_outcome(coef_l_wavelet, m, meta,
                                      intercept = FALSE)
  }
  out
}

# ---- fitted --------------------------------------------------

#' Fitted response curves on the training X
#'
#' Projects the running per-outcome wavelet-domain fit
#' `fit$fitted[[m]]` (left on the fit by the IBSS loop and reused
#' by warm-starts) back to the original position grid via the
#' inverse wavelet transform. No `X` needed: the running fit is
#' the per-individual training prediction in wavelet form.
#'
#' @param object an `mfsusie` fit.
#' @param ... ignored.
#' @return list of length `M`; each element a numeric matrix
#'   `n x T_m` of fitted curves on the training X.
#' @export
fitted.mfsusie <- function(object, ...) {
  meta <- object$dwt_meta
  M    <- meta$M
  out  <- vector("list", M)
  for (m in seq_len(M)) {
    out[[m]] <- mf_invert_per_outcome(object$fitted[[m]], m, meta)
  }
  out
}

# ---- print --------------------------------------------------

#' Print method for `mfsusie` fits
#'
#' Compact, user-facing one-screen summary: dimensions, convergence,
#' top PIPs, credible-set membership counts.
#'
#' @param x an `mfsusie` fit.
#' @param ... ignored.
#' @return `invisible(x)`.
#' @export
print.mfsusie <- function(x, ...) {
  meta <- x$dwt_meta
  cat("mfsusie fit\n")
  cat(sprintf("  p (predictors): %d\n", ncol(x$alpha)))
  cat(sprintf("  L (effects):    %d\n", nrow(x$alpha)))
  cat(sprintf("  M (outcomes): %d\n", meta$M))
  T_str <- paste(meta$T_basis, collapse = ", ")
  cat(sprintf("  T_basis:       (%s)\n", T_str))
  cat(sprintf("  iterations:     %d %s\n", x$niter %||% 0L,
              if (isTRUE(x$converged)) "(converged)" else "(NOT converged)"))
  if (!is.null(x$elbo)) {
    cat(sprintf("  ELBO (last):    %.4f\n", x$elbo[length(x$elbo)]))
  }
  if (!is.null(x$sets$cs)) {
    n_cs <- length(x$sets$cs)
    cs_sizes <- if (n_cs > 0L) sapply(x$sets$cs, length) else integer(0)
    cat(sprintf("  credible sets:  %d (sizes: %s)\n",
                n_cs, paste(cs_sizes, collapse = ", ")))
  }
  if (!is.null(x$pip)) {
    top <- order(x$pip, decreasing = TRUE)[seq_len(min(5L, length(x$pip)))]
    cat("  top PIPs:\n")
    for (j in top) {
      jname <- if (!is.null(names(x$pip))) names(x$pip)[j] else as.character(j)
      cat(sprintf("    %-12s %.4f\n", jname, x$pip[j]))
    }
  }
  invisible(x)
}

# ---- summary --------------------------------------------------

#' Summary method for `mfsusie` fits
#'
#' Returns a list with the aggregate fit metadata and per-CS
#' summaries; can be printed for human inspection.
#'
#' @param object an `mfsusie` fit.
#' @param ... ignored.
#' @return an object of class `summary.mfsusie` (a list).
#' @export
summary.mfsusie <- function(object, ...) {
  meta <- object$dwt_meta
  alpha <- object$alpha
  pip   <- object$pip
  sets  <- object$sets

  cs_table <- if (!is.null(sets$cs) && length(sets$cs) > 0L) {
    do.call(rbind, lapply(seq_along(sets$cs), function(i) {
      cs   <- sets$cs[[i]]
      data.frame(
        cs_index    = i,
        size        = length(cs),
        variables   = paste(cs, collapse = ","),
        purity      = if (!is.null(sets$purity)) sets$purity[i, "min.abs.corr"] else NA_real_,
        coverage    = if (!is.null(sets$coverage)) sets$coverage[i] else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    NULL
  }

  # Summarize null mass across all (l, m, s) cells. pi_V is per-effect:
  # list[L] of list[M] of S_m x K matrix; column 1 is the null component.
  pi_null_summary <- if (!is.null(object$pi_V)) {
    null_mass <- unlist(lapply(object$pi_V, function(pi_l)
      lapply(pi_l, function(piVm) piVm[, 1L])))
    if (length(null_mass) > 0L) {
      list(min    = min(null_mass),
           median = stats::median(null_mass),
           max    = max(null_mass))
    } else NULL
  } else NULL

  out <- list(
    n_effects   = nrow(alpha),
    n_variables      = ncol(alpha),
    n_outcomes = meta$M,
    T_basis    = meta$T_basis,
    converged   = isTRUE(object$converged),
    n_iter      = object$niter %||% 0L,
    elbo_final  = if (!is.null(object$elbo)) object$elbo[length(object$elbo)] else NA_real_,
    pip         = pip,
    pi_null     = pi_null_summary,
    cs          = cs_table
  )
  class(out) <- "summary.mfsusie"
  out
}

#' @export
print.summary.mfsusie <- function(x, ...) {
  cat(sprintf("mfsusie summary: p=%d, L=%d, M=%d, %s in %d iter\n",
              x$n_variables, x$n_effects, x$n_outcomes,
              if (x$converged) "converged" else "NOT converged",
              x$n_iter))
  cat(sprintf("  T_basis per outcome: (%s)\n",
              paste(x$T_basis, collapse = ", ")))
  cat(sprintf("  Final ELBO: %.4f\n", x$elbo_final))
  if (!is.null(x$pi_null)) {
    cat(sprintf("  Mixture null-mass across (m, s): min=%.3f, median=%.3f, max=%.3f\n",
                x$pi_null$min, x$pi_null$median, x$pi_null$max))
  }
  if (!is.null(x$cs) && nrow(x$cs) > 0L) {
    cat("  Credible sets:\n")
    for (i in seq_len(nrow(x$cs))) {
      cat(sprintf("    CS %d: size=%d, purity=%.3f, variables=%s\n",
                  x$cs$cs_index[i], x$cs$size[i], x$cs$purity[i],
                  x$cs$variables[i]))
    }
  } else {
    cat("  No credible sets.\n")
  }
  invisible(x)
}

# ---- summarize_effects (affected-region accessor) ------------

#' Per-(CS, outcome) regions where the credible band excludes zero
#'
#' Given an `mfsusie` fit that has been post-smoothed via
#' `mf_post_smooth()`, returns one row per `(CS, outcome,
#' contiguous-run)` triple identifying the position ranges
#' where the smoothed credible band excludes zero. Useful as a
#' compact diagnostic for "where on the curve does this CS act
#' on this outcome?".
#'
#' @param fit an `mfsusie` / `fsusie()` fit. The fit must
#'   carry post-smoothed credible bands; call
#'   `mf_post_smooth(fit, method = ...)` first.
#' @param smooth_method optional name of the smoother to read
#'   from `fit$smoothed`. When `NULL` (default) the smoother
#'   priority order is used: `"TI" > "smash" > "HMM" >
#'   "scalewise"`.
#' @return A data frame with columns
#'   \describe{
#'     \item{`cs_index`}{integer, the credible-set index
#'       (1..length(fit$sets$cs_index)).}
#'     \item{`outcome`}{integer, the outcome index `m`.}
#'     \item{`start`, `end`}{integer position indices
#'       defining the contiguous run on the
#'       `T_basis[m]`-position grid where the band excludes
#'       zero.}
#'     \item{`n_positions`}{integer, run length
#'       (`end - start + 1`).}
#'   }
#'   When no CS has any band-excludes-zero run, returns an
#'   empty data frame with the same columns.
#' @examples
#' \donttest{
#' set.seed(1L)
#' n <- 100; p <- 20; T_m <- 32L
#' X <- matrix(rnorm(n * p), n)
#' beta <- numeric(p); beta[3] <- 1.2
#' shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
#' Y <- X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
#'        matrix(rnorm(n * T_m, sd = 0.3), n)
#' fit <- fsusie(Y, X, L = 1, verbose = FALSE)
#' fit_s <- mf_post_smooth(fit, method = "TI")
#' mf_summarize_effects(fit_s)
#' }
#' @export
mf_summarize_effects <- function(fit, smooth_method = NULL) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  picked <- .pick_smooth_method(fit, smooth_method)
  if (is.null(picked)) {
    stop("`fit` has no post-smoothed credible bands. Run ",
         "`mf_post_smooth(fit, method = ...)` first.")
  }
  smoothed <- fit$smoothed[[picked]]
  bands    <- smoothed$credible_bands
  cs_idx   <- fit$sets$cs_index %||% integer(0L)
  M        <- length(bands)

  empty <- data.frame(
    cs_index    = integer(0L),
    outcome     = integer(0L),
    start       = integer(0L),
    end         = integer(0L),
    n_positions = integer(0L)
  )
  if (length(cs_idx) == 0L) return(empty)

  rows <- list()
  for (m in seq_len(M)) {
    bands_m <- bands[[m]]
    for (i in seq_along(cs_idx)) {
      l <- cs_idx[i]
      band <- bands_m[[l]]
      if (is.null(band)) next
      runs <- credibly_nonzero_runs(band)
      for (run in runs) {
        rows[[length(rows) + 1L]] <- data.frame(
          cs_index    = as.integer(i),
          outcome     = as.integer(m),
          start       = as.integer(run[1L]),
          end         = as.integer(run[2L]),
          n_positions = as.integer(run[2L] - run[1L] + 1L)
        )
      }
    }
  }
  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}

# Post-processing of effect curves on an `mfsusie` fit.
#
# `mf_post_smooth(fit)` returns the fit with three new slots:
#   $effect_curves[[m]][[l]]  : length-T_basis[m] smoothed curve
#   $credible_bands[[m]][[l]] : T_basis[m] x 2 [lower, upper]
#   $lfsr_curves[[m]][[l]]    : length-T_basis[m] in [0, 1]
#
# `method = "HMM"` derives lfsr from its mixture posterior. All
# other methods derive lfsr from the posterior (mean, sd) under
# the Gaussian approximation:
#   lfsr(t) = pnorm(-|mean(t)| / sd(t)).
#
# All methods operate on the fit alone when it carries the
# smoothing inputs (`fit$Y_grid` and `fit$X_eff`, attached by
# default at fit time). When the fit was built with
# `attach_smoothing_inputs = FALSE`, `mf_post_smooth(fit, X, Y)`
# attaches them on the fly from the supplied data.

# Internal: per-position local false sign rate from a Gaussian
# posterior summary. Returns `pnorm(-|mean| / sd)`, with sd
# floored at machine epsilon to avoid divide-by-zero.
lfsr_from_gaussian <- function(mean, sd) {
  sd <- pmax(sd, .Machine$double.eps)
  stats::pnorm(-abs(mean) / sd)
}

# Internal: rebuild and attach `Y_grid` and `X_eff` to a fit that
# was created with `attach_smoothing_inputs = FALSE`. Mirrors the
# attachment block at the end of `mfsusie()`.
.attach_smoothing_inputs <- function(fit, X, Y) {
  meta <- fit$dwt_meta
  if (!is.matrix(X) || ncol(X) != ncol(fit$alpha)) {
    stop("`X` must be the original n x p training genotype matrix.")
  }
  if (!is.list(Y) || length(Y) != meta$M) {
    stop(sprintf("`Y` must be a length-%d list of response matrices.",
                 meta$M))
  }
  fit$Y_grid <- vector("list", meta$M)
  for (m in seq_len(meta$M)) {
    Y_m <- as.matrix(Y[[m]])
    if (ncol(Y_m) != meta$T_basis[m]) {
      Y_m <- remap_data(Y_m, meta$pos[[m]], verbose = FALSE,
                        max_scale = log2(meta$T_basis[m]))$Y
    }
    fit$Y_grid[[m]] <- Y_m
  }
  L <- nrow(fit$alpha)
  fit$X_eff <- vector("list", L)
  X_centered <- if (any(meta$X_center != 0))
    sweep(X, 2, meta$X_center, "-") else X
  for (l in seq_len(L)) {
    # alpha-weighted aggregate of (raw) X columns
    fit$X_eff[[l]] <- as.numeric(X_centered %*% fit$alpha[l, ])
  }
  fit
}

#' Post-smooth a fit's per-effect curves and add credible bands
#'
#' Four smoothing methods are dispatched by `method`:
#'
#' - `"scalewise"` -- per-scale soft-thresholding of the lead
#'   variable's wavelet posterior mean. Fast, no iterations,
#'   uses only the wavelet posterior moments. Suitable for quick
#'   visual cleanup.
#' - `"TI"` (default) -- cycle-spinning translation-invariant
#'   denoising (Coifman & Donoho 1995). For each effect, isolates
#'   the per-effect residual response (in position space),
#'   regresses onto the alpha-weighted X aggregate, applies the
#'   stationary wavelet transform row-by-row, scalewise
#'   `ashr::ash` shrinkage on wavelet coefficients, and inverts
#'   via cycle-spinning average. Produces tighter credible bands
#'   than scalewise.
#' - `"HMM"` -- hidden Markov denoising on per-position regression
#'   coefficients. Yields a posterior mean curve plus a mixture-
#'   posterior local false sign rate.
#' - `"smash"` -- `smashr::smash.gaus` empirical-Bayes wavelet
#'   shrinkage on the per-position regression estimate. Requires
#'   the `smashr` package (Suggests).
#'
#' All four methods populate `lfsr_curves`. For TI, scalewise,
#' and smash the lfsr is `pnorm(-|mean| / sd)` under the
#' Gaussian-posterior approximation; for HMM it comes from the
#' mixture posterior directly.
#'
#' By default the smoother reads `fit$Y_grid` (post-remap Y) and
#' `fit$X_eff` (per-effect alpha-weighted X aggregate), both
#' attached by `mfsusie()` when `attach_smoothing_inputs = TRUE`.
#' When the fit was built with `attach_smoothing_inputs = FALSE`
#' you can pass the original `X` and `Y` here and the smoother
#' will compute and attach them on the fly.
#'
#' @param fit a fit returned by `mfsusie()` or `fsusie()`.
#' @param method one of `"scalewise"`, `"TI"`, `"HMM"`.
#' @param level numeric in (0, 1), credible-band coverage.
#'   Default 0.95.
#' @param threshold_factor numeric, multiplier on the universal
#'   `sqrt(2 log T)` threshold for `method = "scalewise"`. Ignored
#'   by other methods.
#' @param wavelet_filter integer, wavelet filter number for the TI
#'   stationary-wavelet transform. Default 1 (Haar).
#' @param wavelet_family character, wavelet family for the TI
#'   stationary-wavelet transform. Default `"DaubExPhase"`.
#' @param halfK integer, half-grid size for the HMM `fit_hmm`
#'   helper. Default 20.
#' @param X optional numeric matrix `n x p`, the original genotype
#'   matrix. Required when the fit was built with
#'   `attach_smoothing_inputs = FALSE`; ignored otherwise.
#' @param Y optional list of length `M` of numeric matrices, the
#'   original per-outcome response matrices. Required when the
#'   fit was built with `attach_smoothing_inputs = FALSE`; ignored
#'   otherwise.
#' @param ... extra arguments forwarded to the underlying
#'   shrinkage tool for the chosen method. Routing:
#'   `method = "ash"` and `method = "TI"` and `method = "HMM"`
#'   forward to `ashr::ash` (e.g., `nullweight = 10` for milder
#'   null shrinkage). `method = "smash"` forwards to
#'   `smashr::smash.gaus`. `method = "TI"` also accepts
#'   `scaling = c("per_scale", "uniform")` (default `"per_scale"`)
#'   to pick between scalewise per-scale ash and a single uniform
#'   ash on raw Y. `method = "scalewise"` ignores `...`.
#' @return the input fit with `$effect_curves`, `$credible_bands`,
#'   and `$lfsr_curves` populated. Scalar outcomes
#'   (`T_basis[m] = 1`) skip the wavelet step (smoothing is a
#'   no-op there).
#' @export
mf_post_smooth <- function(fit,
                           method           = c("TI", "scalewise",
                                                "HMM", "smash", "ash"),
                           level            = 0.95,
                           threshold_factor = 1,
                           wavelet_filter   = 1L,
                           wavelet_family   = "DaubExPhase",
                           halfK            = 20L,
                           overwrite_previous = FALSE,
                           X = NULL, Y = NULL,
                           ...) {
  if (!inherits(fit, "mfsusie")) {
    stop("`fit` must be an `mfsusie` (or `fsusie`) fit object.")
  }
  method <- match.arg(method)
  if (level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).")
  }

  # `scalewise` operates on the wavelet posterior alone and needs
  # neither Y_grid nor X_eff. The other smoothers need both. If
  # they're missing on the fit, attach them from the supplied
  # `X` / `Y`.
  if (method != "scalewise" &&
      (is.null(fit$Y_grid) || is.null(fit$X_eff))) {
    if (is.null(X) || is.null(Y)) {
      stop(sprintf(
        paste0("method = '%s' needs the smoothing inputs ",
               "(`fit$Y_grid`, `fit$X_eff`). The fit was built with ",
               "`attach_smoothing_inputs = FALSE`; pass `X` and `Y` ",
               "explicitly, or refit with the default."),
        method))
    }
    fit <- .attach_smoothing_inputs(fit, X = X, Y = Y)
  }
  if (method == "smash" && !requireNamespace("smashr", quietly = TRUE)) {
    stop("method = 'smash' requires the `smashr` package. ",
         "Install via remotes::install_github(\"stephenslab/smashr\") ",
         "or pass `method = \"ash\"` for the smashr-free alternative.")
  }

  if (!isTRUE(overwrite_previous) &&
      !is.null(fit$smoothed) && !is.null(fit$smoothed[[method]])) {
    stop(sprintf(
      "Smoothing for method = '%s' is already on this fit. ",
      method),
      "Pass `overwrite_previous = TRUE` to recompute and replace it.")
  }

  payload <- switch(method,
    "scalewise" = .post_smooth_scalewise(fit, level, threshold_factor),
    "TI"        = .post_smooth_ti(fit, level, wavelet_filter,
                                  wavelet_family, ...),
    "HMM"       = .post_smooth_hmm(fit, level, halfK, ...),
    "smash"     = .post_smooth_smash(fit, level, flavor = "smashr", ...),
    "ash"       = .post_smooth_smash(fit, level, flavor = "ash", ...))

  if (is.null(fit$smoothed)) fit$smoothed <- list()
  fit$smoothed[[method]] <- payload
  fit
}

# Priority order used by `mfsusie_plot()` and friends when the
# fit carries multiple smoothings and the user has not named
# one explicitly.
.smoother_priority <- c("TI", "smash", "ash", "HMM", "scalewise")

# Pick the highest-priority smoothing on `fit$smoothed` and emit
# a hint listing the others. Returns NULL when no smoothing is
# present.
.pick_smooth_method <- function(fit, requested = NULL) {
  if (!is.null(requested)) {
    if (is.null(fit$smoothed) || is.null(fit$smoothed[[requested]])) {
      avail <- if (is.null(fit$smoothed)) character(0L)
               else names(fit$smoothed)
      stop(sprintf(
        "Smoothing for method = '%s' is not on this fit. ", requested),
        if (length(avail) > 0L)
          sprintf("Applied methods: %s.",
                  paste(sQuote(avail), collapse = ", "))
        else "No smoothings have been applied.")
    }
    return(requested)
  }
  if (is.null(fit$smoothed) || length(fit$smoothed) == 0L) return(NULL)
  # `mf_post_smooth(method = ...)` constrains method names via
  # `match.arg` to the known set, so every name in `fit$smoothed`
  # is in `.smoother_priority` and the intersection is non-empty.
  applied <- intersect(.smoother_priority, names(fit$smoothed))
  picked  <- applied[1L]
  if (length(applied) > 1L) {
    others <- setdiff(applied, picked)
    msg <- sprintf(
      "Plotting smoothing '%s'. Other smoothings on this fit: %s. Pass `smooth_method = '<name>'` to plot a different one.",
      picked, paste(sQuote(others), collapse = ", "))
    warning_message(msg, style = "hint")
  }
  picked
}

# Higher-order loop shared by every smoother. Initializes the four
# payload lists, dispatches the per-(l, m) `kernel`, falls back to
# scalewise on scalar outcomes (T_m == 1) when `method_name` is
# given, and computes `clfsr_curves` (the per-variant lfsr derived
# from `fit$mu` / `fit$mu2`) in one place.
#
# `kernel(fit, l, m, T_m, level)` returns a list with
# `effect_estimate`, `credible_band`, `lfsr`. `method_name` and
# `kind` are used only to format the scalar-fallback warning.
.smoother_loop <- function(fit, level, kernel,
                           method_name = NULL, kind = NULL) {
  meta <- fit$dwt_meta
  M    <- length(meta$T_basis)
  L    <- nrow(fit$alpha)

  effect_curves  <- vector("list", M)
  credible_bands <- vector("list", M)
  lfsr_curves    <- vector("list", M)
  clfsr_curves   <- vector("list", M)
  scalar_fallback <- NULL

  for (m in seq_len(M)) {
    T_m <- meta$T_basis[m]
    effect_curves[[m]]  <- vector("list", L)
    credible_bands[[m]] <- vector("list", L)
    lfsr_curves[[m]]    <- vector("list", L)
    clfsr_curves[[m]]   <- vector("list", L)

    if (T_m == 1L && !is.null(method_name)) {
      warning_message(sprintf(
        "method = '%s' is a %s smoother and adds no power for outcome %d (T_m = 1, scalar). Falling back to method = 'scalewise' for that outcome.",
        method_name, kind, m), style = "hint")
      if (is.null(scalar_fallback))
        scalar_fallback <- .post_smooth_scalewise(fit, level,
                                                  threshold_factor = 1)
      effect_curves[[m]]  <- scalar_fallback$effect_curves[[m]]
      credible_bands[[m]] <- scalar_fallback$credible_bands[[m]]
      lfsr_curves[[m]]    <- scalar_fallback$lfsr_curves[[m]]
      clfsr_curves[[m]]   <- scalar_fallback$clfsr_curves[[m]]
      next
    }

    for (l in seq_len(L)) {
      out <- kernel(fit, l, m, T_m, level)
      effect_curves[[m]][[l]]  <- out$effect_estimate
      credible_bands[[m]][[l]] <- out$credible_band
      lfsr_curves[[m]][[l]]    <- out$lfsr
      clfsr_curves[[m]][[l]]   <- lfsr_from_gaussian(
        fit$mu[[l]][[m]],
        sqrt(pmax(fit$mu2[[l]][[m]] - fit$mu[[l]][[m]]^2, 0)))
    }
  }

  list(effect_curves  = effect_curves,
       credible_bands = credible_bands,
       lfsr_curves    = lfsr_curves,
       clfsr_curves   = clfsr_curves)
}

# ---- scalewise -----------------------------------------------------

.post_smooth_scalewise <- function(fit, level, threshold_factor) {
  z_crit <- stats::qnorm((1 + level) / 2)
  meta   <- fit$dwt_meta
  kernel <- function(fit, l, m, T_m, level) {
    # Alpha-weighted aggregate across variables (no lead-variant
    # tie-break). Mean: sum_j alpha[l, j] * mu[l, j, t]; second
    # moment: sum_j alpha[l, j] * mu2[l, j, t]. Variance via the
    # law of total variance.
    mean_w <- as.numeric(fit$alpha[l, ] %*% fit$mu[[l]][[m]])
    mu2_w  <- as.numeric(fit$alpha[l, ] %*% fit$mu2[[l]][[m]])
    var_w  <- pmax(mu2_w - mean_w^2, 0)

    shrunk_w <- if (T_m == 1L) mean_w else
      scalewise_soft_threshold(mean_w, sd = sqrt(var_w),
        scale_index = meta$scale_index[[m]],
        T_basis     = T_m,
        factor      = threshold_factor)

    effect <- dwt_invert_packed(matrix(shrunk_w, nrow = 1L), meta, m)
    # Position-space variance: Var(pos[t]) = sum_k W^T_{t,k}^2 *
    # var_w[k]. `abs(invert_dwt(sqrt(var_w)))` would mistake a
    # linear combination of sds for the position sd.
    sd_pos <- sqrt(mf_invert_variance_curve(
      var_w,
      T_basis       = T_m,
      filter_number = meta$wavelet_filter %||% 10L,
      family        = meta$wavelet_family %||% "DaubLeAsymm",
      wavelet_scale = meta$wavelet_scale[[m]]))
    list(effect_estimate = effect,
         credible_band   = cbind(effect - z_crit * sd_pos,
                                 effect + z_crit * sd_pos),
         lfsr            = lfsr_from_gaussian(effect, sd_pos))
  }
  .smoother_loop(fit, level, kernel)
}

# ---- TI: cycle-spinning translation-invariant wavelet denoising ----

.post_smooth_ti <- function(fit, level, wavelet_filter, wavelet_family,
                            ...) {
  z_crit <- stats::qnorm((1 + level) / 2)
  Y_pos_cache <- vector("list", length(fit$dwt_meta$T_basis))
  kernel <- function(fit, l, m, T_m, level) {
    if (is.null(Y_pos_cache[[m]]))
      Y_pos_cache[[m]] <<- .iso_response_pos(fit, m)
    iso <- Y_pos_cache[[m]] - .other_effects_pos(fit, m, exclude = l)
    out <- univariate_ti_regression(iso, fit$X_eff[[l]],
                                    wavelet_filter, wavelet_family,
                                    z_crit, ...)
    list(effect_estimate = out$effect_estimate,
         credible_band   = cbind(out$cred_band[2L, ], out$cred_band[1L, ]),
         lfsr            = lfsr_from_gaussian(out$effect_estimate,
                                              out$fitted_sd))
  }
  .smoother_loop(fit, level, kernel,
                 method_name = "TI", kind = "wavelet")
}

# Position-space response. Reads the post-remap `Y_grid[[m]]`
# attached at fit time (n x T_basis, raw Y units, padded grid).
# When `Y_grid` is absent (the user fit with
# `attach_smoothing_inputs = FALSE`), `mf_post_smooth(fit, Y, ...)`
# attaches it on the fly before reaching here.
.iso_response_pos <- function(fit, m) {
  Y_pos <- fit$Y_grid[[m]]
  pos_m <- fit$dwt_meta$pos[[m]]
  if (ncol(Y_pos) > length(pos_m))
    Y_pos <- Y_pos[, seq_along(pos_m), drop = FALSE]
  Y_pos
}

# Sum_{l != exclude} X_eff[l] * (alpha-weighted inverse-DWT of mu_l)
# in raw-Y / raw-X units. `X_eff[[l]] = X_raw %*% alpha[l, ]`
# replaces the previous lead-variant column. The aggregate effect
# curve for outcome m, effect l is the alpha-weighted sum of the
# per-variable wavelet posteriors `sum_j alpha[l, j] * mu[l, j, ·]`,
# inverse-DWT'd, then multiplied by `csd_Y` (per-position) so the
# result lives in raw-Y units.
.other_effects_pos <- function(fit, m, exclude) {
  meta <- fit$dwt_meta
  L    <- nrow(fit$alpha)
  T_pos <- length(meta$pos[[m]])
  csd_m <- meta$column_scale[[m]]
  out  <- matrix(0, nrow = length(fit$X_eff[[1]]), ncol = T_pos)
  for (l in seq_len(L)) {
    if (l == exclude) next
    # Alpha-weighted aggregate per-position wavelet coefficient
    # across variables. `mu[[l]][[m]]` is p x T_basis.
    mu_w_agg <- as.numeric(fit$alpha[l, ] %*% fit$mu[[l]][[m]])
    eff_pos  <- dwt_invert_packed(matrix(mu_w_agg, nrow = 1L),
                                      meta, m)
    if (length(eff_pos) > T_pos) eff_pos <- eff_pos[seq_len(T_pos)]
    eff_pos_raw <- eff_pos * csd_m[seq_along(eff_pos)]
    out <- out + outer(fit$X_eff[[l]], eff_pos_raw)
  }
  out
}

# Stationary-wavelet regression of one outcome's isolated
# response on a per-effect predictor (an alpha-weighted X
# aggregate in production, or any single-column regressor), with
# scalewise ash shrinkage of the wavelet coefficients and a
# cycle-spinning average for the point estimate. Returns
# `effect_estimate` (length T_pos) and `cred_band` (2 x T_pos:
# row 1 "up", row 2 "low").
univariate_ti_regression <- function(Y_pos, x_eff,
                                      filter_number, family,
                                      z_crit,
                                      scaling = c("per_scale", "uniform"),
                                      ...) {
  scaling <- match.arg(scaling)

  # `per_scale`: column-scale Y and run a per-scale ash loop with
  # `nullweight = 30`. `uniform`: skip Y scaling, run a single ash
  # over all D coefficients with `nullweight = 3`. Both modes
  # forward `...` to `ashr::ash` so callers can override defaults
  # (e.g., `nullweight`, `mixcompdist`).
  if (scaling == "per_scale") {
    Y_sc  <- col_scale(Y_pos, center = TRUE, scale = TRUE)
    csd_Y <- attr(Y_sc, "scaled:scale")
  } else {
    Y_sc  <- Y_pos
    csd_Y <- rep(1, ncol(Y_pos))
  }
  x_mat <- col_scale(matrix(x_eff, ncol = 1L),
                     center = TRUE, scale = TRUE)
  csd_X <- attr(x_mat, "scaled:scale")
  x_sc  <- as.numeric(x_mat)

  # Stationary wavelet transform of each row of Y_sc.
  T_pos <- ncol(Y_sc)
  dummy <- wd(Y_sc[1L, ], type = "station",
              filter.number = filter_number, family = family)
  Y_f <- do.call(rbind, lapply(seq_len(nrow(Y_sc)), function(i)
    wd(Y_sc[i, ], type = "station",
       filter.number = filter_number, family = family)$D))
  Y_c <- do.call(rbind, lapply(seq_len(nrow(Y_sc)), function(i)
    wd(Y_sc[i, ], type = "station",
       filter.number = filter_number, family = family)$C))

  # Univariate regression of each wavelet column on x_sc.
  reg_d  <- compute_marginal_bhat_shat(matrix(x_sc, ncol = 1L), Y_f)
  reg_c  <- compute_marginal_bhat_shat(matrix(x_sc, ncol = 1L), Y_c)
  bhat_d <- as.numeric(reg_d$Bhat)
  shat_d <- as.numeric(reg_d$Shat)
  bhat_c <- as.numeric(reg_c$Bhat)
  shat_c <- as.numeric(reg_c$Shat)

  fl       <- dummy$fl.dbase$first.last.d
  n_basis  <- 2L^nlevelsWT(dummy)
  K        <- nrow(fl)
  wd_post  <- numeric(length(bhat_d))
  wd_var   <- numeric(length(bhat_d))
  if (scaling == "per_scale") {
    ash_d_args <- utils::modifyList(
      list(nullweight = 30, mixcompdist = "normal"), list(...))
    for (s in seq_len(K)) {
      first_s   <- fl[s, 1L]
      offset_s  <- fl[s, 3L]
      idx <- (offset_s + 1L - first_s):(offset_s + n_basis - first_s)
      t_ash <- do.call(ash,
                       c(list(bhat_d[idx], shat_d[idx]), ash_d_args))
      wd_post[idx] <- t_ash$result$PosteriorMean
      wd_var[idx]  <- t_ash$result$PosteriorSD^2
    }
  } else {
    ash_d_args <- utils::modifyList(
      list(nullweight = 3, mixcompdist = "normal"), list(...))
    t_ash <- do.call(ash, c(list(bhat_d, shat_d), ash_d_args))
    wd_post <- t_ash$result$PosteriorMean
    wd_var  <- t_ash$result$PosteriorSD^2
  }
  ash_c_args <- utils::modifyList(
    list(nullweight = 3, mixcompdist = "normal"), list(...))
  t_ash_c <- do.call(ash, c(list(bhat_c, shat_c), ash_c_args))

  # Cycle-spinning average via av.basis.
  dummy$D <- wd_post
  dummy$C <- t_ash_c$result$PosteriorMean
  mywst   <- convert(dummy)
  fitted_func_sc <- av.basis(mywst,
                             level = dummy$nlevels - 1L,
                             ix1 = 0, ix2 = 1,
                             filter = mywst$filter)

  # Exact pointwise variance via squared-filter wd / wst /
  # av-basis pipeline. The variance basis uses filter 10 /
  # DaubLeAsymm regardless of the smoothing filter.
  var_wd <- wd_variance(rep(0, T_pos))
  var_wd$D <- wd_var
  fitted_var_sc <- av_basis_variance(wst_variance(var_wd))

  # Unscale: recover the position-space estimate on the original
  # Y / X scale.
  fitted_func <- fitted_func_sc * csd_Y / csd_X
  fitted_var  <- fitted_var_sc  * (csd_Y / csd_X)^2
  fitted_sd   <- sqrt(pmax(fitted_var, 0))

  cred_band <- rbind(up  = fitted_func + z_crit * fitted_sd,
                     low = fitted_func - z_crit * fitted_sd)
  list(effect_estimate = fitted_func, cred_band = cred_band,
       fitted_sd = fitted_sd)
}

# ---- HMM denoising -------------------------------------------------

.post_smooth_hmm <- function(fit, level, halfK, ...) {
  z_crit <- stats::qnorm((1 + level) / 2)
  Y_pos_cache <- vector("list", length(fit$dwt_meta$T_basis))
  kernel <- function(fit, l, m, T_m, level) {
    if (is.null(Y_pos_cache[[m]]))
      Y_pos_cache[[m]] <<- .iso_response_pos(fit, m)
    iso <- Y_pos_cache[[m]] - .other_effects_pos(fit, m, exclude = l)
    out <- mf_univariate_hmm_regression(
      Y     = iso,
      X     = matrix(fit$X_eff[[l]], ncol = 1L),
      halfK = halfK, ...)
    list(effect_estimate = out$effect_estimate,
         credible_band   = cbind(
           out$effect_estimate - z_crit * out$effect_sd,
           out$effect_estimate + z_crit * out$effect_sd),
         lfsr            = out$lfsr)
  }
  .smoother_loop(fit, level, kernel,
                 method_name = "HMM", kind = "position-space")
}

# HMM helpers live in R/post_smooth_hmm.R: mf_fit_hmm() and
# mf_univariate_hmm_regression().

# --- helpers --------------------------------------------------------

# Soft-threshold the packed wavelet vector per scale at
# `factor * sd * sqrt(2 log T)`.
scalewise_soft_threshold <- function(coef_vec, sd, scale_index,
                                      T_basis, factor) {
  out <- coef_vec
  T_eff <- T_basis
  for (idx in scale_index) {
    if (length(idx) == 0L) next
    sigma <- mean(sd[idx], na.rm = TRUE)
    if (!is.finite(sigma) || sigma == 0) next
    thr <- factor * sigma * sqrt(2 * log(T_eff))
    x <- out[idx]
    out[idx] <- sign(x) * pmax(abs(x) - thr, 0)
  }
  out
}

# Use the existing mf_invert_dwt machinery to bring a length-T_basis
# packed wavelet vector back to position space.
dwt_invert_packed <- function(D_packed_row, meta, m) {
  inverted <- mf_invert_dwt(
    D_packed      = D_packed_row,
    column_center = rep(0, ncol(D_packed_row)),
    column_scale  = rep(1, ncol(D_packed_row)),
    wavelet_center = rep(0, ncol(D_packed_row)),
    wavelet_scale  = meta$wavelet_scale[[m]] %||% rep(1, ncol(D_packed_row)),
    filter_number = meta$wavelet_filter %||% 10L,
    family        = meta$wavelet_family %||% "DaubLeAsymm")
  as.numeric(inverted)
}
