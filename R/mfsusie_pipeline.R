# Wavelet-level association screen and multi-region pipeline for mfSuSiE.
#
# Lifecycle: user-facing pipeline layer. Called after ordinary upstream
# data QC. Runs per-region wavelet association tests (mfsusie_screen),
# then dispatches passing regions to mfsusie() and returns a unified
# result list plus an analysis funnel.
#
# OpenSpec: add-wavelet-pipeline-and-multiregion

# ── Helpers ──────────────────────────────────────────────────────────────────

#' Test whether an mfsusie object is an empty (screen-rejected / failed) result
#'
#' @param object an object with class `mfsusie`.
#' @return logical scalar.
#' @keywords internal
#' @noRd
mfsusie_is_empty <- function(object) isTRUE(object[["is_empty"]])


# ── mfsusie_empty_result ─────────────────────────────────────────────────────

#' Structured empty mfsusie result
#'
#' Returns an `mfsusie` object whose model fields are NULL / empty but whose
#' class and top-level field names are identical to a fitted result. Downstream
#' code that checks `$sets$cs`, `$pip`, or calls S3 methods registered on
#' `mfsusie` will not error (they return zero-length or `invisible(NULL)`).
#'
#' @param region_id character scalar; the region identifier, or `NULL`.
#' @param status character scalar; one of `"screen_rejected"` or
#'   `"fit_failed"`.
#' @param min_pval numeric scalar; the minimum wavelet association p-value
#'   observed during screening, or `NA_real_` when the screen was not run.
#' @return An object of class `c("mfsusie", "susie")` with empty model fields
#'   and a `$screen_result` sub-list recording `region_id`, `status`, and
#'   `min_pval`.
#' @export
mfsusie_empty_result <- function(region_id = NULL,
                                 status    = "screen_rejected",
                                 min_pval  = NA_real_) {
  structure(
    list(
      is_empty     = TRUE,
      alpha        = NULL,
      mu           = NULL,
      mu2          = NULL,
      lbf          = NULL,
      lbf_variable = NULL,
      V            = NULL,
      sigma2       = NULL,
      pi           = NULL,
      pip          = numeric(0),
      sets         = list(
        cs       = list(),
        purity   = list(),
        cs_index = integer(0)
      ),
      converged    = NA,
      niter        = 0L,
      screen_result = list(
        region_id = region_id,
        status    = status,
        min_pval  = min_pval
      )
    ),
    class = c("mfsusie", "susie")
  )
}


# ── mfsusie_screen ───────────────────────────────────────────────────────────

#' Wavelet-level marginal association screen for a single region
#'
#' Applies a per-outcome wavelet transform (optionally followed by quantile
#' normalization) and computes the minimum marginal association p-value over
#' all (wavelet coefficient, outcome, SNP) triples. The computation is
#' equivalent to running `lm(wavelet_Y[, j] ~ X[, k])` for every triple and
#' taking the minimum p-value, but uses the closed-form Pearson-correlation
#' t-statistic for efficiency.
#'
#' Scalar outcomes (`T_m = 1`) are excluded from the screen because they have
#' no wavelet expansion. When all outcomes are scalar, `min_pval` is `NA` and
#' `passed` is `TRUE` (the region is allowed through).
#'
#' @param X numeric matrix `n x p`; genotype dosages or other predictors.
#'   Need not be pre-centered.
#' @param Y list of length `M`; each element a numeric matrix `n x T_m` (or
#'   a length-`n` vector when `T_m = 1`).
#' @param pos list of length `M`; each element a numeric vector of length
#'   `T_m` recording genomic positions. If `NULL`, defaults to
#'   `seq_len(T_m)` per outcome.
#' @param wavelet_qnorm logical. When `TRUE` (default), applies
#'   `mf_quantile_normalize()` to the packed wavelet matrix before computing
#'   correlations. Must match the `wavelet_qnorm` argument used in
#'   `mfsusie()` to keep the tested and fitted representations consistent.
#' @param p_cutoff numeric; regions with `min_pval < p_cutoff` are marked as
#'   passing. Default `1e-5`.
#' @param region_id character scalar or `NULL`; passed through to the return
#'   value for tracking.
#' @param filter_number integer; wavelet filter number forwarded to `mf_dwt`.
#'   Default `10`.
#' @param family character; wavelet family forwarded to `mf_dwt`. Default
#'   `"DaubLeAsymm"`.
#' @return A named list with elements:
#'   \describe{
#'     \item{`region_id`}{The `region_id` argument, passed through.}
#'     \item{`min_pval`}{Minimum marginal p-value over all tested triples, or
#'       `NA` when no functional outcomes are present.}
#'     \item{`n_tested`}{Integer; total (wavelet column x outcome x SNP)
#'       triples evaluated.}
#'     \item{`passed`}{Logical; `TRUE` when `min_pval < p_cutoff`.}
#'   }
#' @export
mfsusie_screen <- function(X, Y, pos         = NULL,
                            wavelet_qnorm     = TRUE,
                            p_cutoff          = 1e-5,
                            region_id         = NULL,
                            filter_number     = 10L,
                            family            = "DaubLeAsymm") {
  if (!is.matrix(X)) stop("`X` must be a numeric matrix.")
  if (!is.list(Y))   stop("`Y` must be a list of matrices.")

  n <- nrow(X)
  X_ctr <- scale(X, center = TRUE, scale = FALSE)

  all_pvals <- numeric(0)
  n_tested  <- 0L

  for (m in seq_along(Y)) {
    Y_m <- Y[[m]]
    if (!is.matrix(Y_m)) Y_m <- as.matrix(Y_m)
    T_m <- ncol(Y_m)
    if (T_m == 1L) next  # scalar outcomes have no wavelet expansion

    pos_m <- if (!is.null(pos) && length(pos) >= m) pos[[m]] else seq_len(T_m)

    dwt_res <- mf_dwt(Y_m, pos_m,
                      filter_number = filter_number,
                      family        = family)
    D_m <- dwt_res$D  # n x T_basis

    if (wavelet_qnorm) D_m <- mf_quantile_normalize(D_m)

    # Pearson r: T_basis x p matrix
    r <- tryCatch(
      cor(D_m, X_ctr, use = "pairwise.complete.obs"),
      error = function(e) NULL
    )
    if (is.null(r)) next

    # t-stat and two-sided p-value (df = n - 2)
    r2     <- r * r
    t_stat <- r * sqrt(pmax(0, (n - 2) / pmax(1e-15, 1 - r2)))
    pv     <- 2 * pt(-abs(t_stat), df = n - 2)

    all_pvals <- c(all_pvals, as.numeric(pv))
    n_tested  <- n_tested + length(pv)
  }

  if (length(all_pvals) == 0L) {
    return(list(region_id = region_id, min_pval = NA_real_,
                n_tested = 0L, passed = TRUE))
  }

  min_pval <- min(all_pvals, na.rm = TRUE)
  list(
    region_id = region_id,
    min_pval  = min_pval,
    n_tested  = n_tested,
    passed    = min_pval < p_cutoff
  )
}


# ── mfsusie_pipeline ─────────────────────────────────────────────────────────

#' Multi-region mfSuSiE pipeline with optional wavelet association pre-screen
#'
#' Accepts one or more genomic regions, optionally screens each with
#' `mfsusie_screen()`, fits `mfsusie()` on passing regions, and returns a
#' unified result list plus an analysis funnel recording how many regions and
#' cell types entered, passed, and produced credible sets.
#'
#' Each element of `regions` is a named list with fields:
#' \describe{
#'   \item{`X`}{numeric matrix `n x p`; genotype dosages.}
#'   \item{`Y`}{list of `M` numeric matrices `n x T_m`.}
#'   \item{`pos`}{list of `M` numeric vectors of length `T_m` (genomic
#'     positions). May be omitted; defaults to `seq_len(T_m)` per outcome.}
#'   \item{`region_id`}{optional character scalar; used for tracking and as
#'     the name in the returned `results` list. Falls back to the element
#'     name in `regions`, then to `"region_<index>"`.}
#' }
#'
#' The analysis funnel (`$funnel`) has one row per region and columns:
#' `region_id`, `status` (`"screen_rejected"` | `"fit_completed"` |
#' `"fit_failed"`), `min_pval`, `n_cs`, `n_hp_cs`, `cell_types_entered`,
#' `cell_types_active`.
#'
#' @param regions named or unnamed list; each element describes one region
#'   (see above).
#' @param use_wavelet_filter logical. When `TRUE` (default), runs
#'   `mfsusie_screen()` before fitting and skips regions that fail.
#' @param wavelet_qnorm logical. Forwarded to both `mfsusie_screen()` and
#'   `mfsusie()`. Default `TRUE`. Must be identical in both to ensure the
#'   tested and fitted wavelet representations are consistent.
#' @param p_cutoff numeric; association p-value threshold passed to
#'   `mfsusie_screen()`. Regions with `min_pval >= p_cutoff` are rejected.
#'   Default `1e-5`.
#' @param n_cores integer. Number of parallel workers. Values `> 1` use
#'   `parallel::mclapply` (forking; not supported on Windows). Default `1L`.
#' @param verbose logical. When `TRUE`, emits a message for each region
#'   reporting screen outcome and fit status. Default `TRUE`.
#' @param min_abs_corr numeric; purity threshold used only for counting HP-CS
#'   in the returned funnel. Default `0.8`. To control which CS mfsusie()
#'   reports at all, pass `min_abs_corr = <value>` via `...` (mfsusie default
#'   is `0.5`).
#' @param ... additional arguments forwarded to `mfsusie()` for all regions.
#' @return A named list with two elements:
#'   \describe{
#'     \item{`results`}{Named list (by `region_id`) of per-region outputs.
#'       Passing regions are `mfsusie` fit objects. Rejected or failed regions
#'       are `mfsusie_empty_result()` objects of the same class, with an
#'       additional `$screen_result` sub-list.}
#'     \item{`funnel`}{data.frame; analysis funnel. One row per region.}
#'   }
#' @export
mfsusie_pipeline <- function(regions,
                              use_wavelet_filter = TRUE,
                              wavelet_qnorm      = TRUE,
                              p_cutoff           = 1e-5,
                              n_cores            = 1L,
                              verbose            = TRUE,
                              min_abs_corr       = 0.8,
                              ...) {
  if (!is.list(regions) || length(regions) == 0L)
    stop("`regions` must be a non-empty list.")

  region_names <- names(regions)
  get_id <- function(i) {
    r <- regions[[i]]
    if (!is.null(r[["region_id"]])) return(as.character(r[["region_id"]]))
    if (!is.null(region_names) && nzchar(region_names[i])) return(region_names[i])
    paste0("region_", i)
  }
  ids <- vapply(seq_along(regions), get_id, character(1L))

  process_one <- function(i) {
    r   <- regions[[i]]
    rid <- ids[i]
    X   <- r[["X"]]
    Y   <- r[["Y"]]
    pos <- r[["pos"]]

    # Step 1 — optional wavelet filter
    screen_res <- if (isTRUE(use_wavelet_filter)) {
      tryCatch(
        mfsusie_screen(X, Y, pos,
                       wavelet_qnorm = wavelet_qnorm,
                       p_cutoff      = p_cutoff,
                       region_id     = rid),
        error = function(e)
          list(region_id = rid, min_pval = NA_real_,
               n_tested = 0L, passed = FALSE,
               screen_error = conditionMessage(e))
      )
    } else {
      list(region_id = rid, min_pval = NA_real_, n_tested = 0L, passed = TRUE)
    }

    if (!isTRUE(screen_res[["passed"]])) {
      if (verbose)
        message(sprintf("[%s] screen rejected (min_pval = %s)",
                        rid, format(screen_res[["min_pval"]], digits = 3)))
      return(list(
        result = mfsusie_empty_result(region_id = rid,
                                      status    = "screen_rejected",
                                      min_pval  = screen_res[["min_pval"]]),
        screen = screen_res,
        status = "screen_rejected"
      ))
    }

    if (verbose)
      message(sprintf("[%s] screen passed (min_pval = %s), fitting ...",
                      rid, format(screen_res[["min_pval"]], digits = 3)))

    # Step 2 — fit mfsusie
    fit <- tryCatch(
      mfsusie(X = X, Y = Y, pos = pos,
              wavelet_qnorm = wavelet_qnorm,
              ...),
      error = function(e) {
        if (verbose)
          message(sprintf("[%s] mfsusie() failed: %s", rid, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(fit)) {
      return(list(
        result = mfsusie_empty_result(region_id = rid,
                                      status    = "fit_failed",
                                      min_pval  = screen_res[["min_pval"]]),
        screen = screen_res,
        status = "fit_failed"
      ))
    }

    if (verbose)
      message(sprintf("[%s] fit complete (%d CS)", rid,
                      length(fit[["sets"]][["cs"]])))

    list(result = fit, screen = screen_res, status = "fit_completed")
  }

  # Run serial or parallel
  if (n_cores > 1L) {
    if (.Platform$OS.type == "windows") {
      warning("n_cores > 1 is not supported on Windows; falling back to n_cores = 1.")
      outs <- lapply(seq_along(regions), process_one)
    } else {
      outs <- parallel::mclapply(seq_along(regions), process_one,
                                 mc.cores = n_cores)
    }
  } else {
    outs <- lapply(seq_along(regions), process_one)
  }
  names(outs) <- ids

  # mclapply returns try-error for workers killed by SIGKILL / OOM.
  # Convert those to fit_failed empty results so the funnel loop never sees them.
  for (i in seq_along(outs)) {
    if (inherits(outs[[i]], "try-error")) {
      msg <- as.character(outs[[i]])
      if (verbose)
        message(sprintf("[%s] worker killed (OOM or SIGKILL): %s", ids[i], msg))
      outs[[i]] <- list(
        result = mfsusie_empty_result(region_id = ids[i],
                                      status    = "fit_failed",
                                      min_pval  = NA_real_),
        screen = list(region_id = ids[i], min_pval = NA_real_,
                      n_tested = 0L, passed = NA),
        status = "fit_failed"
      )
    }
  }

  # Assemble results
  results        <- lapply(outs, `[[`, "result")
  names(results) <- ids

  # Build funnel
  funnel_rows <- lapply(seq_along(outs), function(i) {
    out <- outs[[i]]
    fit <- out[["result"]]
    scr <- out[["screen"]]

    n_cs     <- 0L
    n_hp_cs  <- 0L
    ct_active <- NA_integer_

    if (out[["status"]] == "fit_completed" && !mfsusie_is_empty(fit)) {
      cs  <- fit[["sets"]][["cs"]]
      n_cs <- length(cs)

      pur <- fit[["sets"]][["purity"]]
      if (is.data.frame(pur)) {
        n_hp_cs <- sum(pur[["min.abs.corr"]] >= min_abs_corr, na.rm = TRUE)
      } else if (is.list(pur) && length(pur) > 0L) {
        n_hp_cs <- sum(vapply(seq_along(pur), function(k) {
          p <- pur[[k]]
          is.matrix(p) && nrow(p) > 0L && p[1L, 1L] >= min_abs_corr
        }, logical(1L)))
      }

      meta <- fit[["dwt_meta"]]
      if (!is.null(meta)) ct_active <- meta[["M"]]
    }

    data.frame(
      region_id          = ids[[i]],
      status             = out[["status"]],
      min_pval           = scr[["min_pval"]] %||% NA_real_,
      n_cs               = n_cs,
      n_hp_cs            = n_hp_cs,
      cell_types_entered = length(regions[[i]][["Y"]]),
      cell_types_active  = ct_active,
      stringsAsFactors   = FALSE
    )
  })
  funnel <- do.call(rbind, funnel_rows)

  list(results = results, funnel = funnel)
}
