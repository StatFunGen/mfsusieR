# `mf_individual` data class.
#
# The internal data class threaded through `susieR::susie_workhorse`
# via S3 method dispatch. Holds the SNP matrix `X`, per-modality
# remapped responses, the per-modality DWT cache, scale-index lists,
# and (when `save_residuals = TRUE`) the per-modality residuals
# filled in by `ibss_finalize.mf_individual`.
#
# Univariate (`T_m = 1`), single-modality (`M = 1`), and ragged
# multi-modality inputs all flow through this one class with shape
# branches, mirroring the `mvsusieR/refactor-s3` paradigm where one
# `mv_individual` covers `R = 1` and `R > 1`.
#
# `create_mf_individual` is internal. `mfsusie()` is the only caller.

#' Construct the internal `mf_individual` data class
#'
#' Validates inputs, centers and scales `X`, runs the per-modality
#' wavelet pipeline (`mf_dwt`), and packages everything the IBSS
#' loop needs into one S3 object of class
#' `c("mf_individual", "individual")`.
#'
#' @param X numeric matrix `n x J`.
#' @param Y list of length M; each element a numeric matrix
#'   `n x T_m`. Ragged `T_m` across modalities is permitted; an
#'   element may have `T_m = 1` (univariate, no DWT).
#' @param pos optional list of length M; each element a numeric
#'   vector of length `T_m` recording sampling positions. Defaults
#'   to `seq_len(T_m)` per modality.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length per modality. Default 10.
#' @param wavelet_filter_number integer, see
#'   `wavethresh::filter.select`.
#' @param wavelet_family character, see `wavethresh::wd`.
#' @param standardize logical, scale `X` columns to unit variance.
#' @param intercept logical, center `X` columns to mean zero.
#' @param save_residuals logical, allocate `residuals` storage on
#'   the object so post-processors can read residuals from the fit
#'   without the caller passing `(X, Y)` again. Set to `FALSE` for
#'   memory-bounded fits; the post-processors then take `(X, Y)`
#'   explicitly. See `inst/notes/refactor-discipline.md` for the
#'   rationale.
#' @param verbose logical, forwarded to `mf_dwt`'s remap step.
#' @return list of class `c("mf_individual", "individual")`.
#' @references
#' Manuscript: `methods/online_method.tex`
#' (per-modality wavelet decomposition + IBSS dispatch).
#' @keywords internal
#' @noRd
create_mf_individual <- function(X,
                                 Y,
                                 pos                   = NULL,
                                 max_padded_log2       = 10,
                                 wavelet_filter_number = 10,
                                 wavelet_family        = "DaubLeAsymm",
                                 standardize           = TRUE,
                                 intercept             = TRUE,
                                 save_residuals        = TRUE,
                                 verbose               = FALSE) {
  # ---- Input validation --------------------------------------------------
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  if (!is.list(Y) || length(Y) == 0L) {
    stop("`Y` must be a non-empty list of numeric matrices, ",
         "one per modality.")
  }
  M <- length(Y)
  n <- nrow(X)
  J <- ncol(X)

  Y <- lapply(Y, function(y) {
    if (is.null(dim(y))) y <- matrix(y, ncol = 1)
    if (!is.numeric(y)) stop("Every modality matrix in `Y` must be numeric.")
    y
  })

  for (m in seq_len(M)) {
    if (nrow(Y[[m]]) != n) {
      stop(sprintf(
        "`Y[[%d]]` has %d rows but `X` has %d rows; they must match.",
        m, nrow(Y[[m]]), n))
    }
  }

  if (is.null(pos)) {
    pos <- lapply(Y, function(y) seq_len(ncol(y)))
  } else {
    if (!is.list(pos) || length(pos) != M) {
      stop(sprintf(
        "`pos` must be a list of length %d (one position vector per modality).",
        M))
    }
    for (m in seq_len(M)) {
      if (length(pos[[m]]) != ncol(Y[[m]])) {
        stop(sprintf(
          "`pos[[%d]]` has length %d but `Y[[%d]]` has %d columns.",
          m, length(pos[[m]]), m, ncol(Y[[m]])))
      }
    }
  }

  # ---- Center and scale X ------------------------------------------------
  X_center <- if (intercept) colMeans(X) else rep(0, J)
  X_centered <- if (intercept) sweep(X, 2, X_center, "-") else X

  if (standardize) {
    X_scale <- matrixStats::colSds(X_centered, center = rep(0, J))
    X_scale[X_scale == 0] <- 1
    X_processed <- sweep(X_centered, 2, X_scale, "/")
  } else {
    X_scale <- rep(1, J)
    X_processed <- X_centered
  }

  # ---- Per-modality wavelet pipeline -------------------------------------
  Y_remapped    <- vector("list", M)
  D             <- vector("list", M)
  scale_index   <- vector("list", M)
  T_padded      <- integer(M)
  pos_grid      <- vector("list", M)
  column_center <- vector("list", M)
  column_scale  <- vector("list", M)

  for (m in seq_len(M)) {
    out <- mf_dwt(Y[[m]],
                  pos[[m]],
                  max_padded_log2 = max_padded_log2,
                  filter_number   = wavelet_filter_number,
                  family          = wavelet_family,
                  verbose         = verbose)
    D[[m]]             <- out$D
    scale_index[[m]]   <- out$scale_index
    T_padded[m]        <- out$T_padded
    pos_grid[[m]]      <- out$pos
    column_center[[m]] <- out$column_center
    column_scale[[m]]  <- out$column_scale
    # Reconstruct the post-remap, post-colScale Y from the wavelet
    # cache lazily when needed. For now, store a remapped Y at the
    # padded grid for residual computation downstream.
    Y_remapped[[m]] <- if (T_padded[m] == 1L) {
      Y[[m]]
    } else {
      remap_data(Y[[m]], pos[[m]],
                 verbose   = FALSE,
                 max_scale = max_padded_log2)$Y
    }
  }

  # ---- Assemble ----------------------------------------------------------
  obj <- list(
    X            = X_processed,
    Y            = Y_remapped,
    pos          = pos_grid,
    D            = D,
    scale_index  = scale_index,
    T_padded     = T_padded,
    csd_X        = X_scale,
    predictor_weights = colSums(X_processed^2),  # X'X diagonal; cached
                                                 # per the susieR / mvsusieR
                                                 # pattern, used by SER stats
                                                 # in compute_ser_statistics
                                                 # without per-iteration
                                                 # recompute.
    n            = n,
    J            = J,
    p            = J,    # alias for susieR helpers that read `data$p`
    M            = M,
    residuals    = if (save_residuals) vector("list", M) else NULL,
    intercept    = intercept,
    standardize  = standardize,
    wavelet_meta = list(
      filter_number = wavelet_filter_number,
      family        = wavelet_family,
      column_center = column_center,
      column_scale  = column_scale,
      X_center      = X_center,
      X_scale       = X_scale
    )
  )
  structure(obj, class = c("mf_individual", "individual"))
}
