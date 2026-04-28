# `mf_individual` data class. Threaded through
# `susie_workhorse` via S3 dispatch. One class covers
# univariate (`T_m = 1`), single-outcome (`M = 1`), and ragged
# multi-outcome inputs.

#' Construct the internal `mf_individual` data class
#'
#' Validates inputs, centers and scales `X`, runs the per-outcome
#' wavelet pipeline (`mf_dwt`), and packages everything the IBSS
#' loop needs into one S3 object of class
#' `c("mf_individual", "individual")`.
#'
#' @param X numeric matrix `n x p`.
#' @param Y list of length M; each element a numeric matrix
#'   `n x T_m`. Ragged `T_m` across outcomes is permitted; an
#'   element may have `T_m = 1` (univariate, no DWT).
#' @param pos optional list of length M; each element a numeric
#'   vector of length `T_m` recording sampling positions. Defaults
#'   to `seq_len(T_m)` per outcome.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length per outcome. Default 10.
#' @param wavelet_filter_number integer, see
#'   `filter.select`.
#' @param wavelet_family character, see `wd`.
#' @param standardize logical, scale `X` columns to unit variance.
#' @param intercept logical, center `X` columns to mean zero.
#' @param verbose logical, forwarded to `mf_dwt`'s remap step.
#' @return list of class `c("mf_individual", "individual")`.
#' @references
#' Manuscript: `methods/online_method.tex`
#' (per-outcome wavelet decomposition + IBSS dispatch).
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
                                 low_count_filter      = 0,
                                 quantile_norm         = FALSE,
                                 verbose               = FALSE) {
  # ---- Input validation --------------------------------------------------
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }
  if (!is.list(Y) || length(Y) == 0L) {
    stop("`Y` must be a non-empty list of numeric matrices, ",
         "one per outcome.")
  }
  M <- length(Y)
  n <- nrow(X)
  p <- ncol(X)

  Y <- lapply(Y, function(y) {
    if (is.null(dim(y))) y <- matrix(y, ncol = 1)
    if (!is.numeric(y)) stop("Every outcome matrix in `Y` must be numeric.")
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
        "`pos` must be a list of length %d (one position vector per outcome).",
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
  X_center <- if (intercept) colMeans(X) else rep(0, p)
  X_centered <- if (intercept) sweep(X, 2, X_center, "-") else X

  if (standardize) {
    X_scale <- colSds(X_centered, center = rep(0, p))
    X_scale[X_scale == 0] <- 1
    X_processed <- sweep(X_centered, 2, X_scale, "/")
  } else {
    X_scale <- rep(1, p)
    X_processed <- X_centered
  }

  # ---- Per-outcome wavelet pipeline -------------------------------------
  Y_remapped    <- vector("list", M)
  D             <- vector("list", M)
  scale_index   <- vector("list", M)
  T_basis      <- integer(M)
  pos_grid      <- vector("list", M)
  column_center <- vector("list", M)
  column_scale  <- vector("list", M)
  lowc_idx      <- vector("list", M)

  for (m in seq_len(M)) {
    out <- mf_dwt(Y[[m]],
                  pos[[m]],
                  max_padded_log2 = max_padded_log2,
                  filter_number   = wavelet_filter_number,
                  family          = wavelet_family,
                  verbose         = verbose)
    D[[m]]             <- out$D
    scale_index[[m]]   <- out$scale_index
    T_basis[m]        <- out$T_basis
    pos_grid[[m]]      <- out$pos
    column_center[[m]] <- out$column_center
    column_scale[[m]]  <- out$column_scale

    # Optional preprocessing: column-wise rank-based normal quantile
    # transform on the wavelet-domain response. Applied BEFORE the
    # low-count index computation so the median check sees the
    # transformed scale.
    if (quantile_norm && T_basis[m] > 1L) {
      D[[m]] <- mf_quantile_normalize(D[[m]])
    }

    # Optional preprocessing: low-count column indices. With the
    # default threshold 0 the set is empty when every column has
    # at least one nonzero absolute value across samples.
    if (T_basis[m] > 1L) {
      lowc_idx[[m]] <- mf_low_count_indices(D[[m]],
                                            threshold = low_count_filter)
      if (length(lowc_idx[[m]]) > 0L) {
        D[[m]][, lowc_idx[[m]]] <- 0
      }
    } else {
      lowc_idx[[m]] <- integer(0)
    }

    # Cache the post-remap Y on the padded grid for downstream
    # residual computation; the wavelet representation lives in `D`.
    Y_remapped[[m]] <- if (T_basis[m] == 1L) {
      Y[[m]]
    } else {
      remap_data(Y[[m]], pos[[m]],
                 verbose   = FALSE,
                 max_scale = max_padded_log2)$Y
    }
  }

  # ---- Per-outcome complete-case indices ---------------------------------
  na_idx <- lapply(seq_len(M), function(m) which(complete.cases(Y[[m]])))
  xtx_diag_list <- lapply(seq_len(M), function(m)
    colSums(X_processed[na_idx[[m]], , drop = FALSE]^2))

  # ---- Assemble ----------------------------------------------------------
  obj <- list(
    X             = X_processed,
    Y             = Y_remapped,
    pos           = pos_grid,
    D             = D,
    scale_index   = scale_index,
    lowc_idx      = lowc_idx,
    T_basis       = T_basis,
    csd           = X_scale,
    xtx_diag      = colSums(X_processed^2),    # global; zero-predictor check
    xtx_diag_list = xtx_diag_list,
    na_idx        = na_idx,
    n             = n,
    p             = p,
    M             = M,
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
