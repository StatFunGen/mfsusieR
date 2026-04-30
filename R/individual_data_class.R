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
#' @param wavelet_basis_order integer, selects the wavelet basis
#'   member within `wavelet_family` (number of vanishing moments
#'   for Daubechies families). Forwarded to `wavethresh::wd`'s
#'   `filter.number`. See `filter.select`.
#' @param wavelet_family character, see `wd`.
#' @param standardize logical, scale `X` columns to unit variance.
#' @param intercept logical, center `X` columns to mean zero.
#' @param wavelet_magnitude_cutoff non-negative numeric. Wavelet
#'   columns whose `median(|column|)` is at or below this cutoff
#'   are zeroed and treated as uninformative.
#' @param wavelet_standardize logical. When `TRUE`, centers and
#'   scales each packed wavelet coefficient column before
#'   downstream steps. The stored wavelet centers and scales are
#'   used by inverse-DWT methods to return fitted values and
#'   coefficients to the original response scale.
#' @param wavelet_qnorm logical. When `TRUE`, applies a column-wise
#'   rank-based normal quantile transform to the wavelet-domain
#'   response before downstream steps.
#' @param verbose logical, forwarded to `mf_dwt`'s remap step.
#' @return list of class `c("mf_individual", "individual")`.
#' @references
#' Manuscript: `methods/online_method.tex`
#' (per-outcome wavelet decomposition + IBSS dispatch).
#' @keywords internal
#' @noRd
create_mf_individual <- function(X,
                                 Y,
                                 pos                      = NULL,
                                 max_padded_log2          = 10,
                                 wavelet_basis_order      = 10,
                                 wavelet_family           = "DaubLeAsymm",
                                 standardize              = TRUE,
                                 intercept                = TRUE,
                                 wavelet_magnitude_cutoff = 0,
                                 wavelet_standardize      = TRUE,
                                 wavelet_qnorm            = FALSE,
                                 verbose                  = FALSE) {
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
  wavelet_center <- vector("list", M)
  wavelet_scale  <- vector("list", M)
  lowc_idx      <- vector("list", M)

  for (m in seq_len(M)) {
    out <- mf_dwt(Y[[m]],
                  pos[[m]],
                  max_padded_log2 = max_padded_log2,
                  filter_number   = wavelet_basis_order,
                  family          = wavelet_family,
                  wavelet_standardize = wavelet_standardize,
                  verbose         = verbose)
    D[[m]]             <- out$D
    scale_index[[m]]   <- out$scale_index
    T_basis[m]        <- out$T_basis
    pos_grid[[m]]      <- out$pos
    column_center[[m]] <- out$column_center
    column_scale[[m]]  <- out$column_scale
    wavelet_center[[m]] <- out$wavelet_center
    wavelet_scale[[m]]  <- out$wavelet_scale

    # Functional-only preprocessing. Quantile-norm, low-count
    # masking, and position remap each operate on a multi-column
    # wavelet matrix and short-circuit to no-op for scalar
    # outcomes (`T_basis[m] == 1L`): the scalar `D[[m]]` is
    # already its own representation, has no inter-column
    # ordering to threshold, and lives on the original (single-
    # position) grid.
    lowc_idx[[m]]   <- integer(0)
    Y_remapped[[m]] <- Y[[m]]
    if (T_basis[m] > 1L) {
      if (wavelet_qnorm) D[[m]] <- mf_quantile_normalize(D[[m]])
      lowc_idx[[m]] <- mf_low_count_indices(D[[m]],
                                            threshold = wavelet_magnitude_cutoff)
      if (length(lowc_idx[[m]]) > 0L) {
        D[[m]][, lowc_idx[[m]]] <- 0
      }
      Y_remapped[[m]] <- remap_data(Y[[m]], pos[[m]],
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
    wavelet_standardize = wavelet_standardize,
    wavelet_meta = list(
      filter_number = wavelet_basis_order,
      family        = wavelet_family,
      column_center = column_center,
      column_scale  = column_scale,
      wavelet_center = wavelet_center,
      wavelet_scale  = wavelet_scale,
      X_center      = X_center,
      X_scale       = X_scale
    )
  )
  structure(obj, class = c("mf_individual", "individual"))
}
