# Per-outcome forward and inverse DWT wrappers.
#
# `mf_dwt` runs the full pipeline for one outcome: position remap to
# a power-of-two grid (when needed), per-position centering and
# scaling, row-wise DWT, and packing of the resulting D and C
# coefficients into a single `n x T_basis` matrix. The univariate
# (`T_m = 1`) case short-circuits the wavelet machinery and returns
# the centered / scaled column directly.
#
# `mf_invert_dwt` is the inverse: take a packed wavelet representation
# and a column-center / column-scale pair and reconstruct the
# position-space curve(s). Used by `predict.mfsusie`, `coef.mfsusie`,
# and the post-processors. Exported nowhere; the public surface goes
# through the S3 methods.
#
# Manuscript reference: methods/online_method.tex
# (wavelet-based functional regression on a per-outcome basis).

#' Forward DWT pipeline for one outcome
#'
#' Applies, in order: `remap_data` (regrid to a power-of-two grid
#' when the input is irregular or has a non-power-of-two length),
#' `col_scale` (center + scale each position column), `dwt_matrix`
#' (row-wise wavelet transform), and packs the resulting detail (D)
#' and smoothing (C) coefficients into a single `n x T_basis`
#' matrix with C in the last column. Returns enough metadata for
#' `mf_invert_dwt` to recover the position-space curve.
#'
#' For univariate inputs (`T_m = 1`), the wavelet machinery is
#' skipped: the single column is centered and scaled, returned as a
#' length-1 "wavelet coefficient" with `scale_index = 1L`, and
#' `T_basis = 1L`. The inverse helper reverses the centering /
#' scaling without any wavelet operation.
#'
#' @param Y_m numeric matrix `n x T_m`, one curve per row.
#' @param pos_m numeric vector of length `T_m`, sampling positions.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length. Forwarded to `remap_data` and to `dwt_matrix`.
#' @param filter_number integer, see `filter.select`.
#' @param family character, wavelet family name. See
#'   `wd`.
#' @param wavelet_standardize logical. When `TRUE`, centers and
#'   scales each packed wavelet coefficient column after the DWT.
#'   The returned `wavelet_center` and `wavelet_scale` make this
#'   linear preprocessing invertible.
#' @param verbose logical, forwarded to `remap_data`.
#' @return list with elements `D` (packed wavelet matrix `n x
#'   T_basis`), `scale_index` (integer vector from
#'   `gen_wavelet_indx`), `T_basis` (integer scalar), `pos`
#'   (post-remap position grid, length `T_basis`),
#'   `column_center`, `column_scale` (vectors of length `T_basis`,
#'   for inverse DWT), `wavelet_center`, `wavelet_scale`,
#'   `family`, `filter_number`.
#' @keywords internal
#' @noRd
mf_dwt <- function(Y_m,
                   pos_m,
                   max_padded_log2  = 10,
                   filter_number    = 10,
                   family           = "DaubLeAsymm",
                   wavelet_standardize = FALSE,
                   verbose          = FALSE) {
  if (!is.matrix(Y_m)) {
    Y_m <- as.matrix(Y_m)
  }
  T_m <- ncol(Y_m)

  # Position remap is only meaningful for `T_m > 1`. For a scalar
  # outcome the wavelet representation IS the signal, sampled at
  # the original (single) position.
  if (T_m == 1L) {
    Y_remapped  <- Y_m
    outing_grid <- pos_m
  } else {
    remap <- remap_data(Y_m, pos_m,
                        verbose   = verbose,
                        max_scale = max_padded_log2)
    Y_remapped  <- remap$Y
    outing_grid <- remap$outing_grid
  }

  # Per-column center + sd scale. col_scale floors zero-variance
  # columns at csd = 1, so the scalar / constant-column case is
  # handled the same way as any zero-variance functional column.
  Y_scaled      <- col_scale(Y_remapped, center = TRUE, scale = TRUE)
  column_center <- attr(Y_scaled, "scaled:center")
  column_scale  <- attr(Y_scaled, "scaled:scale")

  T_basis <- ncol(Y_remapped)

  # Forward DWT runs only for `T_basis > 1`; `wavethresh::wd`
  # rejects length-1 inputs. The scalar case takes the identity
  # path: `D` is the centered/scaled signal and `scale_index` is
  # the trivial one-group list, matching the multi-column shape.
  # Strip `col_scale`'s bookkeeping attributes from the scalar
  # `D_packed` so the field shape matches `cbind(W$D, W$C)` from
  # the functional branch (no attributes beyond `dim`).
  if (T_basis == 1L) {
    D_packed    <- matrix(as.numeric(Y_scaled),
                          nrow = nrow(Y_scaled), ncol = 1L)
    scale_index <- list(1L)
  } else {
    W <- dwt_matrix(Y_scaled,
                    filter_number = filter_number,
                    family        = family,
                    max_scale     = max_padded_log2)
    D_packed    <- cbind(W$D, W$C)
    scale_index <- gen_wavelet_indx(log2(T_basis))
  }

  wavelet_center <- rep(0, T_basis)
  wavelet_scale  <- rep(1, T_basis)
  if (wavelet_standardize && T_basis > 1L) {
    D_scaled       <- col_scale(D_packed, center = TRUE, scale = TRUE)
    wavelet_center <- attr(D_scaled, "scaled:center")
    wavelet_scale  <- attr(D_scaled, "scaled:scale")
    D_packed       <- unclass(D_scaled)
  }

  list(
    D             = D_packed,
    scale_index   = scale_index,
    T_basis       = T_basis,
    pos           = outing_grid,
    column_center = column_center,
    column_scale  = column_scale,
    wavelet_center = wavelet_center,
    wavelet_scale  = wavelet_scale,
    family        = family,
    filter_number = filter_number
  )
}

#' Inverse of `mf_dwt` for one outcome
#'
#' Takes a packed wavelet representation `n x T_basis` (or a
#' length-`T_basis` vector for a single curve) and the
#' `column_center` / `column_scale` recorded by `mf_dwt`, and
#' returns the corresponding position-space curve(s). For
#' univariate inputs the inverse is a single multiply-add. For
#' functional inputs, the inverse threads through `wr`.
#'
#' @param D_packed numeric matrix `n x T_basis` or numeric vector
#'   of length `T_basis`. The detail coefficients occupy the
#'   first `T_basis - 1` columns; the smoothing coefficient
#'   occupies the last column (the convention `mf_dwt` produces).
#' @param column_center numeric of length `T_basis`, the
#'   per-position centers from `mf_dwt`.
#' @param column_scale numeric of length `T_basis`, the
#'   per-position scales from `mf_dwt`.
#' @param filter_number integer, must match the value passed to
#'   `mf_dwt`.
#' @param family character, must match the value passed to
#'   `mf_dwt`.
#' @return numeric matrix of curves with the same shape as
#'   `D_packed`, in the position space and on the same grid as the
#'   input to `mf_dwt`.
#' @keywords internal
#' @noRd
mf_invert_dwt <- function(D_packed,
                          column_center,
                          column_scale,
                          wavelet_center = NULL,
                          wavelet_scale  = NULL,
                          filter_number = 10,
                          family        = "DaubLeAsymm") {
  if (is.null(dim(D_packed))) {
    D_packed <- matrix(D_packed, nrow = 1)
  }
  T_basis <- ncol(D_packed)
  n <- nrow(D_packed)
  if (is.null(wavelet_center)) wavelet_center <- rep(0, T_basis)
  if (is.null(wavelet_scale))  wavelet_scale  <- rep(1, T_basis)
  D_packed <- sweep(sweep(D_packed, 2, wavelet_scale, "*"),
                    2, wavelet_center, "+")

  # `wavethresh::wr` rejects length-1 inputs, so the scalar
  # case takes the identity path. The reverse-center+scale step
  # below would also handle this branch — we keep it explicit
  # because building a `wd` template at `T_basis = 1` errors.
  if (T_basis == 1L) {
    Y_curves <- D_packed
  } else {
    # Build a `wd` skeleton at the right length (filled with zeros);
    # we inject the D and C coefficients per row and call wr().
    template <- wd(rep(0, T_basis),
                               filter.number = filter_number,
                               family        = family,
                               min.scale     = log2(T_basis))

    Y_curves <- matrix(0, nrow = n, ncol = T_basis)
    for (i in seq_len(n)) {
      w <- template
      w$D <- D_packed[i, -T_basis]
      # wavethresh stores the C coefficients as a vector spanning all
      # scales; the coarsest (level 0) C lives at the last position.
      w$C[length(w$C)] <- D_packed[i, T_basis]
      Y_curves[i, ] <- wr(w)
    }
  }

  # Reverse the per-position centering and scaling.
  sweep(sweep(Y_curves, 2, column_scale, "*"),
        2, column_center, "+")
}
