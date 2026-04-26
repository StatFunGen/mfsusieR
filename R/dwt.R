# Per-modality forward and inverse DWT wrappers.
#
# `mf_dwt` runs the full pipeline for one modality: position remap to
# a power-of-two grid (when needed), per-position centering and
# scaling, row-wise DWT, and packing of the resulting D and C
# coefficients into a single `n x T_padded` matrix. The univariate
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
# (wavelet-based functional regression on a per-modality basis).

#' Forward DWT pipeline for one modality
#'
#' Applies, in order: `remap_data` (regrid to a power-of-two grid
#' when the input is irregular or has a non-power-of-two length),
#' `col_scale` (center + scale each position column), `dwt_matrix`
#' (row-wise wavelet transform), and packs the resulting detail (D)
#' and smoothing (C) coefficients into a single `n x T_padded`
#' matrix with C in the last column. Returns enough metadata for
#' `mf_invert_dwt` to recover the position-space curve.
#'
#' For univariate inputs (`T_m = 1`), the wavelet machinery is
#' skipped: the single column is centered and scaled, returned as a
#' length-1 "wavelet coefficient" with `scale_index = 1L`, and
#' `T_padded = 1L`. The inverse helper reverses the centering /
#' scaling without any wavelet operation.
#'
#' @param Y_m numeric matrix `n x T_m`, one curve per row.
#' @param pos_m numeric vector of length `T_m`, sampling positions.
#' @param max_padded_log2 integer, log2 cap on the post-remap grid
#'   length. Forwarded to `remap_data` and to `dwt_matrix`.
#' @param filter_number integer, see `filter.select`.
#' @param family character, wavelet family name. See
#'   `wd`.
#' @param verbose logical, forwarded to `remap_data`.
#' @return list with elements `D` (packed wavelet matrix `n x
#'   T_padded`), `scale_index` (integer vector from
#'   `gen_wavelet_indx`), `T_padded` (integer scalar), `pos`
#'   (post-remap position grid, length `T_padded`),
#'   `column_center`, `column_scale` (vectors of length `T_padded`,
#'   for inverse DWT), `family`, `filter_number`.
#' @keywords internal
#' @noRd
mf_dwt <- function(Y_m,
                   pos_m,
                   max_padded_log2  = 10,
                   filter_number    = 10,
                   family           = "DaubLeAsymm",
                   verbose          = FALSE) {
  if (!is.matrix(Y_m)) {
    Y_m <- as.matrix(Y_m)
  }
  T_m <- ncol(Y_m)

  # Univariate short-circuit. The "wavelet representation" of a
  # length-1 signal is the signal itself; we only need the
  # per-column center + sd scale (matching the functional path's
  # convention so the two paths are interchangeable at T_m = 1).
  # The branch exists because `wd` does not accept
  # length-1 inputs.
  if (T_m == 1) {
    cm  <- mean(Y_m, na.rm = TRUE)
    csd <- sd(as.numeric(Y_m), na.rm = TRUE)
    if (!is.finite(csd) || csd == 0) csd <- 1
    Y_scaled <- (Y_m - cm) / csd
    return(list(
      D             = Y_scaled,
      scale_index   = list(1L),
      T_padded      = 1L,
      pos           = pos_m,
      column_center = cm,
      column_scale  = csd,
      family        = family,
      filter_number = filter_number
    ))
  }

  # Functional path.
  remap <- remap_data(Y_m, pos_m,
                      verbose   = verbose,
                      max_scale = max_padded_log2)
  Y_remapped  <- remap$Y
  outing_grid <- remap$outing_grid

  # Per-column center + sd scale. The C2 / C3 fidelity tests pin
  # bit-identical wavelet output against the port-source pipeline.
  # For the C1 degeneracy contract against susieR, the test
  # pre-scales y on both sides so the mfsusieR-internal sd-scaling
  # is a no-op.
  Y_scaled      <- col_scale(Y_remapped, center = TRUE, scale = TRUE)
  column_center <- attr(Y_scaled, "scaled:center")
  column_scale  <- attr(Y_scaled, "scaled:scale")

  T_padded   <- ncol(Y_remapped)
  log2_T_pad <- log2(T_padded)
  W <- dwt_matrix(Y_scaled,
                  filter_number = filter_number,
                  family        = family,
                  max_scale     = max_padded_log2)
  D_packed    <- cbind(W$D, W$C)
  scale_index <- gen_wavelet_indx(log2_T_pad)

  list(
    D             = D_packed,
    scale_index   = scale_index,
    T_padded      = T_padded,
    pos           = outing_grid,
    column_center = column_center,
    column_scale  = column_scale,
    family        = family,
    filter_number = filter_number
  )
}

#' Inverse of `mf_dwt` for one modality
#'
#' Takes a packed wavelet representation `n x T_padded` (or a
#' length-`T_padded` vector for a single curve) and the
#' `column_center` / `column_scale` recorded by `mf_dwt`, and
#' returns the corresponding position-space curve(s). For
#' univariate inputs the inverse is a single multiply-add. For
#' functional inputs, the inverse threads through `wr`.
#'
#' @param D_packed numeric matrix `n x T_padded` or numeric vector
#'   of length `T_padded`. The detail coefficients occupy the
#'   first `T_padded - 1` columns; the smoothing coefficient
#'   occupies the last column (the convention `mf_dwt` produces).
#' @param column_center numeric of length `T_padded`, the
#'   per-position centers from `mf_dwt`.
#' @param column_scale numeric of length `T_padded`, the
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
                          filter_number = 10,
                          family        = "DaubLeAsymm") {
  if (is.null(dim(D_packed))) {
    D_packed <- matrix(D_packed, nrow = 1)
  }
  T_padded <- ncol(D_packed)
  n <- nrow(D_packed)

  # Univariate short-circuit.
  if (T_padded == 1) {
    return(D_packed * column_scale + column_center)
  }

  # Build a `wd` skeleton at the right length (filled with zeros);
  # we inject the D and C coefficients per row and call wr().
  template <- wd(rep(0, T_padded),
                             filter.number = filter_number,
                             family        = family,
                             min.scale     = log2(T_padded))

  Y_curves <- matrix(0, nrow = n, ncol = T_padded)
  for (i in seq_len(n)) {
    w <- template
    w$D <- D_packed[i, -T_padded]
    # wavethresh stores the C coefficients as a vector spanning all
    # scales; the coarsest (level 0) C lives at the last position.
    w$C[length(w$C)] <- D_packed[i, T_padded]
    Y_curves[i, ] <- wr(w)
  }

  # Reverse the per-position centering and scaling.
  sweep(sweep(Y_curves, 2, column_scale, "*"),
        2, column_center, "+")
}
