# Wavelet-pathway utilities (internal).
#
# These helpers prepare functional response matrices for the
# discrete wavelet transform: column-wise centering and scaling,
# remapping of unevenly-spaced or non-power-of-two functional inputs
# onto a power-of-two grid via the Kovac-Silverman lifting scheme,
# scale-index generation, and a row-wise DWT wrapper around
# `wd`. None of the helpers are exported.

#' Test whether a numeric value is integer-valued
#'
#' @param x numeric vector.
#' @param tol numeric tolerance.
#' @return logical vector, `TRUE` where `x[i]` is within `tol` of an integer.
#' @keywords internal
#' @noRd
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

#' Column-wise centering and scaling
#'
#' Centers each column of `x` to mean zero and (optionally) scales to
#' unit standard deviation, missing values handled per column. Zero-
#' variance columns are left untouched (their `csd` entry is forced
#' to 1). The scaling factors are returned as attributes
#' (`scaled:center`, `scaled:scale`) so an inverse transform is
#' available.
#'
#' @param x numeric matrix.
#' @param center logical, center columns to mean zero.
#' @param scale logical, scale columns to unit variance.
#' @param add_attr logical, return scaling factors as attributes.
#' @param rows optional integer vector, subset rows before scaling.
#' @param cols optional integer vector, subset columns before scaling.
#' @return numeric matrix; same dimensions as the (subset) input.
#' @keywords internal
#' @noRd
col_scale <- function(x,
                      center   = TRUE,
                      scale    = TRUE,
                      add_attr = TRUE,
                      rows     = NULL,
                      cols     = NULL) {
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  cm <- colMeans(x, na.rm = TRUE)
  csd <- if (scale) {
    colSds(x, center = cm, na.rm = TRUE)
  } else {
    rep(1, length(cm))
  }
  csd[csd == 0] <- 1

  if (!center) {
    cm <- rep(0, length(cm))
  }
  x <- t((t(x) - cm) / csd)

  if (add_attr) {
    if (center) attr(x, "scaled:center") <- cm
    if (scale)  attr(x, "scaled:scale")  <- csd
    # The `d` attribute is algebraically (n - 1) per column once
    # zero-variance columns have had their csd replaced with 1
    # (which has happened above). The port source computes the
    # same quantity via a two-step `(a + b - a) / b` expression
    # that introduces ULP-level rounding noise; we compute the
    # simplified form directly. The C2 fidelity test on the `d`
    # attribute uses a relaxed tolerance (`1e-12`) to accept the
    # resulting machine-epsilon difference.
    attr(x, "d") <- rep(nrow(x) - 1, ncol(x))
  }
  x
}

#' Power-of-two grid size from a sample count
#'
#' Returns the next power of two above `n - 1`, capped at
#' `2^max_scale`. Used to size the regular grid that
#' `makegrid` interpolates onto.
#'
#' @param n integer, sample count along one curve.
#' @param max_scale integer, log2 cap on grid length.
#' @return integer grid length, a power of two.
#' @keywords internal
#' @noRd
power_of_two_gridn <- function(n, max_scale = 10) {
  min(2^max_scale, 2^(floor(log(n - 1, base = 2)) + 1))
}

#' Kovac-Silverman lifting-scheme interpolation onto a 2^k grid
#'
#' Interpolates a single observed curve `y` sampled at scaled
#' positions `bp` (in `[0, 1]`) onto a regular grid of length
#' `power_of_two_gridn(length(bp), max_scale)`. The interpolation
#' is the lifting scheme of Kovac and Silverman (2000), implemented
#' in `makegrid`.
#'
#' @param y numeric vector of observed values along one curve.
#' @param bp numeric vector of scaled positions in `[0, 1]` matching
#'   `y` element-wise.
#' @param gridn integer, regular-grid length. Computed by
#'   `power_of_two_gridn(length(bp), max_scale)` in `interpol_mat`
#'   and passed in to avoid recomputing once per row.
#' @return numeric vector on the regular grid.
#' @keywords internal
#' @noRd
interpol_ks <- function(y, bp, gridn) {
  makegrid(t = bp, y = y, gridn = gridn)$gridy
}

#' Interpolate a matrix of curves onto a power-of-two grid
#'
#' Applies `interpol_ks` row by row to remap a matrix `Y` of curves
#' sampled at irregular positions `pos` onto a regular grid of length
#' that is a power of two and at most `2^max_scale`. Returns the
#' interpolated `Y` and the regular grid expressed in the original
#' position units, ready for downstream use.
#'
#' Behaviour-preserving simplification of the port source's
#' `interpol_mat` plus its separate position-rescaling step. The
#' original chain rescaled the wavethresh `gridt` (in `[0, 1]`)
#' to a `[0, length(grid)]` integer scale, then back to position
#' units; we collapse the two rescalings into one linear
#' interpolation that produces the same numbers in fewer passes
#' over `pos`.
#'
#' @param Y numeric matrix `n x T`, one curve per row.
#' @param pos numeric vector of length `T`, sampling positions.
#' @param max_scale integer, log2 cap on grid length.
#' @return list with elements `Y` (interpolated matrix) and
#'   `outing_grid` (regular grid in `pos` units, length a power of
#'   two and at most `2^max_scale`).
#' @keywords internal
#' @noRd
interpol_mat <- function(Y, pos, max_scale = 10) {
  start_pos <- min(pos)
  end_pos   <- max(pos)
  bp <- (pos - start_pos) / (end_pos - start_pos)
  gridn <- power_of_two_gridn(length(pos), max_scale)
  Y_new <- t(apply(Y, 1, interpol_ks, bp = bp, gridn = gridn))

  raw_grid <- makegrid(t = bp, y = seq_len(ncol(Y)),
                                   gridn = gridn)$gridt
  outing_grid <- start_pos +
    (end_pos - start_pos) * (raw_grid - min(raw_grid)) /
    (max(raw_grid) - min(raw_grid))

  list(Y = Y_new, outing_grid = outing_grid)
}

#' Remap a matrix of curves to a power-of-two-length regular grid
#'
#' Detects whether `Y` already has a power-of-two number of columns
#' AND a regularly spaced `pos`. If so, returns `Y` unchanged with
#' `outing_grid = pos`. Otherwise interpolates onto a power-of-two
#' regular grid using `interpol_mat`. NA rows are zero-filled before
#' the interpolation and restored to NA on output.
#'
#' @param Y numeric matrix `n x T`.
#' @param pos numeric vector of length `T`, sampling positions.
#' @param verbose logical; when `TRUE`, emits a `message()` if
#'   interpolation is performed.
#' @param max_scale integer, log2 cap on grid length passed through
#'   to `interpol_mat`.
#' @return list with elements `Y` (possibly remapped matrix) and
#'   `outing_grid` (the grid `Y` is now sampled at).
#' @keywords internal
#' @noRd
remap_data <- function(Y, pos, verbose = TRUE, max_scale = 10) {
  na_rows <- which(!complete.cases(Y))
  if (length(na_rows) > 0) {
    Y[is.na(Y)] <- 0
  }

  log2_ncol <- log2(ncol(Y))
  needs_remap <-
    !is_wholenumber(log2_ncol) ||
    sum(duplicated(diff(pos))) != (length(pos) - 2) ||
    log2_ncol > max_scale

  if (needs_remap) {
    inter <- interpol_mat(Y, pos, max_scale = max_scale)
    Y <- inter$Y
    outing_grid <- inter$outing_grid
    if (verbose) {
      warning_message(
        "ncol(Y) is not 2^J or positions are unevenly spaced; interpolated to a regular dyadic grid.",
        style = "hint")
    }
  } else {
    outing_grid <- pos
  }

  if (length(na_rows) > 0) {
    Y[na_rows, ] <- NA
  }
  list(Y = Y, outing_grid = outing_grid)
}

#' Generate the wavelet-coefficient scale-index list
#'
#' Returns a list mapping wavelet scales to the column indices of a
#' `wd` `$D` vector at level of resolution `lev_res`.
#' Element 1 is scale 0, element 2 is scale 1, ..., element
#' `lev_res` is the finest detail scale, and the final element is
#' the position of the smoothing (C) coefficient.
#'
#' @param lev_res integer, log2 of the signal length. Must be at
#'   most 15 (the upper bound the `wavethresh` indexing convention
#'   handles cleanly).
#' @return list of integer vectors, length `lev_res + 1`.
#' @keywords internal
#' @noRd
#' @examples
#' # The example uses `wavethresh` directly; commented out so the
#' # roxygen build does not require it as a Suggests entry.
#' # library(wavethresh)
#' # tem_func <- rnorm(2^8)
#' # twav <- wd(tem_func)
#' # indx_lst <- gen_wavelet_indx(8)
#' # plot(accessD(twav, level = 6),
#' #      twav$D[unlist(indx_lst[6 + 1])])
gen_wavelet_indx <- function(lev_res) {
  if (lev_res > 15) {
    stop("lev_res must be at most 15; use a smaller level of resolution")
  }
  indx_lst <- list()
  indx_lst[[1]] <- 2^lev_res - 1
  for (s in seq_len(lev_res - 1)) {
    indx_lst[[s + 1]] <- 2^lev_res - (2^(s + 1) - 1):(2^s)
  }
  indx_lst[[length(indx_lst) + 1]] <- 2^lev_res
  indx_lst
}

#' Row-wise discrete wavelet transform of a matrix of curves
#'
#' Applies `wd` to each row of `data`, collecting the D
#' (detail) coefficients into a matrix and the C (smooth) coefficient
#' into a vector. NA rows are zero-filled before the transform and
#' restored to NA on output. The number of columns of `data` MUST be
#' a power of two.
#'
#' `max_scale` is passed through to `wd` as its
#' `min.scale` argument (wavethresh's convention names the smallest
#' decomposition scale, which equals our largest decomposition
#' depth). Setting `max_scale` equal to or above `log2(ncol(data))`
#' performs the full decomposition; smaller values truncate the
#' deepest detail levels. The default `10` matches the cap used by
#' `remap_data`, so a `data` matrix produced by `remap_data` is
#' always decomposed to the full available depth.
#'
#' @param data numeric matrix `n x J`, one curve per row, `J` a
#'   power of two.
#' @param filter_number integer, see `filter.select`.
#'   Default 10 (Daubechies least-asymmetric, 10 vanishing moments).
#' @param family character, wavelet family name passed to
#'   `wd`. Default `"DaubLeAsymm"`.
#' @param max_scale integer, log2 cap shared with `remap_data`.
#'   Forwarded to `wd` as `min.scale`. Default `10`.
#' @return list of class `"DWT"` with elements `C` (length-`n`
#'   vector of smoothing coefficients), `D` (`n x (J - 1)` matrix of
#'   detail coefficients), `J` (log2 of input length),
#'   `filter_number`, `family`, and `max_scale`.
#' @keywords internal
#' @noRd
dwt_matrix <- function(data,
                       filter_number = 10,
                       family        = "DaubLeAsymm",
                       max_scale     = 10) {
  na_rows <- which(!complete.cases(data))
  if (length(na_rows) > 0) {
    data[is.na(data)] <- 0
  }
  J <- ncol(data)
  n <- nrow(data)
  D <- matrix(NA_real_, nrow = n, ncol = J - 1)
  C <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    w <- wd(data[i, ],
                        filter.number = filter_number,
                        family        = family,
                        min.scale     = max_scale)
    D[i, ] <- w$D
    C[i]   <- accessC(w, level = 0)
  }
  if (length(na_rows) > 0) {
    D[na_rows, ] <- NA
    C[na_rows]   <- NA
  }
  out <- list(
    C             = C,
    D             = D,
    J             = log2(J),
    filter_number = filter_number,
    family        = family,
    max_scale     = max_scale
  )
  class(out) <- "DWT"
  out
}
