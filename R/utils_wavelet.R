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

#' Index of low-count wavelet-coefficient columns
#'
#' Returns the column indices `t` of `Y_wd` where
#' `median(|Y_wd[, t]|) <= threshold`. With the default
#' `threshold = 0`, only columns whose absolute-value median is
#' exactly zero are flagged; raising the threshold widens the
#' set of columns considered uninformative.
#'
#' @param Y_wd numeric matrix; the wavelet-domain response
#'   `cbind(W$D, W$C)`.
#' @param threshold non-negative numeric.
#' @return integer vector of column indices; `integer(0)` if no
#'   column meets the threshold.
#' @references
#' Manuscript: methods/online_method.tex
#' (low-count masking for sparse-coverage wavelet columns).
#' @importFrom stats median
#' @keywords internal
#' @noRd
mf_low_count_indices <- function(Y_wd, threshold = 0) {
  if (threshold < 0) {
    stop("`low_count_filter` must be non-negative.")
  }
  col_med <- apply(abs(Y_wd), 2L, median)
  which(col_med <= threshold)
}

#' Column-wise rank-based normal quantile transform
#'
#' Applies `qnorm(rank(., ties.method = "random") / (n + 1))`
#' column-by-column, after seeding the RNG to `set.seed(1)` for
#' tie-break reproducibility. Each input column is mapped to the
#' standard normal scale via its empirical rank.
#'
#' @param Y_wd numeric matrix.
#' @return numeric matrix; same dimensions as `Y_wd`.
#' @references
#' Manuscript: methods/online_method.tex
#' (column-wise rank-INT for non-Gaussian wavelet coefficients).
#' @importFrom stats qqnorm
#' @keywords internal
#' @noRd
mf_quantile_normalize <- function(Y_wd) {
  if (is.null(dim(Y_wd))) {
    set.seed(1L)
    return(qqnorm(rank(Y_wd, ties.method = "random"),
                  plot.it = FALSE)$x)
  }
  out <- matrix(0, nrow(Y_wd), ncol(Y_wd))
  for (j in seq_len(ncol(Y_wd))) {
    set.seed(1L)
    out[, j] <- qqnorm(rank(Y_wd[, j], ties.method = "random"),
                       plot.it = FALSE)$x
  }
  out
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

# =====================================================================
# Wavelet-variance helpers for translation-invariant smoothing.
#
# Three internal routines used by the TI post-smoother to obtain
# pointwise credible bands. Specialised to station-type + periodic
# + real filters (the only configuration TI uses).
# =====================================================================

#' Stationary wavelet decomposition with squared filter coefficients
#'
#' Returns an object with the same shape as `wd(...,
#' type = "station")` but holding per-coefficient variances rather
#' than means. For a linear filter `W` and input `x` with diagonal
#' covariance `diag(s^2)`, `Var((Wx)_j) = sum_i W_{ji}^2 s_i^2`,
#' so applying the squared filter coefficients to the input
#' yields the per-coefficient variance.
#'
#' @references
#' Manuscript: methods/online_method.tex (TIWT section, Nason 2008
#'   Chapter 2 reference for stationary wavelet variance).
#' @keywords internal
#' @noRd
wd_variance <- function(data, filter.number = 10, family = "DaubLeAsymm") {
  if (is.complex(data))
    stop("`data` must be real; complex input is not supported.")
  data_length <- length(data)
  nlevels <- nlevelsWT(data)
  if (is.na(nlevels))
    stop("Data length is not a power of two")

  # Route through the complex C path with both filter slots set
  # to `H^2`. Imaginary input is zero so the C / D variance
  # outputs read off the real part.
  base_filter <- filter.select(filter.number, family)
  H_sq <- base_filter$H^2
  filt <- list(name = base_filter$name, family = family,
               filter.number = filter.number,
               H = H_sq, G = H_sq)

  fl_dbase <- first.last(LengthH = length(filt$H),
                                     DataLength = data_length,
                                     type = "station", bc = "periodic")
  C <- numeric(fl_dbase$ntotal); C[seq_len(data_length)] <- data
  zero_C <- numeric(fl_dbase$ntotal)
  zero_D <- numeric(fl_dbase$ntotal.d)

  decomp <- .C("comwd",
               CR = as.double(C), CI = as.double(zero_C),
               LengthC = as.integer(fl_dbase$ntotal),
               DR = as.double(zero_D), DI = as.double(zero_D),
               LengthD = as.integer(fl_dbase$ntotal.d),
               HR = as.double(filt$H),
               HI = as.double(numeric(length(filt$H))),
               GR = as.double(filt$G),
               GI = as.double(numeric(length(filt$G))),
               LengthH = as.integer(length(filt$H)),
               nlevels = as.integer(nlevels),
               firstC = as.integer(fl_dbase$first.last.c[, 1]),
               lastC  = as.integer(fl_dbase$first.last.c[, 2]),
               offsetC = as.integer(fl_dbase$first.last.c[, 3]),
               firstD = as.integer(fl_dbase$first.last.d[, 1]),
               lastD  = as.integer(fl_dbase$first.last.d[, 2]),
               offsetD = as.integer(fl_dbase$first.last.d[, 3]),
               ntype = 2L,    # station
               nbc   = 1L,    # periodic
               error = 0L,
               PACKAGE = "wavethresh")
  if (decomp$error != 0L)
    stop("comwd returned error code ", decomp$error)

  out <- list(
    C = complex(real = decomp$CR, imaginary = decomp$CI),
    D = complex(real = decomp$DR, imaginary = decomp$DI),
    nlevels = nlevelsWT(decomp),
    fl.dbase = fl_dbase, filter = filt,
    type = "station", bc = "periodic", date = date()
  )
  class(out) <- "wd"
  out
}

#' Convert a station-type variance `wd` to a `wst`
#'
#' Mirrors the rearrangement that `convert()` does on
#' a station-type wd, but on a variance object.
#' @keywords internal
#' @noRd
wst_variance <- function(wd) {
  if (wd$type != "station")
    stop("Object to convert must be of type 'station'")
  nlevels <- nlevelsWT(wd)
  n <- 2L^nlevels
  # Allocate the wst structure directly. `wst(numeric(n), ...)`
  # would also work but runs a full SWT on zeros; the wp / Carray
  # entries get overwritten by the loop below in either case.
  tmp <- structure(list(
    wp      = matrix(0, nrow = nlevels + 1L, ncol = n),
    Carray  = matrix(0, nrow = nlevels + 1L, ncol = n),
    nlevels = nlevels,
    filter  = wd$filter,
    date    = wd$date
  ), class = "wst")
  arrvec <- getarrvec(nlevelsWT(wd),
                                  sort = FALSE)
  for (lev in (nlevelsWT(wd) - 1L):1L) {
    ds <- accessD.wd(wd, level = lev)
    cs <- accessC.wd(wd, level = lev)
    ord <- arrvec[, nlevelsWT(wd) - lev]
    tmp <- putD(tmp, level = lev, v = ds[ord])
    tmp <- putC(tmp, level = lev, v = cs[ord])
  }
  top_C <- accessC(wd, level = wd$nlevels)
  tmp <- putC(tmp, level = nlevelsWT(wd), v = top_C)
  tmp <- putD(tmp, level = nlevelsWT(wd), v = top_C)
  tmp <- putC(tmp, level = 0L,
                          v = accessC(wd, level = 0L))
  inv_arrvec <- sort.list(levarr(seq_len(n),
    levstodo = nlevelsWT(wd)))
  tmp <- putD(tmp, level = 0L,
    v = accessD(wd, level = 0L)[inv_arrvec])
  tmp
}

#' Variance basis-average (squared-filter `AvBasis`) of a wst
#'
#' Specialised real-valued path; calls wavethresh's compiled
#' `av_basisWRAP` with squared filter coefficients. Output is the
#' pointwise variance of `AvBasis(wst)`.
#' @keywords internal
#' @noRd
av_basis_variance <- function(wst) {
  nlevels <- nlevelsWT(wst)
  H <- wst$filter$H
  G <- wst$filter$G
  n <- 2L^nlevels
  zero <- numeric(length(H))
  out <- .C("comAB_WRAP",
            wstR  = as.double(Re(wst$wp)),
            wstI  = as.double(Im(wst$wp)),
            wstCR = as.double(Re(wst$Carray)),
            wstCI = as.double(Im(wst$Carray)),
            LengthData = as.integer(n),
            level  = as.integer(nlevels - 1L),
            HR = as.double(H), HI = as.double(zero),
            GR = as.double(G), GI = as.double(zero),
            LengthH = as.integer(length(H)),
            answerR = as.double(numeric(n)),
            answerI = as.double(numeric(n)),
            error   = 0L,
            PACKAGE = "wavethresh")
  if (out$error != 0L)
    stop("comAB_WRAP returned error code ", out$error)
  out$answerR
}
