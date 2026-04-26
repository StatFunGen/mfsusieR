# Single-outcome wrapper. Canonicalises (Y, X, pos) into
# (X, list(Y), list(pos)) and forwards to `mfsusie()`. No numerics.

#' Single-outcome functional SuSiE
#'
#' Convenience wrapper around `mfsusie()` for the single-outcome
#' case (`M = 1`). The response `Y` may be scalar (`T = 1`,
#' susieR-degenerate path) or functional (`T > 1`, wavelet-basis
#' path). Internally `fsusie(Y, X, pos, ...)` is equivalent to
#' `mfsusie(X, list(Y), list(pos), ...)`.
#'
#' Multi-outcome-only arguments (e.g. `cross_outcome_prior`) are
#' rejected with an explicit error to keep the wrapper honest about
#' its `M = 1` scope.
#'
#' @param Y numeric matrix `n x T` of functional responses on a
#'   regular position grid (or numeric vector for `T = 1`).
#' @param X numeric matrix `n x p` of covariates.
#' @param pos optional numeric vector of length `T` recording
#'   sampling positions. Defaults to `seq_len(ncol(Y))`.
#' @param ... forwarded to `mfsusie()`. See `?mfsusie` for the full
#'   parameter set; arguments that only make sense in the
#'   multi-outcome case (currently just `cross_outcome_prior`)
#'   error out.
#'
#' @return a list of class `c("mfsusie", "susie")`. See `?mfsusie`
#'   for the documented field set.
#' @export
fsusie <- function(Y, X, pos = NULL, ...) {
  args <- list(...)

  # Reject arguments that only make sense for M >= 2.
  forbidden_mv <- c("cross_outcome_prior")
  bad <- intersect(names(args), forbidden_mv)
  if (length(bad) > 0L) {
    stop(sprintf(
      "`fsusie()` is the single-outcome wrapper. The following arguments are only meaningful for multi-outcome fits and SHALL not be passed here: %s. Use `mfsusie()` directly.",
      paste(bad, collapse = ", ")))
  }

  if (!is.matrix(Y)) {
    if (is.numeric(Y)) Y <- matrix(Y, ncol = 1)
    else stop("`Y` must be a numeric matrix or numeric vector.")
  }

  pos_arg <- if (is.null(pos)) NULL else list(pos)
  do.call(mfsusie, c(list(X = X, Y = list(Y), pos = pos_arg), args))
}
