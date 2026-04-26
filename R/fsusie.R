# Thin wrapper for `mfsusie()` matching `fsusieR::susiF`'s argument
# order. Provides a one-line migration path for users moving off
# fsusieR. The wrapper performs no numerics: it only canonicalizes
# (Y, X, pos) into the (X, list(Y), list(pos)) shape that
# `mfsusie()` expects, then forwards.
#
# Per design.md D8d and `mf-public-api/spec.md`. Argument order
# matches `fsusieR::susiF(Y, X, ...)`.

#' Single-modality functional SuSiE (drop-in for `fsusieR::susiF`)
#'
#' Convenience wrapper around `mfsusie()` for the single-modality
#' functional case (`M = 1`, `T_1 > 1`). Argument order matches
#' `fsusieR::susiF(Y, X, ...)` so existing code can switch with a
#' rename of the function. Internally `fsusie(Y, X, pos, ...)` is
#' equivalent to `mfsusie(X, list(Y), list(pos), ...)`.
#'
#' Multi-modality-only arguments (e.g. `cross_modality_prior`) are
#' rejected with an explicit error to keep the wrapper honest about
#' its `M = 1` scope.
#'
#' @param Y numeric matrix `n x T` of functional responses on a
#'   regular position grid.
#' @param X numeric matrix `n x J` of covariates.
#' @param pos optional numeric vector of length `T` recording
#'   sampling positions. Defaults to `seq_len(ncol(Y))`.
#' @param ... forwarded to `mfsusie()`. See `?mfsusie` for the full
#'   parameter set; arguments that only make sense in the
#'   multi-modality case (currently just `cross_modality_prior`)
#'   error out.
#'
#' @return a list of class `c("mfsusie", "susie")`. See `?mfsusie`
#'   for the documented field set.
#' @export
fsusie <- function(Y, X, pos = NULL, ...) {
  args <- list(...)

  # Reject arguments that only make sense for M >= 2.
  forbidden_mv <- c("cross_modality_prior")
  bad <- intersect(names(args), forbidden_mv)
  if (length(bad) > 0L) {
    stop(sprintf(
      "`fsusie()` is the single-modality wrapper. The following arguments are only meaningful for multi-modality fits and SHALL not be passed here: %s. Use `mfsusie()` directly.",
      paste(bad, collapse = ", ")))
  }

  if (!is.matrix(Y)) {
    if (is.numeric(Y)) Y <- matrix(Y, ncol = 1)
    else stop("`Y` must be a numeric matrix or numeric vector.")
  }

  pos_arg <- if (is.null(pos)) NULL else list(pos)
  do.call(mfsusie, c(list(X = X, Y = list(Y), pos = pos_arg), args))
}
