# Cross-modality prior plug-in seam (mf-prior/spec.md).
#
# The dispatch verb is `combine_modality_lbfs` (S3 generic), not
# `apply` (collides with apply). When non-NULL, the
# cross_modality_prior on `mfsusie()` modifies the per-modality
# log-Bayes factors before they are summed into a joint log-BF in
# the SER step. The default in v1 is the trivial independent
# prior (sum of log-BFs unchanged). Future changes may ship
# mash-style cross-modality priors that implement the same S3
# generic.

#' Combine per-modality log-Bayes factors into a joint log-BF
#'
#' S3 generic. Default implementation sums per-modality log-BFs
#' (modality independence). Custom cross-modality priors register
#' an S3 method for their class to adjust the combination.
#'
#' @param prior an mfsusieR cross-modality prior object.
#' @param modality_lbfs list of length M, each entry a numeric
#'   vector of length p (per-effect, per-modality log-BFs).
#' @param model_state the running model state (alpha, mu, etc.);
#'   passed for forward compatibility with priors that condition
#'   on the current fit.
#' @return numeric vector of length p, the joint log-BF.
#' @keywords internal
#' @noRd
combine_modality_lbfs <- function(prior, modality_lbfs, model_state) {
  UseMethod("combine_modality_lbfs")
}

#' @exportS3Method combine_modality_lbfs mf_prior_cross_modality_independent
#' @keywords internal
#' @noRd
combine_modality_lbfs.mf_prior_cross_modality_independent <-
  function(prior, modality_lbfs, model_state) {
    Reduce(`+`, modality_lbfs)
  }

#' Independent (modality-product) cross-modality prior
#'
#' Returns a stub object whose `combine_modality_lbfs` method sums
#' per-modality log-BFs (modality independence). Default in
#' mfsusie().
#'
#' @return list of class `c("mf_prior_cross_modality_independent",
#'   "mf_prior_cross_modality")`.
#' @references
#' Manuscript: methods/online_method.tex line 41.
#' @keywords internal
#' @noRd
cross_modality_prior_independent <- function() {
  obj <- list()
  class(obj) <- c("mf_prior_cross_modality_independent",
                  "mf_prior_cross_modality")
  obj
}
