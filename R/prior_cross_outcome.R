# Cross-outcome prior plug-in seam.
#
# The dispatch verb is `combine_outcome_lbfs` (S3 generic), not
# `apply` (collides with apply). When non-NULL, the
# cross_outcome_prior on `mfsusie()` modifies the per-outcome
# log-Bayes factors before they are summed into a joint log-BF in
# the SER step. The default is the trivial independent
# prior (sum of log-BFs unchanged). Future changes may ship
# mash-style cross-outcome priors that implement the same S3
# generic.

#' Combine per-outcome log-Bayes factors into a joint log-BF
#'
#' S3 generic. Default implementation sums per-outcome log-BFs
#' (outcome independence). Custom cross-outcome priors register
#' an S3 method for their class to adjust the combination.
#'
#' @param prior an mfsusieR cross-outcome prior object.
#' @param outcome_lbfs list of length M, each entry a numeric
#'   vector of length p (per-effect, per-outcome log-BFs).
#' @param model_state the running model state (alpha, mu, etc.);
#'   passed for forward compatibility with priors that condition
#'   on the current fit.
#' @return numeric vector of length p, the joint log-BF.
#' @keywords internal
#' @noRd
combine_outcome_lbfs <- function(prior, outcome_lbfs, model_state) {
  UseMethod("combine_outcome_lbfs")
}

#' @exportS3Method combine_outcome_lbfs mf_prior_cross_outcome_independent
#' @keywords internal
#' @noRd
combine_outcome_lbfs.mf_prior_cross_outcome_independent <-
  function(prior, outcome_lbfs, model_state) {
    Reduce(`+`, outcome_lbfs)
  }

#' Independent (outcome-product) cross-outcome prior
#'
#' Returns a stub object whose `combine_outcome_lbfs` method sums
#' per-outcome log-BFs (outcome independence). Default in
#' mfsusie().
#'
#' @return list of class `c("mf_prior_cross_outcome_independent",
#'   "mf_prior_cross_outcome")`.
#' @references
#' Manuscript: methods/online_method.tex line 41.
#' @keywords internal
#' @noRd
cross_outcome_prior_independent <- function() {
  obj <- list()
  class(obj) <- c("mf_prior_cross_outcome_independent",
                  "mf_prior_cross_outcome")
  obj
}
