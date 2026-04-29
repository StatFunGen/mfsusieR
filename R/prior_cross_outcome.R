# Cross-outcome prior plug-in seam.
#
# Per IBSS effect, the SER step needs a single length-p joint
# log-Bayes factor across the M outcomes to derive the posterior
# `alpha`. The per-outcome log-BFs are first computed independently
# (one length-p vector per outcome `m`, with positions summed
# inside each outcome). They are then combined into a single
# length-p joint log-BF.
#
# The combine step is dispatched through the S3 generic
# `combine_outcome_lbfs(prior, outcome_lbfs, model_state)`. The
# class of `prior` selects the rule: outcome independence sums
# the per-outcome log-BFs (the default), and a non-independence
# rule (e.g. a modality-covariance prior) would register its own
# S3 method on the same generic. The verb is named
# `combine_outcome_lbfs` rather than `apply` to avoid colliding
# with `base::apply`.
#
# Currently only the independence combiner is shipped. The seam
# exists as an extension point so future priors that condition on
# the running model state (passed via `model_state`) can be
# plugged in without touching the SER call site.

#' Combine per-outcome log-Bayes factors into a joint log-BF
#'
#' S3 generic. Each IBSS effect's SER step computes one length-`p`
#' log-BF per outcome (per-position log-BFs already summed inside
#' each outcome); this generic reduces those `M` vectors to a
#' single length-`p` joint log-BF that drives the per-effect
#' `alpha`. The shipped default is outcome independence; a custom
#' rule registers its own S3 method on this generic.
#'
#' @param prior cross-outcome prior object whose class selects the
#'   combine rule. The default
#'   `cross_outcome_prior_independent()` dispatches to outcome-
#'   sum semantics.
#' @param outcome_lbfs list of length `M`. Entry `m` is a numeric
#'   vector of length `p`, the per-outcome log-BF at the current
#'   IBSS effect.
#' @param model_state the running model state (`alpha`, `mu`,
#'   `mu2`, `sigma2`, ...). Unused by the independence combiner;
#'   exposed for combiners that condition on the current fit.
#' @return numeric vector of length `p`, the joint log-BF.
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
#' Constructs the default cross-outcome prior used by
#' `mfsusie()`. Under outcome independence, the joint Bayes factor
#' is the product of per-outcome Bayes factors, so its log is the
#' elementwise sum of the per-outcome log-BFs. The returned object
#' carries no parameters; it acts as an S3 dispatch tag.
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
