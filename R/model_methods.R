# Per-effect getter and setter accessors on the `mfsusie` model.
#
# These methods are used by the IBSS loop and by user-facing
# inspection helpers. The pattern mirrors `mvsusieR/R/model_methods.R`
# so the susieR generic dispatch finds them once `zzz.R::.onLoad`
# registers them on the susieR namespace.

#' @keywords internal
#' @noRd
get_alpha_l.mfsusie <- function(model, l) {
  model$alpha[l, ]
}

#' @keywords internal
#' @noRd
get_posterior_moments_l.mfsusie <- function(model, l) {
  list(
    mu  = model$mu[[l]],
    mu2 = model$mu2[[l]]
  )
}

#' Per-effect posterior mean curves (alpha-weighted mu)
#'
#' Returns a list of length M where element m is a `J x T_padded[m]`
#' matrix; row j of element m is `alpha[l, j] * mu[[l]][[m]][j, ]`,
#' the contribution of SNP j to the posterior mean of effect l on
#' modality m.
#' @keywords internal
#' @noRd
get_posterior_mean_l.mfsusie <- function(model, l) {
  alpha_l <- model$alpha[l, ]
  lapply(model$mu[[l]], function(mu_lm) sweep(mu_lm, 1, alpha_l, "*"))
}

#' Posterior mean curves summed over all L effects
#'
#' Returns a list of length M where element m is a `J x T_padded[m]`
#' matrix giving `sum_l alpha[l, j] * mu[[l]][[m]][j, ]`. Used by
#' `get_posterior_mean_sum` callers in the IBSS dispatch.
#' @keywords internal
#' @noRd
get_posterior_mean_sum.mfsusie <- function(model) {
  L <- model$L
  M <- model$M
  out <- vector("list", M)
  for (m in seq_len(M)) {
    out[[m]] <- matrix(0, nrow = nrow(model$mu[[1]][[m]]),
                          ncol = ncol(model$mu[[1]][[m]]))
  }
  for (l in seq_len(L)) {
    pm_l <- get_posterior_mean_l.mfsusie(model, l)
    for (m in seq_len(M)) {
      out[[m]] <- out[[m]] + pm_l[[m]]
    }
  }
  out
}

#' @keywords internal
#' @noRd
get_prior_variance_l.mfsusie <- function(model, l) {
  model$V[l]
}

#' @keywords internal
#' @noRd
set_prior_variance_l.mfsusie <- function(model, l, V) {
  model$V[l] <- V
  model
}
