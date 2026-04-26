# Scale-mixture-of-normals prior class for mfsusieR.
#
# Two paths per spec mf-prior/spec.md "per-modality init via susieR
# helper + ash fit":
#
#   1. User supplies `prior_variance_grid` -- use it directly,
#      distribute mixture weights with `null_prior_weight`. No ash
#      fit. The susieR-degeneracy contract C1 takes this path.
#   2. `prior_variance_grid = NULL` -- per modality, call
#      `compute_marginal_bhat_shat(X, data$D[[m]])` and
#      `ash` to fit the K-vector grid, then assemble.
#
# Same control-flow for any `T_m`. NO `T_m == 1` branch. Public
# `mfsusie()` warns when any `T_m <= 3`.
#
# Manuscript reference: methods/derivation.tex eq:additive_model.

#' Distribute mixture weights with a null component
#'
#' Given a length-K vector of non-null variances `V_grid` and a
#' null-component weight `null_prior_weight`, returns the (K+1)-
#' length pi vector with `pi[1] = null_prior_weight / (K + 1)` on
#' the null component (V = 0) and the remaining weight split
#' uniformly across the K non-null components.
#'
#' @param K integer, number of non-null components.
#' @param null_prior_weight numeric, weight on the null component.
#' @return numeric vector of length K + 1, sums to 1.
#' @keywords internal
#' @noRd
distribute_mixture_weights <- function(K, null_prior_weight) {
  null_pi <- null_prior_weight / (K + 1)
  c(null_pi, rep((1 - null_pi) / K, K))
}

#' Per-modality data-driven prior init via susieR helper + ash
#'
#' Calls `compute_marginal_bhat_shat(X, Y_m)` to obtain
#' Bhat / Shat, draws a sample, and fits `ash` to obtain
#' the K-vector grid of mixture variances. The seed sequence
#' (`set.seed(1)`) and sample-size caps (5000 for
#' `mixture_normal`, 50000 for `mixture_normal_per_scale`) are
#' set to satisfy the C2 fidelity contract at machine precision;
#' provenance is in `inst/notes/refactor-exceptions.md`.
#'
#' @param Y_m numeric matrix `n x T_padded[m]` of wavelet
#'   coefficients for modality m.
#' @param X numeric matrix `n x p` of column-centred covariates.
#' @param prior_class one of `"mixture_normal"` or
#'   `"mixture_normal_per_scale"`.
#' @param scale_index list of integer vectors from
#'   `gen_wavelet_indx`. Required for `mixture_normal_per_scale`.
#' @param grid_multiplier numeric, forwarded to `ash` as
#'   `gridmult`.
#' @return list with `G_prior` (the ash fit, possibly replicated)
#'   and `tt` (the marginal Bhat / Shat from the susieR helper).
#' @keywords internal
#' @noRd
init_scale_mixture_prior_default <- function(Y_m,
                                             X,
                                             prior_class     = "mixture_normal_per_scale",
                                             groups          = NULL,
                                             grid_multiplier = sqrt(2)) {
  prior_class <- match.arg(prior_class,
                           c("mixture_normal", "mixture_normal_per_scale"))
  if (is.null(groups)) {
    stop("`groups` is required: a list of column-index vectors covered by each prior entry.")
  }

  bs <- compute_marginal_bhat_shat(X, Y_m)

  sample_size <- if (prior_class == "mixture_normal_per_scale") 50000 else 5000
  pool_dim    <- prod(dim(bs$Bhat))
  draw_n      <- min(pool_dim, sample_size)

  set.seed(1)
  betahat <- c(max(abs(bs$Bhat)), sample(bs$Bhat, size = draw_n))
  set.seed(1)
  sdhat <- c(0.01, sample(bs$Shat, size = draw_n))

  t_ash <- ash(betahat, sdhat,
                     mixcompdist = "normal",
                     outputlevel = 0,
                     gridmult    = grid_multiplier)

  K <- length(t_ash$fitted_g$pi)
  t_ash$fitted_g$pi <- c(0.8, rep(0.2 / (K - 1), K - 1))

  G_prior <- lapply(groups, function(idx) {
    rec <- t_ash
    rec$idx <- idx
    rec
  })
  attr(G_prior, "class") <- prior_class

  list(G_prior = G_prior, tt = bs)
}

#' Build an `mf_prior_scale_mixture` object
#'
#' Constructs the per-modality scale-mixture-of-normals prior.
#' When `prior_variance_grid` is supplied, uses the user-given
#' grid and distributes weights via `null_prior_weight`. When
#' `NULL`, runs the data-driven path
#' (`init_scale_mixture_prior_default`) per modality. Same code
#' path for any `T_m`.
#'
#' @param data an `mf_individual` object.
#' @param prior_variance_grid optional length-K vector of mixture
#'   variances (`sigma_k^2`). When supplied, the data-driven ash
#'   fit is skipped and the grid is used directly. Used by the
#'   susieR-degeneracy contract C1 with `length(grid) == 1`.
#' @param prior_variance_scope `"per_scale_modality"` (default,
#'   stores prior per scale per modality) or `"per_modality"`
#'   (collapses the scale dimension; legacy mode).
#' @param null_prior_weight scalar, default 2 per design.md D5/D7.
#' @param grid_multiplier numeric, forwarded to `ash`.
#' @return list of class `"mf_prior_scale_mixture"`.
#' @references
#' Manuscript: methods/derivation.tex eq:additive_model.
#' @keywords internal
#' @noRd
mf_prior_scale_mixture <- function(data,
                                   prior_variance_grid = NULL,
                                   prior_variance_scope = c("per_scale_modality",
                                                            "per_modality"),
                                   null_prior_weight    = 2,
                                   grid_multiplier      = sqrt(2)) {
  prior_variance_scope <- match.arg(prior_variance_scope)
  if (!inherits(data, "mf_individual")) {
    stop("`data` must be an mf_individual object.")
  }

  M        <- data$M
  T_padded <- data$T_padded
  X        <- data$X

  prior_class <- if (prior_variance_scope == "per_scale_modality") {
    "mixture_normal_per_scale"
  } else {
    "mixture_normal"
  }

  G_prior_per_modality <- vector("list", M)
  V_grid               <- vector("list", M)
  pi_weights           <- vector("list", M)

  use_user_grid <- !is.null(prior_variance_grid)

  for (m in seq_len(M)) {
    # Column-index groups covered by each G_prior entry. The IBSS
    # SER kernels (loglik, calculate_posterior_moments,
    # optimize_prior_variance) iterate `G_prior[[m]]` directly and
    # read `G_prior[[m]][[s]]$idx` to slice (Bhat, Shat) -- so
    # `per_modality` (1 entry covering every wavelet column) and
    # `per_scale_modality` (S_m entries, one per wavelet scale)
    # share a single uniform loop body. No mode-specific branching
    # downstream.
    groups_m <- if (prior_variance_scope == "per_scale_modality") {
      data$scale_index[[m]]
    } else {
      list(unlist(data$scale_index[[m]], use.names = FALSE))
    }

    if (use_user_grid) {
      # User-supplied path: same grid replicated per modality (and
      # per scale when scope = per_scale_modality). No ash fit; the
      # G_prior slot holds a minimal ash-shaped record so downstream
      # mixsqp / EM updates have a uniform interface.
      V_grid[[m]] <- prior_variance_grid
      K           <- length(prior_variance_grid)
      pi_kvec     <- distribute_mixture_weights(K, null_prior_weight)

      sd_grid <- c(0, sqrt(prior_variance_grid))    # null + non-null sds
      G_prior_per_modality[[m]] <- lapply(groups_m, function(idx) {
        list(
          fitted_g = list(pi = pi_kvec, sd = sd_grid, mean = rep(0, K + 1)),
          idx      = idx
        )
      })
      attr(G_prior_per_modality[[m]], "class") <- prior_class
      pi_weights[[m]] <- matrix(pi_kvec, nrow = length(groups_m),
                                ncol = K + 1, byrow = TRUE)
    } else {
      # Data-driven path: susieR helper -> ash. Helper
      # returns one G_prior entry per group, with `$idx` attached.
      out <- init_scale_mixture_prior_default(
        Y_m             = data$D[[m]],
        X               = X,
        prior_class     = prior_class,
        groups          = groups_m,
        grid_multiplier = grid_multiplier
      )
      G_prior_per_modality[[m]] <- out$G_prior

      sd_grid <- out$G_prior[[1]]$fitted_g$sd
      V_grid[[m]] <- sd_grid^2
      K_total <- length(sd_grid)
      pi_kvec <- out$G_prior[[1]]$fitted_g$pi
      pi_weights[[m]] <- matrix(pi_kvec, nrow = length(groups_m),
                                ncol = K_total, byrow = TRUE)
    }
  }

  obj <- list(
    G_prior              = G_prior_per_modality,
    V_grid               = V_grid,
    pi                   = pi_weights,
    null_prior_weight    = null_prior_weight,
    prior_variance_scope = prior_variance_scope,
    update_method        = "mixsqp"
  )
  class(obj) <- "mf_prior_scale_mixture"
  obj
}
