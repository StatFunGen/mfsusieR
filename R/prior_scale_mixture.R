# Scale-mixture-of-normals prior class for mfsusieR.
#
# Two paths "per-outcome init via susieR
# helper + ash fit":
#
#   1. User supplies `prior_variance_grid` -- use it directly,
#      distribute mixture weights with `null_prior_init`. No ash
#      fit. Setting `length(grid) == 1` reproduces the
#      single-Gaussian semantics of `susieR::susie()`.
#   2. `prior_variance_grid = NULL` -- per outcome, call
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
#' direct null-component initialization mass `null_prior_init` in
#' `[0, 1]`, returns the (K+1)-length pi vector with
#' `pi[1] = null_prior_init` on the null component (V = 0) and
#' the remaining weight `1 - null_prior_init` split uniformly
#' across the K non-null components.
#'
#' @param K integer, number of non-null components.
#' @param null_prior_init numeric in `[0, 1]`, initial mass on the
#'   null component.
#' @return numeric vector of length K + 1, sums to 1.
#' @keywords internal
#' @noRd
distribute_mixture_weights <- function(K, null_prior_init) {
  null_pi <- null_prior_init
  c(null_pi, rep((1 - null_pi) / K, K))
}

#' Per-outcome data-driven prior init
#'
#' Calls `compute_marginal_bhat_shat(X, Y_m)` to obtain
#' Bhat / Shat, draws a sample, and fits `ash` to obtain
#' the K-vector grid of mixture variances. The internal
#' `set.seed(1)` and sample-size caps (5000 for
#' `mixture_normal`, 50000 for `mixture_normal_per_scale`) make
#' the prior grid deterministic for a given input.
#'
#' @param Y_m numeric matrix `n x T_basis[m]` of wavelet
#'   coefficients for outcome m.
#' @param X numeric matrix `n x p` of column-centred covariates.
#' @param prior_class one of `"mixture_normal"` or
#'   `"mixture_normal_per_scale"`.
#' @param scale_index list of integer vectors from
#'   `gen_wavelet_indx`. Required for `mixture_normal_per_scale`.
#' @param grid_multiplier numeric, forwarded to `ash` as
#'   `gridmult`.
#' @param lowc_idx integer vector of column indices in `Y_m`
#'   masked by `wavelet_magnitude_cutoff`. When non-empty, those
#'   columns are excluded from the ash sampling pool so
#'   masked-zero coefficients do not pull the prior toward a
#'   degenerate spike.
#' @return list with `G_prior` (the ash fit, possibly replicated)
#'   and `tt` (the marginal Bhat / Shat from the susieR helper).
#' @keywords internal
#' @noRd
init_scale_mixture_prior_default <- function(Y_m,
                                             X,
                                             prior_class       = "mixture_normal_per_scale",
                                             groups            = NULL,
                                             grid_multiplier   = sqrt(2),
                                             lowc_idx          = integer(0),
                                             null_prior_init = 0,
                                             na_idx            = NULL) {
  prior_class <- match.arg(prior_class,
                           c("mixture_normal", "mixture_normal_per_scale"))
  if (is.null(groups)) {
    stop("`groups` is required: a list of column-index vectors covered by each prior entry.")
  }

  if (!is.null(na_idx)) {
    Y_m <- Y_m[na_idx, , drop = FALSE]
    X   <- X[na_idx,   , drop = FALSE]
  }
  bs <- compute_marginal_bhat_shat(X, Y_m)

  if (length(lowc_idx) > 0L) {
    pool_Bhat <- bs$Bhat[, -lowc_idx, drop = FALSE]
    pool_Shat <- bs$Shat[, -lowc_idx, drop = FALSE]
  } else {
    pool_Bhat <- bs$Bhat
    pool_Shat <- bs$Shat
  }

  sample_size <- if (prior_class == "mixture_normal_per_scale") 50000 else 5000
  pool_dim    <- prod(dim(pool_Bhat))
  draw_n      <- min(pool_dim, sample_size)

  set.seed(1)
  betahat <- c(max(abs(pool_Bhat)), sample(pool_Bhat, size = draw_n))
  set.seed(1)
  sdhat <- c(0.01, sample(pool_Shat, size = draw_n))

  t_ash <- ash(betahat, sdhat,
                     mixcompdist = "normal",
                     outputlevel = 0,
                     gridmult    = grid_multiplier)

  # Use ash's sd-grid but discard its fitted pi: under LD the iid
  # assumption ash makes is violated. Set the init pi directly
  # from `null_prior_init` (a probability in `[0, 1]`); the EM
  # M-step washes this out within a few iterations, so this is
  # only the cold-start point. The same parameter drives both the
  # ash-driven path and the user-supplied-grid path
  # (`distribute_mixture_weights`).
  K <- length(t_ash$fitted_g$pi)
  pi_null <- null_prior_init
  t_ash$fitted_g$pi <- c(pi_null, rep((1 - pi_null) / (K - 1), K - 1))

  G_prior <- lapply(groups, function(idx) {
    rec <- t_ash
    rec$idx <- idx
    rec
  })
  attr(G_prior, "class") <- prior_class

  list(G_prior = G_prior, tt = bs)
}

#' Per-outcome ebnm-driven prior init for the per-scale point-* paths
#'
#' Mirrors `init_scale_mixture_prior_default()` shape for the
#' two ebnm-backed prior classes. For each group (one wavelet
#' scale `s` with column index set `idx_s`), picks the lead
#' variable from the marginal data via
#' `which.max(rowMeans(Bhat[, idx_s]^2))` and fits a 2-component
#' point-* prior on the lead variable's per-scale slice via the
#' matching `ebnm::ebnm_point_*()` call. The lead picker at init
#' is alpha-independent because `alpha` is uniform `1/p` at IBSS
#' iter 0 and `which.max` would tie-break to variable 1.
#'
#' @param Y_m numeric matrix `n x T_basis[m]` of wavelet
#'   coefficients for outcome m.
#' @param X numeric matrix `n x p` of column-centred covariates.
#' @param prior_class one of `"mixture_point_normal_per_scale"`
#'   or `"mixture_point_laplace_per_scale"`.
#' @param groups list of integer vectors from
#'   `data$scale_index[[m]]`. One group per scale.
#' @param na_idx integer vector of complete-case row indices for
#'   outcome m, or `NULL` for all rows.
#' @return list with `G_prior` (one ebnm fit per scale, `idx`
#'   attached, `lead_init_s` recorded) and `tt` (the marginal
#'   Bhat / Shat from the susieR helper).
#' @keywords internal
#' @importFrom ebnm ebnm_point_normal ebnm_point_laplace
#' @noRd
init_ebnm_prior_per_scale <- function(Y_m,
                                       X,
                                       prior_class,
                                       groups,
                                       na_idx = NULL) {
  prior_class <- match.arg(prior_class,
                           c("mixture_point_normal_per_scale",
                             "mixture_point_laplace_per_scale"))
  if (is.null(groups)) {
    stop("`groups` is required: a list of column-index vectors covered by each prior entry.")
  }
  ebnm_fn <- switch(prior_class,
    mixture_point_normal_per_scale  = ebnm::ebnm_point_normal,
    mixture_point_laplace_per_scale = ebnm::ebnm_point_laplace
  )

  if (!is.null(na_idx)) {
    Y_m <- Y_m[na_idx, , drop = FALSE]
    X   <- X[na_idx,   , drop = FALSE]
  }
  bs <- compute_marginal_bhat_shat(X, Y_m)

  G_prior <- lapply(groups, function(idx) {
    bhat_slice <- bs$Bhat[, idx, drop = FALSE]
    shat_slice <- bs$Shat[, idx, drop = FALSE]
    lead_s     <- which.max(rowMeans(bhat_slice^2))
    fit        <- ebnm_fn(x = bhat_slice[lead_s, ],
                          s = shat_slice[lead_s, ])
    list(
      fitted_g    = fit$fitted_g,
      idx         = idx,
      lead_init_s = lead_s
    )
  })
  attr(G_prior, "class") <- prior_class

  list(G_prior = G_prior, tt = bs)
}

#' Build an `mf_prior_scale_mixture` object
#'
#' Constructs the per-outcome scale-mixture-of-normals prior.
#' When `prior_variance_grid` is supplied, uses the user-given
#' grid and distributes weights via `null_prior_init`. When
#' `NULL`, runs the data-driven path
#' (`init_scale_mixture_prior_default`) per outcome. Same code
#' path for any `T_m`.
#'
#' @param data an `mf_individual` object.
#' @param prior_variance_grid optional length-K vector of mixture
#'   variances (`sigma_k^2`). When supplied, the data-driven ash
#'   fit is skipped and the grid is used directly. Setting
#'   `length(grid) == 1` collapses the mixture to a single
#'   Gaussian (matches `susieR::susie()` semantics).
#' @param prior_variance_scope one of `"per_scale"` (default,
#'   stores prior per scale per outcome), `"per_outcome"`
#'   (collapses the scale dimension), `"per_scale_normal"`
#'   (per-scale point-normal prior fit by `ebnm`), or
#'   `"per_scale_laplace"` (per-scale point-laplace prior fit by
#'   `ebnm`).
#' @param null_prior_init scalar in `[0, 1]`, initial `pi[null]`
#'   for the scale-mixture prior. Default `0` (the EM washes it
#'   out within a few iterations regardless of starting value).
#' @param grid_multiplier numeric, forwarded to `ash`.
#' @return list of class `"mf_prior_scale_mixture"`.
#' @references
#' Manuscript: methods/derivation.tex eq:additive_model.
#' @keywords internal
#' @noRd
mf_prior_scale_mixture <- function(data,
                                   prior_variance_grid = NULL,
                                   prior_variance_scope = c("per_scale",
                                                            "per_outcome",
                                                            "per_scale_normal",
                                                            "per_scale_laplace"),
                                   null_prior_init    = 0,
                                   grid_multiplier      = sqrt(2)) {
  prior_variance_scope <- match.arg(prior_variance_scope)
  if (!inherits(data, "mf_individual")) {
    stop("`data` must be an mf_individual object.")
  }

  M        <- data$M
  T_basis <- data$T_basis
  X        <- data$X

  prior_class <- switch(prior_variance_scope,
    per_outcome       = "mixture_normal",
    per_scale         = "mixture_normal_per_scale",
    per_scale_normal  = "mixture_point_normal_per_scale",
    per_scale_laplace = "mixture_point_laplace_per_scale"
  )
  use_ebnm <- prior_variance_scope %in%
    c("per_scale_normal", "per_scale_laplace")

  G_prior_per_outcome <- vector("list", M)
  V_grid               <- vector("list", M)
  pi_weights           <- vector("list", M)

  use_user_grid <- !is.null(prior_variance_grid)

  for (m in seq_len(M)) {
    # Column-index groups covered by each G_prior entry. The IBSS
    # SER kernels (loglik, calculate_posterior_moments,
    # optimize_prior_variance) iterate `G_prior[[m]]` directly and
    # read `G_prior[[m]][[s]]$idx` to slice (Bhat, Shat) -- so
    # `per_outcome` (1 entry covering every wavelet column) and
    # `per_scale` / `per_scale_*` (S_m entries, one per wavelet
    # scale) share a single uniform loop body. No mode-specific
    # branching downstream.
    groups_m <- if (prior_variance_scope == "per_outcome") {
      list(unlist(data$scale_index[[m]], use.names = FALSE))
    } else {
      data$scale_index[[m]]
    }

    if (use_ebnm) {
      if (use_user_grid) {
        if (length(prior_variance_grid) != 1L) {
          stop("`prior_variance_grid` on `per_scale_normal` / `per_scale_laplace` must have length 1 (a fixed slab variance); got length ",
               length(prior_variance_grid), ".")
        }
        # User-supplied fixed slab path. Builds the matching
        # per-scale `fitted_g` (normalmix for `per_scale_normal`,
        # laplacemix for `per_scale_laplace`) without calling ebnm
        # at init; pairs with `estimate_prior_variance = FALSE` to
        # hold the prior fixed across IBSS (susie-degenerate lock).
        sigma <- sqrt(prior_variance_grid[1L])
        pi_kvec <- c(null_prior_init, 1 - null_prior_init)
        if (prior_class == "mixture_point_normal_per_scale") {
          fg_proto <- structure(
            list(pi = pi_kvec, mean = c(0, 0), sd = c(0, sigma)),
            class = "normalmix", row.names = c(1L, 2L))
        } else {
          # Match the Laplace slab variance to V = sigma^2:
          # Laplace variance = 2 * scale^2, so scale = sigma / sqrt(2).
          fg_proto <- structure(
            list(pi = pi_kvec, mean = c(0, 0),
                 scale = c(0, sigma / sqrt(2))),
            class = "laplacemix", row.names = c(1L, 2L))
        }
        G_prior_per_outcome[[m]] <- lapply(groups_m, function(idx) {
          list(fitted_g = fg_proto, idx = idx, lead_init_s = NA_integer_)
        })
        attr(G_prior_per_outcome[[m]], "class") <- prior_class
        pi_weights[[m]] <- matrix(pi_kvec, nrow = length(groups_m),
                                  ncol = 2L, byrow = TRUE)
        V_grid[[m]] <- rep(prior_variance_grid[1L], length(groups_m))
      } else {
        out <- init_ebnm_prior_per_scale(
          Y_m         = data$D[[m]],
          X           = X,
          prior_class = prior_class,
          groups      = groups_m,
          na_idx      = data$na_idx[[m]]
        )
        G_prior_per_outcome[[m]] <- out$G_prior
        pi_kvec <- out$G_prior[[1L]]$fitted_g$pi
        pi_weights[[m]] <- matrix(pi_kvec, nrow = length(groups_m),
                                  ncol = length(pi_kvec), byrow = TRUE)
        V_grid[[m]] <- vapply(out$G_prior, function(g) {
          sl <- g$fitted_g
          # Slab variance: sd[2]^2 for normalmix, scale[2]^2 * 2 for
          # laplacemix (Laplace variance is 2 * scale^2). Stored for
          # introspection only; not consumed downstream.
          if (!is.null(sl$sd)) sl$sd[2L]^2
          else 2 * sl$scale[2L]^2
        }, numeric(1))
      }
    } else if (use_user_grid) {
      # User-supplied path: same grid replicated per outcome (and
      # per scale when scope = per_scale). No ash fit; the
      # G_prior slot holds a minimal ash-shaped record so downstream
      # mixsqp / EM updates have a uniform interface.
      V_grid[[m]] <- prior_variance_grid
      K           <- length(prior_variance_grid)
      pi_kvec     <- distribute_mixture_weights(K, null_prior_init)

      sd_grid <- c(0, sqrt(prior_variance_grid))    # null + non-null sds
      G_prior_per_outcome[[m]] <- lapply(groups_m, function(idx) {
        list(
          fitted_g = list(pi = pi_kvec, sd = sd_grid, mean = rep(0, K + 1)),
          idx      = idx
        )
      })
      attr(G_prior_per_outcome[[m]], "class") <- prior_class
      pi_weights[[m]] <- matrix(pi_kvec, nrow = length(groups_m),
                                ncol = K + 1, byrow = TRUE)
    } else {
      # Data-driven path: susieR helper -> ash. Helper
      # returns one G_prior entry per group, with `$idx` attached.
      out <- init_scale_mixture_prior_default(
        Y_m               = data$D[[m]],
        X                 = X,
        prior_class       = prior_class,
        groups            = groups_m,
        grid_multiplier   = grid_multiplier,
        lowc_idx          = if (!is.null(data$lowc_idx)) data$lowc_idx[[m]]
                            else integer(0),
        null_prior_init = null_prior_init,
        na_idx            = data$na_idx[[m]]   # complete-case rows for trait m
      )
      G_prior_per_outcome[[m]] <- out$G_prior

      sd_grid <- out$G_prior[[1]]$fitted_g$sd
      V_grid[[m]] <- sd_grid^2
      K_total <- length(sd_grid)
      pi_kvec <- out$G_prior[[1]]$fitted_g$pi
      pi_weights[[m]] <- matrix(pi_kvec, nrow = length(groups_m),
                                ncol = K_total, byrow = TRUE)
    }
  }

  obj <- list(
    G_prior              = G_prior_per_outcome,
    V_grid               = V_grid,
    pi                   = pi_weights,
    null_prior_init    = null_prior_init,
    prior_variance_scope = prior_variance_scope,
    update_method        = if (use_ebnm) "ebnm" else "mixsqp"
  )
  class(obj) <- "mf_prior_scale_mixture"
  obj
}
