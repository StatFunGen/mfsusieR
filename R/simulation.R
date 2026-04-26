# Simulation helpers for example data and vignettes.
#
# `mf_simu_ibss_per_level` samples a smooth random function on a
# length-`2^lev_res` grid by drawing each wavelet scale's detail
# coefficients from a per-scale ash normal mixture and applying
# the inverse DWT. The mixture grid of standard deviations is
# the same across scales; the per-scale mixing weights have a
# null component whose probability follows
# `1 - exp(-prop_decay * scale_index)`. The resulting prior on
# the position-space function is the same one fSuSiE infers in
# its IBSS step, so the bundled simulation matches the prior
# class the model assumes.

#' Sample a smooth random function from the IBSS wavelet prior
#'
#' Draws a length-`2^lev_res` function whose wavelet detail
#' coefficients at scale `s` are independent draws from an
#' ash-style normal mixture
#' `\sum_k \pi_{s,k} N(0, (\sigma_k / 2^{\alpha s})^2)`.
#' The first mixture component is a point mass at zero with
#' weight `pi0[s]`; the remaining components share a common
#' grid of standard deviations on `(0, \infty)`. The inverse
#' DWT (`wavethresh::wr`) maps the sampled coefficient vector
#' to position space.
#'
#' This is the prior class fSuSiE infers in its IBSS step, so
#' the bundled simulation matches the model's assumed
#' generative process. Used by the example datasets and the
#' covariate-adjustment vignette.
#'
#' @param lev_res integer `>= 2`. Resolution level: the
#'   sampled function lives on a length-`2^lev_res` grid.
#'   Default `7` (length 128).
#' @param length_grid integer `>= 2`. Number of components in
#'   the per-scale ash mixture, including the null component
#'   at zero. Default `10`.
#' @param pi0 optional numeric vector of length `lev_res`.
#'   Per-scale null-component weight. When omitted, defaults
#'   to `1 - exp(-prop_decay * 1:lev_res)`.
#' @param alpha non-negative numeric. Scale-dependent variance
#'   shrinkage exponent: at scale `s`, the per-component sd is
#'   divided by `2^{alpha * s}`. Larger `alpha` produces
#'   smoother functions. Default `0.8`.
#' @param prop_decay numeric in `[0, 1]`. Controls how the
#'   default `pi0` decays across scales: larger `prop_decay`
#'   gives a sparser fine-scale signal. Ignored when `pi0` is
#'   supplied explicitly. Default `0.1`.
#'
#' @return A list with components
#' \describe{
#'   \item{`sim_func`}{numeric vector of length `2^lev_res`,
#'     the position-space sample.}
#'   \item{`true_coef`}{numeric vector of the wavelet detail
#'     coefficients used to construct the sample.}
#'   \item{`mix_per_scale`}{length-`lev_res` list of
#'     `ashr::normalmix` objects, one per wavelet scale.}
#'   \item{`emp_pi0`}{length-`lev_res` numeric vector of
#'     empirical null-component fractions per scale (the
#'     fraction of detail coefficients that came back exactly
#'     zero from the mixture sampler).}
#' }
#'
#' @references
#' Manuscript: methods/derivation.tex (IBSS prior derivation).
#'
#' @examples
#' set.seed(1L)
#' f <- mf_simu_ibss_per_level(lev_res = 7L)
#' length(f$sim_func)            # 128
#' plot(f$sim_func, type = "l", xlab = "position", ylab = "f")
#'
#' @importFrom wavethresh wd accessD wr
#' @importFrom ashr normalmix
#' @importFrom stats rchisq runif rnorm
#' @export
mf_simu_ibss_per_level <- function(lev_res     = 7L,
                                   length_grid = 10L,
                                   pi0,
                                   alpha       = 0.8,
                                   prop_decay  = 0.1) {
  if (!is.numeric(lev_res) || length(lev_res) != 1L ||
      !is.finite(lev_res) || lev_res < 2L) {
    stop("`lev_res` must be a single integer >= 2.")
  }
  if (!is.numeric(length_grid) || length(length_grid) != 1L ||
      !is.finite(length_grid) || length_grid < 2L) {
    stop("`length_grid` must be a single integer >= 2.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0) {
    stop("`alpha` must be a non-negative scalar.")
  }
  if (!is.numeric(prop_decay) || length(prop_decay) != 1L ||
      !is.finite(prop_decay) || prop_decay < 0 || prop_decay > 1) {
    stop("`prop_decay` must be a scalar in [0, 1].")
  }
  lev_res     <- as.integer(lev_res)
  length_grid <- as.integer(length_grid)

  if (missing(pi0)) {
    pi0 <- 1 - exp(-prop_decay * seq_len(lev_res))
  } else if (length(pi0) != lev_res ||
             any(!is.finite(pi0)) ||
             any(pi0 < 0) || any(pi0 > 1)) {
    stop("`pi0` must be a length-`lev_res` numeric vector with entries in [0, 1].")
  }

  # Index lists for the wavethresh-style packed wavelet
  # coefficient vector. Scale `s` (1-indexed) occupies
  # positions `(2^(s-1)):(2^s - 1)` in the detail vector.
  indx_lst       <- vector("list", lev_res)
  indx_lst[[1L]] <- 1L
  for (s in seq_len(lev_res - 1L)) {
    indx_lst[[s + 1L]] <- (2L^s):(2L^(s + 1L) - 1L)
  }

  # Shared per-scale grid of mixture sds: a 0-anchored cumsum
  # of chi-squared draws. The first component is the null
  # (sd = 0).
  grid <- c(0, cumsum(rchisq(length_grid - 1L, df = 3)))

  G_level <- vector("list", lev_res)
  for (i in seq_len(lev_res)) {
    tt     <- runif(length_grid - 1L)
    tt     <- (1 - pi0[i]) * tt / sum(tt)
    pi_sim <- c(pi0[i], tt)
    G_level[[i]] <- ashr::normalmix(pi_sim,
                                    rep(0, length_grid),
                                    grid)
  }

  T_grid <- 2L^lev_res
  twav   <- wavethresh::wd(rep(0, T_grid))

  # Resample if every detail coefficient came back zero (the
  # all-null case is a degenerate function and would defeat
  # the purpose of the simulation).
  while (sum(twav$D == 0) == length(twav$D)) {
    for (i in seq_len(lev_res)) {
      for (k in unlist(indx_lst[i])) {
        clust <- sample.int(length_grid, size = 1L,
                            prob = G_level[[i]]$pi)
        twav$D[k] <- rnorm(1, mean = 0,
                           sd = grid[clust] / (2^(alpha * i)))
      }
    }
    # Reverse the detail vector to match wavethresh's coarse-
    # to-fine ordering convention used by `wr()`.
    twav$D <- rev(twav$D)
  }

  sim_func <- wavethresh::wr(twav)
  emp_pi0  <- numeric(lev_res)
  for (i in seq_len(lev_res) - 1L) {
    d_i <- wavethresh::accessD(twav, level = i)
    emp_pi0[i + 1L] <- if (length(d_i) == 0L) 0 else
      sum(d_i == 0) / length(d_i)
  }

  list(sim_func      = sim_func,
       true_coef     = twav$D,
       mix_per_scale = G_level,
       emp_pi0       = emp_pi0)
}
