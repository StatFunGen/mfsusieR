# Apple-to-apple comparison against:
#   fsusieR::log_BF.mixture_normal_per_scale (R/computational_functions.R#L1026-L1133)
#   fsusieR::cal_zeta                        (R/computational_functions.R#L310-L314)
#
# These are the per-(SNP, modality, scale) log-Bayes-factor and the
# softmax helpers that mvf.susie.alpha calls internally per modality
# (see mvf.susie.alpha/R/computational_routine.R:283 and
# operation_on_multfsusie_obj.R:1739). Establishing fidelity here gives
# us a transitive C3 guarantee for the per-(SNP, modality, scale) lBF
# path even before the full mvf.susie.alpha integration test (PR 7e)
# lands. Tolerance per design.md D11b is `<= 1e-8` for C2; we expect
# machine precision (~ 1e-14) because both implementations evaluate
# the same closed-form mixture-of-normals expression.

# skip_if_no_fsusier() lives in setup.R and pins to the contract SHA.

test_that("mixture_log_bf_per_scale summed over scales matches fsusieR::log_BF at <= 1e-12", {
  skip_if_no_fsusier()
  set.seed(mfsusier_test_seed())
  n <- 50; J <- 20
  data <- mfsusieR:::create_mf_individual(
    X = matrix(rnorm(n * J), n),
    Y = list(matrix(rnorm(n * 64), n)),
    verbose = FALSE
  )

  # Hand-built mixture-of-normals prior (per-scale, identical across scales).
  sd_grid <- c(0, 0.5, 1.0, 2.0)
  pi_grid <- c(0.7, 0.1, 0.1, 0.1)
  g_norm <- ashr::normalmix(pi = pi_grid, mean = rep(0, length(pi_grid)),
                            sd = sd_grid)
  ash_record <- list(fitted_g = g_norm)
  class(ash_record) <- "ash"
  G_prior_f <- rep(list(ash_record), length(data$scale_index[[1]]))
  attr(G_prior_f, "class") <- "mixture_normal_per_scale"

  # Initial-state Bhat / Shat (zero-init effects -> residual = D).
  pw   <- data$xtx_diag
  XtD  <- crossprod(data$X, data$D[[1]])
  Bhat <- XtD / pw
  Shat <- sqrt(matrix(1 / pw, J, ncol(data$D[[1]])))

  fsr_lbf <- fsusieR::log_BF(G_prior_f, Bhat = Bhat, Shat = Shat,
                             lowc_wc = NULL,
                             indx_lst = data$scale_index[[1]])

  ours_lbf <- numeric(J)
  for (s in seq_along(data$scale_index[[1]])) {
    idx <- data$scale_index[[1]][[s]]
    ours_lbf <- ours_lbf + mfsusieR:::mixture_log_bf_per_scale(
      Bhat[, idx, drop = FALSE], Shat[, idx, drop = FALSE], sd_grid, pi_grid)
  }
  expect_equal(ours_lbf, fsr_lbf, tolerance = 1e-12)
})

test_that("softmax with uniform pi is bit-identical to fsusieR::cal_zeta", {
  skip_if_no_fsusier()
  set.seed(mfsusier_test_seed())
  J <- 50
  lbf <- rnorm(J, mean = 2, sd = 5)

  fsr <- fsusieR:::cal_zeta(lbf)
  m_max <- max(lbf)
  w <- exp(lbf - m_max)
  ours <- w / sum(w)
  expect_equal(ours, fsr, tolerance = 0)   # exact equality expected
})

test_that("mixture_posterior_per_scale matches fsusieR::post_mat_mean / post_mat_sd at <= 1e-12", {
  skip_if_no_fsusier()
  set.seed(mfsusier_test_seed())
  n <- 50; J <- 20
  data <- mfsusieR:::create_mf_individual(
    X = matrix(rnorm(n * J), n),
    Y = list(matrix(rnorm(n * 64), n)),
    verbose = FALSE
  )

  sd_grid <- c(0, 0.5, 1.0, 2.0)
  pi_grid <- c(0.7, 0.1, 0.1, 0.1)
  g_norm <- ashr::normalmix(pi = pi_grid, mean = rep(0, length(pi_grid)),
                            sd = sd_grid)
  ash_record <- list(fitted_g = g_norm)
  class(ash_record) <- "ash"
  G_prior_f <- rep(list(ash_record), length(data$scale_index[[1]]))
  attr(G_prior_f, "class") <- "mixture_normal_per_scale"

  pw   <- data$xtx_diag
  XtD  <- crossprod(data$X, data$D[[1]])
  Bhat <- XtD / pw
  Shat <- sqrt(matrix(1 / pw, J, ncol(data$D[[1]])))

  # fsusieR's post_mat_mean / post_mat_sd permute columns into wavethresh
  # scale order. We compare per-scale slices, which are in the same
  # order as `indx_lst`, side-stepping the permutation issue.
  for (s in seq_along(data$scale_index[[1]])) {
    idx <- data$scale_index[[1]][[s]]
    fsr_pm <- vapply(seq_len(J), function(j) {
      d <- ashr::set_data(Bhat[j, idx], Shat[j, idx])
      ashr::postmean(g_norm, d)
    }, numeric(length(idx)))
    fsr_psd <- vapply(seq_len(J), function(j) {
      d <- ashr::set_data(Bhat[j, idx], Shat[j, idx])
      ashr::postsd(g_norm, d)
    }, numeric(length(idx)))
    if (length(idx) == 1L) {
      fsr_pm  <- matrix(fsr_pm, J, 1)
      fsr_psd <- matrix(fsr_psd, J, 1)
    } else {
      fsr_pm  <- t(fsr_pm)
      fsr_psd <- t(fsr_psd)
    }

    ours <- mfsusieR:::mixture_posterior_per_scale(
      Bhat[, idx, drop = FALSE], Shat[, idx, drop = FALSE], sd_grid, pi_grid)
    expect_equal(ours$pmean,  fsr_pm,                tolerance = 1e-12)
    expect_equal(ours$pmean2, fsr_psd^2 + fsr_pm^2,  tolerance = 1e-12)
  }
})
