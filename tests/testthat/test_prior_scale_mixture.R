# Apple-to-apple comparison against:
#   fsusieR/R/operation_on_prior.R#L42-L185 (init_prior.default,
#     mixture_normal and mixture_normal_per_scale branches)
#   fsusieR/R/computational_functions.R#L69-L155 (cal_Bhat_Shat,
#     full-Y case, replaced by susieR::compute_marginal_bhat_shat)
#
# Two paths per spec mf-prior/spec.md:
#   1. user-supplied prior_variance_grid -- no ash fit, weights via
#      distribute_mixture_weights.
#   2. data-driven (NULL grid) -- susieR helper + ashr::ash. C2
#      contract: bit-identical to fsusieR::init_prior.default at
#      tolerance 1e-12.

# ---- Helpers -----------------------------------------------------------

make_data <- function(n = 30, J = 8, T_per_modality = c(64L, 128L)) {
  set.seed(mfsusier_test_seed())
  X <- matrix(rnorm(n * J), nrow = n)
  Y <- lapply(T_per_modality,
              function(T_m) matrix(rnorm(n * T_m), nrow = n))
  mfsusieR:::create_mf_individual(X = X, Y = Y, verbose = FALSE)
}

# ---- distribute_mixture_weights ---------------------------------------

test_that("distribute_mixture_weights places null weight on first slot", {
  w <- mfsusieR:::distribute_mixture_weights(K = 3, null_prior_weight = 2)

  expect_length(w, 4L)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_equal(w[1], 2 / 4, tolerance = 1e-12)        # null weight
  expect_equal(w[-1], rep((1 - 2 / 4) / 3, 3), tolerance = 1e-12)
})

test_that("distribute_mixture_weights with null_prior_weight = 0", {
  w <- mfsusieR:::distribute_mixture_weights(K = 1, null_prior_weight = 0)

  expect_length(w, 2L)
  expect_equal(w, c(0, 1), tolerance = 1e-12)
})

# ---- User-supplied prior_variance_grid path ---------------------------

test_that("user-supplied grid is honored, no ash fit runs", {
  data <- make_data()
  v_grid <- c(0.1, 0.5, 1.0)

  prior <- mfsusieR:::mf_prior_scale_mixture(
    data,
    prior_variance_grid = v_grid,
    prior_variance_scope = "per_scale_modality",
    null_prior_weight = 2
  )

  expect_identical(class(prior), "mf_prior_scale_mixture")
  expect_identical(length(prior$V_grid), data$M)
  for (m in seq_len(data$M)) {
    expect_equal(prior$V_grid[[m]], v_grid, tolerance = 0)
  }
})

test_that("user-supplied grid weights match distribute_mixture_weights", {
  data <- make_data()
  v_grid <- c(0.5, 2.0)
  null_w <- 2

  prior <- mfsusieR:::mf_prior_scale_mixture(
    data,
    prior_variance_grid = v_grid,
    prior_variance_scope = "per_scale_modality",
    null_prior_weight = null_w
  )

  expected_pi <- mfsusieR:::distribute_mixture_weights(
    K = length(v_grid), null_prior_weight = null_w
  )
  for (m in seq_len(data$M)) {
    # Each row of pi[[m]] (one row per scale) equals expected_pi.
    for (s in seq_len(nrow(prior$pi[[m]]))) {
      expect_equal(prior$pi[[m]][s, ], expected_pi, tolerance = 1e-12)
    }
  }
})

test_that("susieR-degeneracy contract C1 inputs produce single-component prior", {
  # Per design.md D11a: prior_variance_grid of length 1, null = 0.
  data <- make_data(T_per_modality = 1L)
  prior <- mfsusieR:::mf_prior_scale_mixture(
    data,
    prior_variance_grid = 0.5,           # length-1 grid
    prior_variance_scope = "per_modality",
    null_prior_weight = 0                 # no null component
  )

  expect_identical(prior$V_grid[[1]], 0.5)
  # pi shape: 1 x 2 (null + 1 non-null component).
  expect_identical(dim(prior$pi[[1]]), c(1L, 2L))
  expect_equal(prior$pi[[1]][1, ], c(0, 1), tolerance = 1e-12)
})

# ---- Data-driven path: C2 fidelity vs fsusieR -------------------------

test_that("data-driven init matches fsusieR::init_prior.default at 1e-12", {
  skip_if_not_installed("fsusieR")
  set.seed(2)
  n <- 50; p <- 10; T_padded <- 64
  X <- matrix(rnorm(n * p), nrow = n)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y_m <- matrix(rnorm(n * T_padded), nrow = n)
  scale_index <- mfsusieR:::gen_wavelet_indx(log2(T_padded))

  ours <- mfsusieR:::init_scale_mixture_prior_default(
    Y_m = Y_m, X = X,
    prior_class = "mixture_normal_per_scale",
    groups      = scale_index
  )
  ref <- fsusieR::init_prior(
    Y = Y_m, X = X, prior = "mixture_normal_per_scale",
    v1 = rep(1, n), indx_lst = scale_index, lowc_wc = NULL,
    control_mixsqp = list(verbose = FALSE), nullweight = 0.7
  )

  expect_equal(length(ours$G_prior), length(ref$G_prior))
  expect_equal(ours$G_prior[[1]]$fitted_g$pi,
               ref$G_prior[[1]]$fitted_g$pi, tolerance = 1e-12)
  expect_equal(ours$G_prior[[1]]$fitted_g$sd,
               ref$G_prior[[1]]$fitted_g$sd, tolerance = 1e-12)
  expect_equal(ours$G_prior[[1]]$fitted_g$mean,
               ref$G_prior[[1]]$fitted_g$mean, tolerance = 1e-12)
  expect_equal(ours$tt$Bhat, ref$tt$Bhat, tolerance = 1e-12)
  expect_equal(ours$tt$Shat, ref$tt$Shat, tolerance = 1e-12)
})

# ---- T_m = 1 unification: same code path, no special case --------------

test_that("T_m = 1 with NULL grid runs the data-driven path without crashing", {
  data <- make_data(T_per_modality = 1L)

  prior <- mfsusieR:::mf_prior_scale_mixture(
    data,
    prior_variance_grid = NULL,
    null_prior_weight = 2
  )

  expect_identical(class(prior), "mf_prior_scale_mixture")
  expect_true(length(prior$V_grid[[1]]) >= 1)
})

# ---- Per-modality scope branch ----------------------------------------

test_that("per_modality scope collapses scale dimension", {
  data <- make_data()
  v_grid <- c(0.1, 0.5)

  prior <- mfsusieR:::mf_prior_scale_mixture(
    data,
    prior_variance_grid = v_grid,
    prior_variance_scope = "per_modality",
    null_prior_weight = 2
  )

  for (m in seq_len(data$M)) {
    # pi shape: 1 x (K + 1) under per_modality.
    expect_identical(dim(prior$pi[[m]]), c(1L, length(v_grid) + 1L))
  }
  expect_identical(prior$prior_variance_scope, "per_modality")
})

# ---- Cross-modality independent prior (combine_modality_lbfs) ---------

test_that("cross_modality_prior_independent has the right class hierarchy", {
  prior <- mfsusieR:::cross_modality_prior_independent()
  expect_identical(class(prior),
                   c("mf_prior_cross_modality_independent",
                     "mf_prior_cross_modality"))
})

test_that("combine_modality_lbfs sums per-modality lbfs (independent default)", {
  prior <- mfsusieR:::cross_modality_prior_independent()
  modality_lbfs <- list(c(1, 2, 3), c(0.1, 0.2, 0.3), c(-1, -2, -3))

  out <- mfsusieR:::combine_modality_lbfs(prior, modality_lbfs,
                                           model_state = NULL)

  expect_equal(out, c(0.1, 0.2, 0.3), tolerance = 1e-12)
})

test_that("combine_modality_lbfs with single modality returns the input", {
  prior <- mfsusieR:::cross_modality_prior_independent()
  out <- mfsusieR:::combine_modality_lbfs(prior, list(c(1, 2, 3)),
                                           model_state = NULL)
  expect_equal(out, c(1, 2, 3), tolerance = 0)
})
