# Audit follow-up B-5: smoke test that mfsusieR's CS construction
# (susieR backbone path) and fsusieR's CS construction (cumsum +
# purity filter) agree on the variant set per CS for a controlled
# fixture. The two paths use different algorithms internally;
# this test pins the OBSERVABLE per-CS variant set.

test_that("mfsusieR::fsusie CS variant set matches fsusieR::susiF on a small two-effect fixture", {
  skip_if_no_fsusier()
  set.seed(31L)
  N <- 80L; J <- 30L; T_m <- 64L
  X <- matrix(rnorm(N * J), N, J)
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 6^2))
  beta <- matrix(0, J, T_m)
  beta[5L, ]  <- 1.2 * shape
  beta[20L, ] <- 0.8 * shape
  Y <- X %*% beta + matrix(rnorm(N * T_m, sd = 0.4), N, T_m)
  pos <- seq_len(T_m)

  fit_ours <- suppressWarnings(mfsusieR::fsusie(
    Y = Y, X = X, pos = pos, L = 5L,
    max_iter = 100L, verbose = FALSE,
    wavelet_qnorm = FALSE, wavelet_standardize = FALSE))
  fit_ref  <- suppressWarnings(fsusieR::susiF(
    Y = Y, X = X, pos = pos, L = 5L, maxit = 100L, verbose = FALSE))

  # The variant SETS in each CS should agree (both packages
  # identify the same causal variants on this fixture). Order /
  # naming may differ; compare as sorted sets.
  cs_ours <- lapply(fit_ours$sets$cs, function(x) sort(unname(x)))
  cs_ref  <- lapply(fit_ref$cs,        function(x) sort(unname(x)))

  expect_equal(length(cs_ours), length(cs_ref))
  ours_sorted <- cs_ours[order(vapply(cs_ours, `[`, integer(1L), 1L))]
  ref_sorted  <- cs_ref[order(vapply(cs_ref, `[`, integer(1L), 1L))]
  for (i in seq_along(ours_sorted)) {
    expect_equal(ours_sorted[[i]], ref_sorted[[i]])
  }
})
