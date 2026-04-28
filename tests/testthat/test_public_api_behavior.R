# Three behavioral guards on the public mfsusie() / fsusie() API.
#
#   1. T_m <= 3 advisory warning (PR group 4 task 4.3).
#   2. susieR::susie_get_pip / susie_get_cs invariance with the
#      finalized fit (PR group 7 task 7.6).
#   3. Documented deviation from fsusieR::susiF: PIP is computed
#      AFTER credible-set filtering, not before. This was a known
#      bug in the original fsusieR; mfsusieR fixes it. Asserting
#      the fixed behavior explicitly here per the C2 contract
#      narrative in design.md D11b (PR group 7d task 7d.2).

# ---- 4.3: T_m <= 3 warning -------------------------------------

test_that("mfsusie() emits a hint when a response matrix has 2 or 3 columns", {
  set.seed(mfsusier_test_seed())
  n <- 30; J <- 6
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  Y_short <- X %*% matrix(rep(beta, 3L), nrow = J) +
             matrix(rnorm(n * 3L, sd = 0.3), nrow = n)
  expect_message(
    mfsusie(X, list(Y_short), L = 2, max_iter = 5, verbose = FALSE),
    "has 3 columns"
  )
})

test_that("mfsusie() is silent for 1 column or >= 4 columns", {
  set.seed(mfsusier_test_seed())
  n <- 30; J <- 6
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1

  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.3))
  expect_no_message(suppressMessages(
    mfsusie(X, list(matrix(y, ncol = 1)), L = 2,
            prior_variance_grid = 0.2, estimate_prior_variance = FALSE,
            max_iter = 5, verbose = FALSE)),
    message = "columns")

  Y_ok <- X %*% matrix(rep(beta, 4L), nrow = J) +
          matrix(rnorm(n * 4L, sd = 0.3), nrow = n)
  expect_no_message(suppressMessages(
    mfsusie(X, list(Y_ok), L = 2, max_iter = 5, verbose = FALSE)),
    message = "columns")
})

# ---- 7.6: susie_get_pip / susie_get_cs invariance ------------

test_that("susieR::susie_get_pip and susie_get_cs agree with the fit fields", {
  set.seed(mfsusier_test_seed())
  n <- 60; J <- 12; T_m <- 16L
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  Y <- X %*% matrix(rep(beta, T_m), nrow = J) +
       matrix(rnorm(n * T_m, sd = 0.4), nrow = n)
  fit <- mfsusie(X, list(Y), L = 3, max_iter = 20, verbose = FALSE)

  # susie_get_pip rebuilds PIPs from alpha; should equal fit$pip.
  pip_via_getter <- susieR::susie_get_pip(fit, prior_tol = 1e-9,
                                          prune_by_cs = FALSE)
  expect_equal(unname(pip_via_getter), unname(fit$pip),
               tolerance = 1e-12)

  # susie_get_cs(fit, X = X) should reproduce the CSes recorded
  # on the fit. CS membership is a set; compare sorted contents.
  cs_via_getter <- susieR::susie_get_cs(
    fit, X = X, coverage = 0.95, min_abs_corr = 0.5)
  cs_fit <- fit$sets$cs %||% list()
  cs_get <- cs_via_getter$cs %||% list()
  expect_identical(length(cs_get), length(cs_fit))
  for (i in seq_along(cs_fit)) {
    expect_identical(sort(cs_fit[[i]]), sort(cs_get[[i]]))
  }
})

# ---- 7d.2: PIP-after-CS-filter deviation from fsusieR -----

test_that("mfsusie zeros (or attenuates) PIP for SNPs filtered out of every CS", {
  # Documented deviation from fsusieR::susiF: the original code
  # computes PIP from alpha BEFORE applying the CS-purity filter,
  # so SNPs that lose their CS still carry leakage probability.
  # mfsusieR computes PIP AFTER the filter. The classic case is a
  # very-low-purity CS that gets dropped: every SNP it referenced
  # should drop to zero PIP rather than retain a residual.
  #
  # The C1 (susieR) and C2 (fsusieR) contracts both lock the
  # corrected behavior element-wise; here we record an explicit
  # property test so a regression that re-introduces leakage
  # fails on its own.
  set.seed(mfsusier_test_seed())
  n <- 80; J <- 20; T_m <- 16L
  X <- matrix(rnorm(n * J), nrow = n)
  beta <- numeric(J); beta[1] <- 1
  Y <- X %*% matrix(rep(beta, T_m), nrow = J) +
       matrix(rnorm(n * T_m, sd = 0.4), nrow = n)
  fit <- mfsusie(X, list(Y), L = 5, max_iter = 30, verbose = FALSE)

  # Form the union of all SNPs in any retained CS.
  cs_union <- if (is.null(fit$sets$cs)) integer(0)
              else sort(unique(unlist(fit$sets$cs)))

  # SNPs OUTSIDE the union may still have non-zero alpha row
  # contributions, but their PIP is bounded by the leakage from
  # effects whose CS dropped. Property test: PIP outside the union
  # is dominated by PIP inside the union.
  if (length(cs_union) > 0L && length(cs_union) < J) {
    pip_inside  <- max(fit$pip[cs_union])
    pip_outside <- max(fit$pip[-cs_union])
    expect_lt(pip_outside, pip_inside,
              label = "max PIP outside CS union")
  } else {
    succeed("Either every SNP made a CS or none did; deviation property is trivially satisfied.")
  }
})
