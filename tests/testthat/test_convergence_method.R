# Tests for `mfsusie(convergence_method = ...)`. Both branches
# (`"elbo"` and `"pip"`) ride on susieR's `check_convergence.default`
# (mfsusieR no longer overrides it as of §3e). At convergence the two
# branches should agree on the fit at float precision -- they only
# differ in WHEN the IBSS loop declares done, not WHAT the fit looks
# like once it has done all the same updates.

# Helper: small high-signal fixture that converges in well under
# max_iter for both convergence rules.
.cm_toy_fit <- function(method) {
  set.seed(11L)
  n <- 80L; p <- 12L; M <- 2L; T_m <- 16L
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[c(2L, 7L)] <- c(1.5, -1.0)
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
  Y <- lapply(seq_len(M), function(m) {
    X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
      matrix(rnorm(n * T_m, sd = 0.3), n)
  })
  # `tol = 1e-10` drives both branches to a near-stationary fixed
  # point so the test's bit-equivalence assertion (`tolerance =
  # 1e-8` on alpha / pip / mu / lbf_variable) is meaningful. With a
  # looser IBSS tol each branch can stop at a different iteration
  # within its own convergence criterion, leaving a residual gap
  # ~tol that breaks tight comparison.
  mfsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X, Y, L = 3,
          max_iter = 200, tol = 1e-10,
          convergence_method = method, verbose = FALSE)
}

test_that("convergence_method = 'elbo' and 'pip' produce numerically equivalent fits", {
  fit_e <- .cm_toy_fit("elbo")
  fit_p <- .cm_toy_fit("pip")

  # Both should converge before max_iter on this clean fixture.
  expect_true(fit_e$converged, info = "ELBO branch should converge")
  expect_true(fit_p$converged, info = "PIP branch should converge")

  # Once both have converged, alpha / lbf_variable / sigma2 are bit-
  # equivalent up to float precision: same SER step, same posterior
  # update, same residual update; the only thing the convergence rule
  # changes is at WHICH iteration the loop stops. With the same tol
  # both branches converge on the same final state.
  expect_equal(fit_e$alpha, fit_p$alpha, tolerance = 1e-8,
               info = "alpha should match across convergence_method")
  expect_equal(fit_e$lbf_variable, fit_p$lbf_variable, tolerance = 1e-8,
               info = "lbf_variable should match across convergence_method")
  expect_equal(fit_e$pip, fit_p$pip, tolerance = 1e-8,
               info = "PIP should match across convergence_method")

  # Posterior moments: list[L] of list[M] of p x T_basis[m] matrices.
  # Walk the structure and compare element-wise.
  expect_equal(length(fit_e$mu), length(fit_p$mu))
  for (l in seq_along(fit_e$mu)) {
    for (m in seq_along(fit_e$mu[[l]])) {
      expect_equal(fit_e$mu[[l]][[m]], fit_p$mu[[l]][[m]],
                   tolerance = 1e-5,
                   info = sprintf("mu[[l=%d]][[m=%d]]", l, m))
    }
  }
})

test_that("convergence_method = 'pip' converges via PIP-difference rather than ELBO", {
  fit <- .cm_toy_fit("pip")
  expect_true(fit$converged)
  # Sanity: niter >= 2 (otherwise the PIP-difference check never
  # fires; the fixture is designed so a couple of sweeps are needed).
  expect_gte(fit$niter, 2L)
})

test_that("convergence_method = 'pip' respects pip_stall_window", {
  # Force a stall by setting tol unreachably tight; pip_stall_window
  # = 2 should declare convergence at iter 2 + 2 = 4 even when the
  # PIP keeps wiggling below tol.
  set.seed(11L)
  n <- 80L; p <- 12L; M <- 2L; T_m <- 16L
  X <- matrix(rnorm(n * p), n, p)
  beta <- numeric(p); beta[c(2L, 7L)] <- c(1.5, -1.0)
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
  Y <- lapply(seq_len(M), function(m) {
    X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
      matrix(rnorm(n * T_m, sd = 0.3), n)
  })
  fit <- suppressWarnings(
    mfsusie(wavelet_qnorm = FALSE, wavelet_standardize = FALSE, X, Y, L = 3, max_iter = 200, tol = 1e-30,
            convergence_method = "pip", pip_stall_window = 2L,
            verbose = FALSE)
  )
  # The fit should converge via the stall path (PIP-difference
  # plateaus when the per-effect priors settle), not by hitting tol
  # (unreachable) and not by hitting max_iter.
  expect_true(fit$converged)
  expect_lt(fit$niter, 200L)
})
