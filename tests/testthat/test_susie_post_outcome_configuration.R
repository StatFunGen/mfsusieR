# Tests for `susieR::susie_post_outcome_configuration()` on mfsusie fits.
#
# Two algorithms exposed by the susieR entry point are exercised:
#
#   1. SuSiEx (combinatorial 2^N enumeration).
#   2. coloc pairwise (H0/H1/H2/H3/H4 ABF).
#
# The SuSiEx branch is cross-checked against a verbatim re-implementation
# of the legacy reference at
# `mvf.susie.alpha::posthoc_multfsusie()` (`R/operation_on_multfsusie_obj.R`,
# core algorithm at lines ~1542-1574). The reference is inlined here so the
# test has no run-time dependency on `mvf.susie.alpha`.
#
# The coloc branch was independently verified bit-for-bit against
# `coloc::coloc.bf_bf` during development; we don't re-run that check here
# (would add a Suggests dependency on coloc for nothing -- the math is a
# verbatim port of `coloc:::combine.abf`).

# -----------------------------------------------------------------------------
# Verbatim re-implementation of the legacy posthoc_multfsusie core, per CS.
#
# Reference algorithm (mvf.susie.alpha::posthoc_multfsusie):
#   logBF_l   <- as.vector(logBF_trait_snp %*% alpha_l)        # length-S
#   logBF_conf <- as.vector(configs %*% logBF_l)               # length-2^S
#   prob_conf  <- exp(logBF_conf - max(logBF_conf))
#   prob_conf  <- prob_conf / sum(prob_conf)
#   posthoc    <- colSums(configs * prob_conf)                 # length-S
# -----------------------------------------------------------------------------
legacy_posthoc_per_cs <- function(alpha_l, logBF_trait_snp, prob_thresh = 0.8) {
  S       <- nrow(logBF_trait_snp)
  configs <- as.matrix(expand.grid(rep(list(c(0L, 1L)), S)))
  colnames(configs) <- paste0("trait", seq_len(S))

  logBF_l    <- as.vector(logBF_trait_snp %*% alpha_l)
  logBF_conf <- as.vector(configs %*% logBF_l)
  maxlog     <- max(logBF_conf)
  prob_conf  <- exp(logBF_conf - maxlog)
  prob_conf  <- prob_conf / sum(prob_conf)
  posthoc    <- colSums(configs * prob_conf)

  list(
    logBF_trait = logBF_l,
    posthoc     = posthoc,
    active      = posthoc >= prob_thresh,
    configs     = configs,
    config_prob = prob_conf
  )
}

# -----------------------------------------------------------------------------
# Synthetic-fixture parity test: feed the same (alpha, logBF_trait_snp) to
# both implementations and check exact numerical agreement.
# -----------------------------------------------------------------------------
test_that("susiex_configurations matches legacy posthoc_multfsusie on a hand-built fixture", {
  set.seed(1L)
  L <- 2L; J <- 12L; M <- 3L

  # Build a synthetic (alpha, lbf) by hand. Each row of alpha sums to 1.
  alpha <- matrix(0, L, J)
  alpha[1, ] <- c(0.4, 0.3, 0.1, 0.05, 0.05, 0.05, 0.03, 0.02, 0, 0, 0, 0)
  alpha[2, ] <- c(0,   0,   0,   0,    0,    0,    0,    0,    0.5, 0.3, 0.15, 0.05)

  # Per-(L, J, M) log BFs sampled from a wide normal.
  lbf_arr <- array(rnorm(L * J * M, mean = 4, sd = 6), c(L, J, M))

  # Synthesize an mfsusie-shaped fit object directly; no IBSS run needed for
  # this test (we're checking the susieR algorithm, not the IBSS).
  fit <- structure(list(
    alpha                = alpha,
    lbf_variable         = matrix(0, L, J),    # unused by the SuSiEx branch
    lbf_variable_outcome = lbf_arr,
    sets                 = list(cs = list(L1 = which(alpha[1, ] > 1e-3),
                                          L2 = which(alpha[2, ] > 1e-3))),
    L = L
  ), class = c("mfsusie", "susie"))

  out <- susieR::susie_post_outcome_configuration(
    fit, by = "outcome",
    method = "susiex",
    cs_only = FALSE)

  expect_length(out$susiex, L)
  expect_s3_class(out, "susie_post_outcome_configuration")
  expect_identical(attr(out, "method"), "susiex")

  # Reference per CS: pull (alpha_l, t(lbf_arr[l, , ])) so that the per-trait
  # rows of `logBF_trait_snp` are length-J SNP vectors per outcome.
  for (l in seq_len(L)) {
    logBF_trait_snp <- t(lbf_arr[l, , ])      # M x J
    ref <- legacy_posthoc_per_cs(alpha[l, ], logBF_trait_snp)

    got <- out$susiex[[l]]
    expect_equal(unname(got$logBF_trait), unname(ref$logBF_trait),
                 tolerance = 1e-12,
                 info = sprintf("logBF_trait, CS %d", l))
    expect_equal(got$config_prob, ref$config_prob, tolerance = 1e-12,
                 info = sprintf("config_prob, CS %d", l))
    # `marginal_prob` in the susieR API is what the legacy reference
    # called `posthoc`: per-trait marginal P(active) summed across the
    # 2^N configuration ensemble.
    expect_equal(unname(got$marginal_prob), unname(ref$posthoc),
                 tolerance = 1e-12,
                 info = sprintf("marginal_prob, CS %d", l))
    # The configs grid is a column-permutation invariant; both should be
    # bit-equal because both build it via expand.grid().
    expect_equal(unname(got$configs), unname(ref$configs),
                 info = sprintf("configs, CS %d", l))
  }
})

# -----------------------------------------------------------------------------
# Real fit: run mfsusie() end-to-end and check that the configuration output
# has the documented shape and that the posthoc marginal sums dominate the
# planted causal outcomes.
# -----------------------------------------------------------------------------
test_that("susie_post_outcome_configuration on a real mfsusie fit returns the documented shape", {
  set.seed(2026L)
  n <- 60L; p <- 10L; M <- 3L
  T_m <- 16L
  X <- matrix(rnorm(n * p), n)

  # Plant a causal at SNP 2 (shared across outcomes) and SNP 7
  # (also shared, but with a different functional shape per outcome).
  beta <- numeric(p); beta[c(2L, 7L)] <- c(1.4, -0.9)
  shape <- exp(-((seq_len(T_m) - T_m / 2)^2) / (2 * 4^2))
  Y <- lapply(seq_len(M), function(m) {
    X %*% (matrix(beta, p, 1) %*% matrix(shape, 1, T_m)) +
      matrix(rnorm(n * T_m, sd = 0.3), n)
  })

  fit <- mfsusie(X, Y, L = 3, max_iter = 30, verbose = FALSE)
  expect_false(is.null(fit$lbf_variable_outcome),
               info = "lbf_variable_outcome should be attached by default.")
  expect_equal(dim(fit$lbf_variable_outcome), c(3L, p, M))

  # SuSiEx run: at least one CS with all-outcomes-causal as the dominant
  # config. With single-method output, only `$susiex` is present.
  out_susiex <- susieR::susie_post_outcome_configuration(
    fit, by = "outcome", method = "susiex")
  expect_named(out_susiex, "susiex")
  expect_null(out_susiex$coloc_pairwise)
  expect_true(length(out_susiex$susiex) >= 1L)
  any_all_causal <- vapply(out_susiex$susiex, function(e) {
    all(e$marginal_prob >= 0.5)
  }, logical(1L))
  expect_true(any(any_all_causal),
              info = "At least one CS should mark all M outcomes as causal under the planted shared signal.")

  # Coloc pairwise run: PP.H4 (shared causal) should dominate for matched
  # (l, l) pairs. With M = 3 outcomes and 2 surviving CSs, choose(3, 2) x
  # (#CS_pairs) rows total.
  out_coloc <- susieR::susie_post_outcome_configuration(
    fit, by = "outcome", method = "coloc_pairwise")
  expect_named(out_coloc, "coloc_pairwise")
  expect_null(out_coloc$susiex)
  expect_s3_class(out_coloc$coloc_pairwise, "data.frame")
  expect_named(out_coloc$coloc_pairwise,
               c("trait1", "trait2", "l1", "l2", "hit1", "hit2",
                 "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"),
               ignore.order = TRUE)
  matched <- subset(out_coloc$coloc_pairwise, l1 == l2)
  expect_true(all(matched$PP.H4 > 0.5),
              info = "Diagonal coloc-pairs should be H4-dominated for the planted shared causal.")
})

# -----------------------------------------------------------------------------
# Opt-out path: when attach_lbf_variable_outcome = FALSE, by = "outcome"
# should error with a clear directive.
# -----------------------------------------------------------------------------
test_that("susie_post_outcome_configuration errors cleanly when lbf_variable_outcome is missing", {
  set.seed(3L)
  n <- 40L; p <- 8L; M <- 2L
  X <- matrix(rnorm(n * p), n)
  Y <- lapply(seq_len(M), function(m) matrix(rnorm(n * 8L), n))

  fit <- mfsusie(X, Y, L = 2, max_iter = 5, verbose = FALSE,
                 attach_lbf_variable_outcome = FALSE)
  expect_null(fit$lbf_variable_outcome)

  # Pass cs_only = FALSE so the sets-check doesn't fire first; the next
  # validation step then raises the lbf_variable_outcome error.
  expect_error(
    susieR::susie_post_outcome_configuration(fit, by = "outcome",
                                              cs_only = FALSE),
    "lbf_variable_outcome"
  )
})
