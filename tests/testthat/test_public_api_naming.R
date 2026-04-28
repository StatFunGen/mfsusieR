# PR group 7 task 7.9: forbidden public-API parameter names.
#
# Per `mf-public-api/spec.md` and CLAUDE.md naming rules, the public
# `mfsusie()` argument list MUST NOT use any of the legacy / abbreviated
# names from `mvf.susie.alpha` and `fsusieR`. Every renamed argument
# follows the principle "describe what the parameter IS, not what it
# does", with snake_case throughout.
#
# This test guards against future regressions: if anyone reintroduces
# an abbreviated name in `formalArgs(mfsusie)`, the test fails.

test_that("mfsusie() arguments do not use any forbidden / abbreviated name", {
  forbidden <- c(
    # Legacy abbreviations from mvf.susie.alpha / fsusieR.
    "max_SNP_EM",
    "null_weight",
    "max_scale",
    "min_scale",
    "init_pi0_w",
    "espsilon",
    "max_step",
    "min_purity",
    "low_trait",
    "lowc_wc",
    "df",
    "init_pi0_w",
    # Generic R footguns we forbid in our public API.
    "x", "y",
    "no_residual",
    "no_intercept",
    "no_standardize"
  )
  args <- formalArgs(mfsusie)
  bad  <- intersect(args, forbidden)
  expect_identical(bad, character(0),
                   info = sprintf("Forbidden args found: %s",
                                  paste(bad, collapse = ", ")))
})

test_that("mfsusie() argument names are snake_case (apart from SuSiE lineage)", {
  # SuSiE lineage allows uppercase X / Y (matrices) and L (effect
  # count). Everything else MUST be snake_case per CLAUDE.md naming
  # rules.
  lineage <- c("X", "Y", "L", "L_greedy")
  args <- setdiff(formalArgs(mfsusie), c("...", lineage))
  bad  <- args[!grepl("^[a-z][a-z0-9_]*$", args)]
  expect_identical(bad, character(0),
                   info = sprintf("Non-snake_case args: %s",
                                  paste(bad, collapse = ", ")))
})

test_that("mfsusie() argument list contains the canonical public-API names", {
  required <- c("X", "Y", "pos", "L", "prior_variance_grid",
                "null_prior_weight", "residual_variance_scope",
                "max_iter", "tol",
                "L_greedy", "greedy_lbf_cutoff",
                "estimate_prior_variance",
                "estimate_residual_variance",
                "convergence_method", "pip_stall_window")
  args <- formalArgs(mfsusie)
  missing_args <- setdiff(required, args)
  expect_identical(missing_args, character(0),
                   info = sprintf("Missing canonical args: %s",
                                  paste(missing_args, collapse = ", ")))
})

