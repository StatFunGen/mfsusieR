# Bit-identity tests for the wavelet-variance helpers used by the
# TI post-smoother. Tolerance is machine precision (max diff = 0
# exact).

skip_if_no_fsusier_var <- function() {
  testthat::skip_if_not_installed("fsusieR")
}

test_that("wd_variance is bit-identical across (filter.number, family)", {
  skip_if_no_fsusier_var()
  combos <- expand.grid(
    fn  = c(1L, 4L, 10L),
    fam = c("DaubExPhase", "DaubLeAsymm"),
    stringsAsFactors = FALSE
  )
  combos <- combos[!(combos$fn == 1L & combos$fam == "DaubLeAsymm"), ]
  set.seed(7)
  for (i in seq_len(nrow(combos))) {
    fn  <- combos$fn[i]; fam <- combos$fam[i]
    x   <- rnorm(128)
    o   <- mfsusieR:::wd_variance(x, fn, fam)
    r   <- fsusieR:::wd.var(x, fn, fam, type = "station", bc = "periodic")
    expect_equal(Re(o$D), Re(r$D), tolerance = 0,
                 info = sprintf("fn=%d fam=%s D-real", fn, fam))
    expect_equal(Im(o$D), Im(r$D), tolerance = 0,
                 info = sprintf("fn=%d fam=%s D-imag", fn, fam))
    expect_equal(Re(o$C), Re(r$C), tolerance = 0,
                 info = sprintf("fn=%d fam=%s C-real", fn, fam))
  }
})

test_that("wst_variance + av_basis_variance compose bit-identically with upstream", {
  skip_if_no_fsusier_var()
  combos <- expand.grid(
    fn  = c(1L, 4L, 10L),
    fam = c("DaubExPhase", "DaubLeAsymm"),
    n_pos = c(64L, 128L),
    stringsAsFactors = FALSE
  )
  combos <- combos[!(combos$fn == 1L & combos$fam == "DaubLeAsymm"), ]
  set.seed(11)
  for (i in seq_len(nrow(combos))) {
    fn   <- combos$fn[i]; fam <- combos$fam[i]
    npos <- combos$n_pos[i]
    x    <- rnorm(npos)
    o_wd <- mfsusieR:::wd_variance(x, fn, fam)
    r_wd <- fsusieR:::wd.var(x, fn, fam, type = "station", bc = "periodic")
    o_w  <- mfsusieR:::wst_variance(o_wd)
    r_w  <- fsusieR:::convert.var(r_wd)
    expect_equal(Re(o_w$wp), Re(r_w$wp), tolerance = 0)
    expect_equal(Im(o_w$wp), Im(r_w$wp), tolerance = 0)
    expect_equal(Re(o_w$Carray), Re(r_w$Carray), tolerance = 0)
    o_av <- mfsusieR:::av_basis_variance(o_w)
    r_av <- fsusieR:::AvBasis.var(r_w)
    expect_equal(o_av, r_av, tolerance = 0,
                 info = sprintf("fn=%d fam=%s n=%d", fn, fam, npos))
  }
})
