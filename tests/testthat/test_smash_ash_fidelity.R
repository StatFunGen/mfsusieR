# C2 fidelity: `mf_smash_ash` is a port of `fsusieR:::smash_lw`.
# Both kernels are deterministic; bit-identity at `tolerance = 1e-14`
# is required on every supported input shape.

test_that("mf_smash_ash matches fsusieR::smash_lw on a power-of-2 length, scalar noise", {
  skip_if_no_fsusier()
  set.seed(1L)
  T_pos <- 128L
  x     <- c(rep(0, 32), rep(1, 32), rep(0, 32), rep(-0.5, 32)) +
           rnorm(T_pos, sd = 0.4)
  ours <- mfsusieR:::mf_smash_ash(noisy_signal = x,
                                  noise_level   = 0.4,
                                  n.shifts      = 50L,
                                  filter.number = 10L,
                                  family        = "DaubExPhase")
  ref  <- fsusieR:::smash_lw(noisy_signal  = x,
                             noise_level   = 0.4,
                             n.shifts      = 50L,
                             filter.number = 10L,
                             family        = "DaubExPhase")
  expect_equal(ours$mu.est,     ref$mu.est,     tolerance = 1e-14)
  expect_equal(ours$mu.est.var, ref$mu.est.var, tolerance = 1e-14)
})

test_that("mf_smash_ash matches fsusieR::smash_lw on a non-power-of-2 length (reflective pad)", {
  skip_if_no_fsusier()
  set.seed(2L)
  T_pos <- 100L
  x     <- sin(seq(0, 2 * pi, length.out = T_pos)) +
           rnorm(T_pos, sd = 0.3)
  ours <- mfsusieR:::mf_smash_ash(noisy_signal = x,
                                  noise_level   = 0.3,
                                  n.shifts      = 25L,
                                  filter.number = 10L,
                                  family        = "DaubExPhase")
  ref  <- fsusieR:::smash_lw(noisy_signal  = x,
                             noise_level   = 0.3,
                             n.shifts      = 25L,
                             filter.number = 10L,
                             family        = "DaubExPhase")
  expect_equal(ours$mu.est,     ref$mu.est,     tolerance = 1e-14)
  expect_equal(ours$mu.est.var, ref$mu.est.var, tolerance = 1e-14)
})

test_that("mf_smash_ash matches fsusieR::smash_lw with a per-position noise_level vector", {
  skip_if_no_fsusier()
  set.seed(3L)
  T_pos <- 128L
  x     <- c(rep(0, 64), rep(1, 64)) + rnorm(T_pos, sd = 0.5)
  noise <- 0.2 + 0.6 * (seq_len(T_pos) > T_pos / 2L)  # heteroscedastic
  ours <- mfsusieR:::mf_smash_ash(noisy_signal = x,
                                  noise_level   = noise,
                                  n.shifts      = 30L,
                                  filter.number = 10L,
                                  family        = "DaubExPhase")
  ref  <- fsusieR:::smash_lw(noisy_signal  = x,
                             noise_level   = noise,
                             n.shifts      = 30L,
                             filter.number = 10L,
                             family        = "DaubExPhase")
  expect_equal(ours$mu.est,     ref$mu.est,     tolerance = 1e-14)
  expect_equal(ours$mu.est.var, ref$mu.est.var, tolerance = 1e-14)
})
