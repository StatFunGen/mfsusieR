# Tests for the simulation helpers `mf_simu_ibss_per_level`
# and `mf_simu_ibss_vanilla`.

test_that("mf_simu_ibss_vanilla returns the documented shape on a length-2^L grid", {
  set.seed(1L)
  out <- mf_simu_ibss_vanilla(lev_res = 7L, length_grid = 8L,
                              pi0 = 0.8)
  expect_type(out, "list")
  expect_named(out, c("sim_func", "true_coef", "true_g", "emp_pi0"),
               ignore.order = TRUE)
  expect_length(out$sim_func, 2L^7L)
  expect_true(all(is.finite(out$sim_func)))
  expect_length(out$true_coef, 2L^7L - 1L)
  expect_s3_class(out$true_g, "normalmix")
  expect_true(out$emp_pi0 >= 0 && out$emp_pi0 <= 1)
})

test_that("mf_simu_ibss_vanilla emp_pi0 tracks pi0 on average", {
  # Monte Carlo sanity: with pi0 = 0.5 over 100 reps the mean
  # empirical null fraction should sit near 0.5.
  emps <- replicate(100L, {
    out <- mf_simu_ibss_vanilla(lev_res = 6L, pi0 = 0.5)
    out$emp_pi0
  })
  expect_true(abs(mean(emps) - 0.5) < 0.05)
})

test_that("mf_simu_ibss_vanilla input validation", {
  expect_error(mf_simu_ibss_vanilla(lev_res = 1L), "lev_res")
  expect_error(mf_simu_ibss_vanilla(length_grid = 1L), "length_grid")
  expect_error(mf_simu_ibss_vanilla(pi0 = -0.1), "pi0")
  expect_error(mf_simu_ibss_vanilla(pi0 = 1.5), "pi0")
})

test_that("mf_simu_ibss_per_level shape matches mf_simu_ibss_vanilla shape", {
  set.seed(1L)
  per_level <- mf_simu_ibss_per_level(lev_res = 6L)
  vanilla   <- mf_simu_ibss_vanilla(lev_res = 6L)
  expect_length(per_level$sim_func, 2L^6L)
  expect_length(vanilla$sim_func,   2L^6L)
  expect_length(per_level$true_coef, 2L^6L - 1L)
  expect_length(vanilla$true_coef,   2L^6L - 1L)
})
