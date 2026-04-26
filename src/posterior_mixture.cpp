// Per-(variable, scale) mixture-of-normals kernels for the mfsusieR per-effect
// SER step. cpp11 implementation; pure-R oracles in
// `R/reference_implementations.R` (`mixture_log_bf_per_scale_R`,
// `mixture_posterior_per_scale_R`) match these at <= 1e-12 in unit tests.
//
// Element-wise / dense-array arithmetic only; no BLAS calls. Per
// design.md D14 (Phase 3 cpp11 stance). Wholesale matrix-algebra
// acceleration is Phase 7 territory (RcppArmadillo).
//
// Loop order is `i` (column) outer, `j` (row) inner so memory access
// on cpp11::doubles_matrix<> (column-major) is stride-1 contiguous.
//
// Manuscript references:
//   methods/derivation.tex eq:post_f_mix
//   methods/derivation.tex eq:post_f2_mix

#include <cpp11.hpp>
#include <cmath>
#include <vector>

using namespace cpp11;

namespace {

constexpr double kLog2Pi = 1.8378770664093454835606594728112352797;

// Log-density of N(x; 0, var); var > 0 enforced by caller.
inline double log_normal_zero_mean(double x, double var) {
  return -0.5 * (kLog2Pi + std::log(var) + x * x / var);
}

}  // namespace

// Per-(variable, scale) log-Bayes factor for the mixture-of-normals prior.
// Returns a length-J vector. See `mixture_log_bf_per_scale_R` for the
// reference R implementation.
[[cpp11::register]]
doubles mixture_log_bf_per_scale_cpp(doubles_matrix<> bhat,
                                     doubles_matrix<> shat,
                                     doubles sd_grid,
                                     doubles pi_grid,
                                     double v_scale) {
  const int J = bhat.nrow();
  const int Ti = bhat.ncol();
  const int K = sd_grid.size();

  if (shat.nrow() != J || shat.ncol() != Ti) {
    stop("`shat_slice` must have the same shape as `bhat_slice`.");
  }
  if (pi_grid.size() != K) {
    stop("`sd_grid` and `pi_grid` must have equal length.");
  }

  // Precompute log(pi_k) and per-k variance V*sd_k^2.
  std::vector<double> log_pi(K), var_k(K);
  for (int k = 0; k < K; ++k) {
    log_pi[k] = std::log(pi_grid[k]);
    var_k[k]  = sd_grid[k] * sd_grid[k] * v_scale;
  }

  // Initialize output to zero.
  writable::doubles out(J);
  for (int j = 0; j < J; ++j) out[j] = 0.0;

  std::vector<double> log_lr_k(K);

  // i (column) outer, j (row) inner: stride-1 access to bhat / shat.
  for (int i = 0; i < Ti; ++i) {
    for (int j = 0; j < J; ++j) {
      const double b = bhat(j, i);
      const double s2 = shat(j, i) * shat(j, i);
      const double log_dens_null = log_normal_zero_mean(b, s2);

      double m = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < K; ++k) {
        if (var_k[k] == 0.0) {
          log_lr_k[k] = log_pi[k];           // null component: log-LR = 0
        } else {
          const double var_alt = var_k[k] + s2;
          log_lr_k[k] = log_pi[k]
                       + log_normal_zero_mean(b, var_alt)
                       - log_dens_null;
        }
        if (log_lr_k[k] > m) m = log_lr_k[k];
      }
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += std::exp(log_lr_k[k] - m);
      out[j] += m + std::log(s);
    }
  }
  return out;
}

// Per-(variable, position) posterior mean and second moment under the
// per-scale mixture-of-normals prior. Returns a list with `pmean` and
// `pmean2`. See `mixture_posterior_per_scale_R` for the R reference.
[[cpp11::register]]
list mixture_posterior_per_scale_cpp(doubles_matrix<> bhat,
                                     doubles_matrix<> shat,
                                     doubles sd_grid,
                                     doubles pi_grid,
                                     double v_scale) {
  const int J = bhat.nrow();
  const int Ti = bhat.ncol();
  const int K = sd_grid.size();

  if (shat.nrow() != J || shat.ncol() != Ti) {
    stop("`shat_slice` must have the same shape as `bhat_slice`.");
  }
  if (pi_grid.size() != K) {
    stop("`sd_grid` and `pi_grid` must have equal length.");
  }

  std::vector<double> log_pi(K), var_k(K);
  for (int k = 0; k < K; ++k) {
    log_pi[k] = std::log(pi_grid[k]);
    var_k[k]  = sd_grid[k] * sd_grid[k] * v_scale;
  }

  writable::doubles_matrix<> pmean(J, Ti);
  writable::doubles_matrix<> pmean2(J, Ti);

  std::vector<double> log_w(K);
  std::vector<double> shrink_k(K);
  std::vector<double> w(K);

  for (int i = 0; i < Ti; ++i) {
    for (int j = 0; j < J; ++j) {
      const double b = bhat(j, i);
      const double s2 = shat(j, i) * shat(j, i);
      const double log_dens_null = log_normal_zero_mean(b, s2);

      double m_max = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < K; ++k) {
        if (var_k[k] == 0.0) {
          log_w[k]    = log_pi[k];
          shrink_k[k] = 0.0;
        } else {
          const double var_alt = var_k[k] + s2;
          log_w[k] = log_pi[k]
                     + log_normal_zero_mean(b, var_alt)
                     - log_dens_null;
          shrink_k[k] = var_k[k] / var_alt;
        }
        if (log_w[k] > m_max) m_max = log_w[k];
      }
      double e_sum = 0.0;
      for (int k = 0; k < K; ++k) {
        w[k] = std::exp(log_w[k] - m_max);
        e_sum += w[k];
      }
      double pm = 0.0, pm2 = 0.0;
      for (int k = 0; k < K; ++k) {
        const double wk = w[k] / e_sum;
        const double post_mean_k = shrink_k[k] * b;
        const double post_var_k  = shrink_k[k] * s2;
        pm  += wk * post_mean_k;
        pm2 += wk * (post_var_k + post_mean_k * post_mean_k);
      }
      pmean(j, i)  = pm;
      pmean2(j, i) = pm2;
    }
  }

  return writable::list({
    "pmean"_nm  = pmean,
    "pmean2"_nm = pmean2
  });
}
