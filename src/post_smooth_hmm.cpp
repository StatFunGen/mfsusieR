// HMM forward / backward / xi kernels for `mf_fit_hmm`. cpp11
// implementations; pure-R oracles in `R/reference_implementations.R`
// (`hmm_forward_R`, `hmm_backward_R`, `hmm_xi_R`) match these at
// <= 1e-12 in `tests/testthat/test_post_smooth_hmm_oracle.R`.
//
// Element-wise / dense-array arithmetic only; no BLAS calls. Per
// design.md D14 (Phase 3 cpp11 stance). The K-by-K transition
// product `alpha %*% P` is small enough (K typically <= 40) that
// a hand-rolled accumulation is faster than calling BLAS.
//
// Loop order convention: `t` outer (sequential dependency on
// alpha[t-1, ] / beta[t+1, ]), state `k` inner. The emission
// matrix is laid out T_pos by K, column-major so column-k access
// is contiguous; we materialize per-row vectors of `emit[t, ]`
// inline.
//
// The forward kernel takes a `t1_normalize` flag. When `false`
// (used in the pre-EM warm-up forward pass) `alpha_hat[1, ]` is
// `pi * emit[1, ]` un-normalized and `G_t[1]` is set to NA_real_,
// matching the upstream `fsusieR::fit_hmm` convention. When `true`
// (used in every EM iteration) `alpha_hat[1, ]` is normalized
// and `G_t[1] = sum(pi * emit[1, ])`.

#include <cpp11.hpp>
#include <cmath>
#include <limits>
#include <vector>

using namespace cpp11;

namespace {

// Multiplies a length-K row vector by a K-by-K matrix, returning
// the length-K result. Used to advance alpha through the
// transition: `alpha_next_unscaled = alpha %*% P`.
inline void row_times_matrix(const std::vector<double>& row,
                             doubles_matrix<>& P,
                             int K,
                             std::vector<double>& out) {
  for (int j = 0; j < K; ++j) {
    double acc = 0.0;
    for (int i = 0; i < K; ++i) acc += row[i] * P(i, j);
    out[j] = acc;
  }
}

}  // namespace

// Forward pass on a precomputed emission matrix.
// Inputs:
//   emit: T_pos by K matrix of emission densities.
//   P:    K by K transition matrix (row-stochastic).
//   pi:   length-K initial state distribution.
//   t1_normalize: if true, normalize alpha_hat[1, ] and set
//                 G_t[1] = sum(pi * emit[1, ]); if false, leave
//                 alpha_hat[1, ] unnormalized and G_t[1] = NA.
// Returns: list(alpha_hat = T_pos x K, G_t = length T_pos).
[[cpp11::register]]
list hmm_forward_cpp(doubles_matrix<> emit,
                     doubles_matrix<> P,
                     doubles pi,
                     bool t1_normalize) {
  const int T_pos = emit.nrow();
  const int K     = emit.ncol();
  if (P.nrow() != K || P.ncol() != K) stop("`P` must be K-by-K.");
  if (pi.size() != K)                  stop("`pi` must have length K.");

  writable::doubles_matrix<> alpha_hat(T_pos, K);
  writable::doubles G_t(T_pos);

  std::vector<double> cur(K), nxt(K);

  // t = 1.
  if (t1_normalize) {
    double s = 0.0;
    for (int k = 0; k < K; ++k) {
      cur[k] = pi[k] * emit(0, k);
      s += cur[k];
    }
    G_t[0] = s;
    for (int k = 0; k < K; ++k) {
      cur[k]            = cur[k] / s;
      alpha_hat(0, k)   = cur[k];
    }
  } else {
    for (int k = 0; k < K; ++k) {
      cur[k]          = pi[k] * emit(0, k);
      alpha_hat(0, k) = cur[k];
    }
    G_t[0] = NA_REAL;
  }

  // t = 2..T_pos.
  for (int t = 1; t < T_pos; ++t) {
    row_times_matrix(cur, P, K, nxt);
    double s = 0.0;
    for (int k = 0; k < K; ++k) {
      nxt[k] *= emit(t, k);
      s += nxt[k];
    }
    G_t[t] = s;
    for (int k = 0; k < K; ++k) {
      cur[k]            = nxt[k] / s;
      alpha_hat(t, k)   = cur[k];
    }
  }

  return writable::list({
    "alpha_hat"_nm = alpha_hat,
    "G_t"_nm       = G_t
  });
}

// Backward pass on a precomputed emission matrix.
// Inputs:
//   emit: T_pos by K matrix of emission densities.
//   P:    K by K transition matrix (row-stochastic).
// Returns: list(beta_hat = T_pos x K, C_t = length T_pos).
//   beta_hat[T_pos, ] = 1; C_t[T_pos] = NA_real_.
//   For t < T_pos: beta_tilde[t, k] = sum_j P[k, j] * beta_hat[t+1, j] * emit[t+1, j]
//                  C_t[t] = max(beta_tilde[t, ])
//                  beta_hat[t, ] = beta_tilde[t, ] / C_t[t]
[[cpp11::register]]
list hmm_backward_cpp(doubles_matrix<> emit,
                      doubles_matrix<> P) {
  const int T_pos = emit.nrow();
  const int K     = emit.ncol();
  if (P.nrow() != K || P.ncol() != K) stop("`P` must be K-by-K.");

  writable::doubles_matrix<> beta_hat(T_pos, K);
  writable::doubles C_t(T_pos);

  std::vector<double> cur(K), prev(K);

  // t = T_pos.
  for (int k = 0; k < K; ++k) {
    cur[k]                 = 1.0;
    beta_hat(T_pos - 1, k) = 1.0;
  }
  C_t[T_pos - 1] = NA_REAL;

  // t = T_pos - 1 downto 1.
  for (int t = T_pos - 2; t >= 0; --t) {
    // weight[j] = beta_hat[t+1, j] * emit[t+1, j]
    std::vector<double> weight(K);
    for (int j = 0; j < K; ++j) weight[j] = cur[j] * emit(t + 1, j);
    // beta_tilde[t, k] = sum_j P[k, j] * weight[j]
    double cmax = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < K; ++k) {
      double acc = 0.0;
      for (int j = 0; j < K; ++j) acc += P(k, j) * weight[j];
      prev[k] = acc;
      if (acc > cmax) cmax = acc;
    }
    C_t[t] = cmax;
    for (int k = 0; k < K; ++k) {
      cur[k]          = prev[k] / cmax;
      beta_hat(t, k)  = cur[k];
    }
  }

  return writable::list({
    "beta_hat"_nm = beta_hat,
    "C_t"_nm      = C_t
  });
}

// xi accumulator: sum_t [outer(alpha[t, ], beta[t+1, ] * emit[t+1, ]) * P] / row_sum_at_t.
// Returns a K-by-K matrix.
[[cpp11::register]]
doubles_matrix<> hmm_xi_cpp(doubles_matrix<> alpha_hat,
                            doubles_matrix<> beta_hat,
                            doubles_matrix<> emit,
                            doubles_matrix<> P) {
  const int T_pos = alpha_hat.nrow();
  const int K     = alpha_hat.ncol();
  if (beta_hat.nrow() != T_pos || beta_hat.ncol() != K)
    stop("`beta_hat` must be T_pos-by-K.");
  if (emit.nrow() != T_pos || emit.ncol() != K)
    stop("`emit` must be T_pos-by-K.");
  if (P.nrow() != K || P.ncol() != K)
    stop("`P` must be K-by-K.");

  writable::doubles_matrix<> xi(K, K);
  for (int i = 0; i < K; ++i)
    for (int j = 0; j < K; ++j) xi(i, j) = 0.0;

  std::vector<double> weight(K);
  std::vector<std::vector<double>> xi_t(K, std::vector<double>(K, 0.0));

  for (int t = 0; t < T_pos - 1; ++t) {
    for (int j = 0; j < K; ++j) weight[j] = beta_hat(t + 1, j) * emit(t + 1, j);
    double s = 0.0;
    for (int i = 0; i < K; ++i) {
      const double a_ti = alpha_hat(t, i);
      for (int j = 0; j < K; ++j) {
        const double v = a_ti * weight[j] * P(i, j);
        xi_t[i][j] = v;
        s += v;
      }
    }
    if (s > 0.0) {
      const double inv_s = 1.0 / s;
      for (int i = 0; i < K; ++i)
        for (int j = 0; j < K; ++j) xi(i, j) += xi_t[i][j] * inv_s;
    }
  }

  return xi;
}
