# Design: objective coherence (Part B)

## Current per-effect order

From `post_loglik_prior_hook.mf_individual`, with
`inner_cap = max_inner_em_steps + 1` (default 6):

```
E (SER: loglik, moments, KL)   <- in susieR's single_effect_regression
repeat k = 1..inner_cap:
  M (mixsqp / ebnm)            -> new prior g_k
  if k == inner_cap: break     <- returned state: q from g_{k-1}, prior g_k
  E (loglik, moments, KL)
  if |d lbf| < tol: break      <- returned state: q and KL from g_k, but the
                                  M-step just ran -> prior g_k, consistent
```

The early-exit branch breaks after an E-step, so (q, KL) match the prior
used in that E-step, and the M-step that preceded it produced that same
prior. That exit is consistent. Only the cap exit (and
`max_inner_em_steps = 0`) leaves an M-step without a matching E-step.

## Options considered

| Option | What | Exact? | Cost | Notes |
|---|---|---|---|---|
| B0 | Document that `fit$elbo` is `F(q, h)` | n/a | 0 | No numerical change; leaves `V[l]`, `pi_V` and `fitted_g_per_effect` describing a prior the posterior never saw |
| B1 | Keep q fixed; compute `KL(q_l \|\| g_new)` directly | yes | high | `q_lj` is a per-cell mixture of normals and `g_new` a mixture prior, so `KL(q \|\| g_new) = KL(q \|\| g_old) + E_q[log g_old(b) - log g_new(b)]` has no closed form (expected log of a mixture density). It needs quadrature over every cell with `alpha > eps`. Bookkeeping only. |
| **B2** | **End the cap exit with an E-step** | yes | +1 loglik/moments/KL pass per effect, only when the cap is hit | q, KL, lbf, fitted and the prior all describe one state; this is standard EM order (E, M, ..., M, E). Changes alpha/mu numerically. |
| B3 | Stop calling the M-step on the last cycle (E, M, E, ..., E) | yes | saves one M-step | Equivalent to B2 with `inner_cap - 1` M-steps; reported `max_inner_em_steps` semantics shift by one. |

**Recommendation: B2.** It keeps the `max_inner_em_steps` semantics (the
number of M-steps) and costs one extra E-step only on the cap exit,
which after the planned C++ kernels is about 15% of one inner cycle. It
also fixes BUG-13: `update_fitted_values` runs after the hook, so
`fitted` then reflects the final q. A trimmed effect (V < `prior_tol`)
has a null prior at its final E-step, so its posterior mean is already
0, and zeroing `mu` in `trim_null_effects` no longer leaves `fitted`
stale.

## Interaction with the FDR work

`inst/notes/sessions/2026-05-*` track FDR inflation on
`per_scale_normal` and `per_outcome`. B2 changes alpha whenever the cap
is hit, which is most likely for spurious effects whose prior is still
moving. So B2 could change FDR in either direction. Run it on the
existing null-simulation harness before it becomes the default. A
transitional option (for example `final_e_step = TRUE/FALSE`) is
possible if the user wants a side-by-side comparison; whether to add
one is an open question.

## Open questions for the user

1. B2 as the new default, or behind an option first for the FDR
   comparison?
2. Withdraw `perf-skip-elbo-on-pip` (recommended), or keep it open with
   a rewritten premise?
