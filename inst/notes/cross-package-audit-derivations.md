# Cross-package audit derivations

Math derivations for the two `fix-now` items from the
2026-04-26 cross-package audit that required mathematical
investigation rather than mechanical edits.

## 1. M-fold scaling of `mixsqp_null_penalty` (audit C-1.3)

**Question.** `mvf.susie.alpha/R/EM.R:58-65` scales `nullweight`
by the number of outcomes $M$ before passing it to the per-
outcome mixsqp M-step:
```r
nullweight_scaled <- nullweight * max(1L, K_f + K_u)
```
mfsusieR's `optimize_prior_variance.mf_individual` did not apply
this scaling. Should it?

**Setup.** The variational EM M-step for the SuSiE
mixture-of-normals prior $\pi$ on effect $l$ solves
$$
\hat{\boldsymbol\pi}
\;=\;
\arg\max_{\boldsymbol\pi}\;
  \sum_j \alpha_{lj} \log \sum_k \pi_k L_{jk}
  + \lambda \log \pi_1,
$$
where:

- $\alpha_{lj}$ is the per-effect SNP posterior weight (the row
  weight `zeta` passed to mixsqp).
- $L_{jk}$ is the marginal Bayes factor of variant $j$ under
  mixture component $k$.
- $\pi_1$ is the null component weight.
- $\lambda$ is the null-component prior penalty (the prepended
  weight `mixsqp_null_penalty * idx_size` in
  `mf_em_m_step_per_scale`).

**M-fold concentration of $\alpha_l$.** mfsusieR's joint
posterior over the SNP indicator $\gamma_l$ is computed from
the joint log-Bayes-factor:
$$
\text{lbf}_{lj}
\;=\; \sum_{m=1}^M \sum_{s=1}^{S_m}
  \log\langle \pi_{m,s}, L_{m,s,j\cdot}\rangle,
\qquad
\alpha_{lj} \propto p_j \exp(\text{lbf}_{lj}).
$$
The per-(outcome, scale) contributions are
near-independent across $j$, so $\operatorname{Var}_j(\text{lbf}_{lj})$
scales linearly with $M$. The softmax over $j$ amplifies the
relative weight of the top SNPs by a factor whose log scales
linearly with $M$, so the maximum value
$\alpha_{l,j^*} \approx M \cdot \alpha_{l,j^*}^{(M=1)}$
in the strong-signal regime, and the effective sample size
$n_{\text{eff}}(\alpha_l) = (\sum_j \alpha_{lj}^2)^{-1}$
shrinks by a factor of $M$.

**Effect on the M-step's first-order condition.** Take the
two-component case for concreteness ($\pi = (1 - \pi_1, \pi_1)$,
with $L_{j0} = 1$, $L_{j1} = \mathrm{BF}_j$):
$$
\sum_j \alpha_{lj} \cdot
  \frac{\mathrm{BF}_j - 1}{1 + \pi_1 (\mathrm{BF}_j - 1)}
\;=\; \frac{\lambda}{1 - \pi_1}.
$$

In the concentrated regime, the LHS sum is dominated by the
top-$\alpha$ SNPs $j^*$, whose contribution scales linearly
with $\alpha_{l,j^*}$:
$$
\text{LHS} \;\approx\; \alpha_{l,j^*} \cdot
  \frac{\mathrm{BF}_{j^*} - 1}{1 + \pi_1 (\mathrm{BF}_{j^*} - 1)}
\;=\; \alpha_{l,j^*} \cdot G(\pi_1, \mathrm{BF}_{j^*}),
$$
with $G(\cdot, \cdot)$ an $O(1)$ function (in the units of LHS).

Consistency of regularization across $M$. Define
$\hat\pi_1$ by setting LHS $=$ RHS:
$$
\hat\pi_1 \;=\; G^{-1}\left(\frac{\lambda}{(1 - \pi_1)\alpha_{l,j^*}}\right).
$$
For $\hat\pi_1$ to be invariant in $M$, the ratio
$\lambda / \alpha_{l,j^*}$ must be invariant. Since
$\alpha_{l,j^*}$ scales linearly with $M$, this requires
$$
\boxed{\;\lambda(M) \;=\; M \cdot \lambda(1).\;}
$$

**Conclusion.** mfsusieR's per-(outcome, scale) M-step uses the
same $\alpha_l$ as mvf.susie.alpha's joint M-step (the joint
posterior is invariant to whether the M-step is per-outcome or
joint). The same $M$-fold concentration argument therefore
applies. To preserve regularization-to-data balance across $M$,
the effective `mixsqp_null_penalty` must scale by $M$.

**Code change.** `R/individual_data_methods.R::optimize_prior_variance.mf_individual`
sets
```r
mixsqp_null_penalty <- (params$mixsqp_null_penalty %||% 0.7) *
                      max(1L, data$M)
```
This mirrors `mvf.susie.alpha/R/EM.R:58-65`'s `nullweight_scaled`
applied per-outcome.

**Backwards compatibility.** Existing fits with $M = 1$ are
unchanged (`max(1L, 1L) = 1`). Fits with $M \geq 2$ now apply
M-fold stronger regularization on the mixture null. This is
expected to tighten PIPs at null variants in multi-outcome
fits and may shift the per-effect $\hat\pi$ at a fixed seed;
the ledger entry records this as a deliberate, derivation-
backed correction.

## 2. ELBO `o` term (audit C-4.4)

**Question.** mvf.susie.alpha's `get_objective` adds an entropy
term $o$ to the ELBO that mfsusieR (following susieR) does not.
Which is correct?

**Standard variational ELBO.** For SuSiE with $L$
single-effect components:
$$
\mathrm{ELBO}
\;=\; \mathbb{E}_q[\log p(Y, b)] - \mathbb{E}_q[\log q(b)]
\;=\; \mathbb{E}_q[\log p(Y \mid b)]
  - \sum_l \mathrm{KL}(q(b_l) \| p(b_l)).
$$
$\mathbb{E}_q[\log p(Y \mid b)]$ is the data term (`Eloglik`).

**Per-effect KL decomposition.** Writing
$b_l = \beta_l \gamma_l$ with $q(\gamma_l = j) = \alpha_{lj}$
and $q(\beta_l \mid \gamma_l = j) = \mathcal{N}(\mu_{lj}, \sigma^2_{lj})$:
$$
\mathrm{KL}(q(b_l) \| p(b_l))
\;=\; \mathrm{KL}(q(\gamma_l) \| p(\gamma_l))
  + \sum_j \alpha_{lj} \mathrm{KL}(q(\beta_l \mid j) \| p(\beta_l \mid j)).
$$
The first term is the categorical KL,
$\sum_j \alpha_{lj} \log(p \cdot \alpha_{lj})$
under uniform prior $1/p$.

**susieR / mfsusieR convention.** susieR computes the full
per-effect KL via the lbf-aggregate identity (see Wang et al.
2020 supplement). The `compute_kl.individual` formula:
$$
\mathrm{KL}_l \;=\; -(\text{lbf}_l + L_{\text{null}}) + E_q[\log p(Y \mid b_l)],
$$
where:

- $\text{lbf}_l = \log \sum_j \pi_j \mathrm{BF}_{lj}$ is the
  per-effect log-marginal-likelihood ratio against the null.
- $L_{\text{null}} = \sum_t \log \mathcal{N}(R_l(t); 0, \sigma^2_t)$
  is the null residual log-likelihood.
- $E_q[\log p(Y \mid b_l)]$ is `SER_posterior_e_loglik`.

**Identity.** Direct calculation (using the softmax form
$\alpha_{lj} = \pi_j \mathrm{BF}_{lj} / \langle \pi, \mathrm{BF}\rangle$):
$$
\text{lbf}_l + L_{\text{null}}
\;=\; -\mathrm{KL}(q(\gamma_l) \| p(\gamma_l))
  + E_q[\log p(Y \mid b_l, \gamma_l)] + L_{\text{null}}.
$$
Rearranging:
$$
\mathrm{KL}_l
\;=\; \mathrm{KL}(q(\gamma_l) \| p(\gamma_l))
  + \underbrace{E_q[\log p(Y \mid b_l)] - E_q[\log p(Y \mid b_l, \gamma_l)] - L_{\text{null}}}_{
        =\;\sum_j \alpha_{lj} \mathrm{KL}(q(\beta_l \mid j) \| p(\beta_l \mid j))
      }.
$$
That is, the susieR formula is **exactly** the standard
per-effect KL decomposition. The categorical entropy is in
$\text{lbf}_l$ via the softmax structure; the Gaussian KL
emerges from $E_q[\log p(Y \mid b_l)] - E_q[\log p(Y \mid b_l, \gamma_l)] - L_{\text{null}}$.

So `Eloglik - sum(KL_l)` is the complete ELBO. **No extra
entropy correction is needed.**

**mvf.susie.alpha convention.** The upstream `cal_KL_l`
returns
$$
\mathrm{KL}_l^{(\text{mvf})}
\;=\; -\,\text{loglik\_SFR}(l) - \text{loglik\_SFR\_post}(l),
$$
where:

- $\text{loglik\_SFR}(l) = \text{lbf}_l^{(\text{model})} + L_{\text{null}}$
  is the joint log-marginal at the prior categorical.
- $\text{loglik\_SFR\_post}(l) = E_q[\log p(Y \mid b_l)]$ is the
  SER posterior expected log-likelihood (a typically negative
  number).

Substituting:
$$
\mathrm{KL}_l^{(\text{mvf})}
\;=\; -\text{lbf}_l - L_{\text{null}} - E_q[\log p(Y \mid b_l)].
$$

Compare to susieR / mfsusieR (re-arranging):
$$
\mathrm{KL}_l
\;=\; -\text{lbf}_l - L_{\text{null}} + E_q[\log p(Y \mid b_l)].
$$

The two differ by the **sign** on $E_q[\log p(Y \mid b_l)]$:
$$
\mathrm{KL}_l^{(\text{mvf})} - \mathrm{KL}_l
\;=\; -\,2\,E_q[\log p(Y \mid b_l)].
$$
Since $E_q[\log p(Y \mid b_l)]$ is a (typically large)
negative number for any non-degenerate Gaussian likelihood,
$\mathrm{KL}_l^{(\text{mvf})}$ exceeds the correct
$\mathrm{KL}_l$ by a large positive amount.

**The `o` term.** mvf's `get_objective` adds
$$
o \;=\; \sum_l \sum_j \alpha_{lj} \log \frac{p}{\max(\alpha_{lj}, 10^{-6})}
\;=\; \sum_l \left[\log p \cdot \sum_j \alpha_{lj} - \sum_j \alpha_{lj}\log\alpha_{lj}\right]
\;=\; \sum_l [\log p + H(\alpha_l)],
$$
which is bounded by $L \log p$. This is much smaller than
$2 \sum_l |E_q[\log p(Y \mid b_l)]|$ (which is of order
$L \cdot n \cdot T_m$), so $o$ does **not** offset the sign error.

**Conclusion.** mvf.susie.alpha's ELBO formula is mathematically
incorrect: the sign on `loglik_SFR_post` should be `+`, not
`-`. The added `o` term is a partial workaround that does not
fully correct the discrepancy. mfsusieR's
`get_objective.mf_individual = Eloglik - sum(KL)` follows the
standard susieR derivation and is the correct ELBO.

**Watertight check.** The ELBO of any correctly-implemented
variational EM algorithm is **monotone non-decreasing** across
iterations. A unit test asserts this property on a multi-iteration
fit; if our ELBO is correct, the value at iteration $t+1$ is
$\geq$ value at iteration $t$ at machine precision (modulo a
small slack for FP summation noise).

**Code (no change to `get_objective`).** The mfsusieR formula
is correct. The audit recommendation is purely a ledger entry
and a tightened test (added in
`tests/testthat/test_variance_and_elbo.R`).
