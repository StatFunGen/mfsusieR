# mf-vignette-why-functional capability

## ADDED Requirements

### Requirement: vignette uses the cosine simulation and compares five methods

`vignettes/fsusie_why_functional.Rmd` SHALL use a single
simulation (cosine effect on positions 41..128 of T = 128 with
a single causal SNP, on `n = 100` rows of
`susieR::N3finemapping$X`) and SHALL compare exactly five
methods: per-CpG SuSiE, top-PC SuSiE for PC = 1..5, mvSuSiE
with default prior, mvSuSiE with the true outer-product prior
and fixed identity residual variance, and fSuSiE.

#### Scenario: simulation recipe

The vignette SHALL produce the response as
`y[i, ] = X[i, 700] * effect + rnorm(T)` where
`effect = 1.2 * cos((1:T) / T * 3 * pi)` truncated to zero on
positions 1..40.

#### Scenario: methods compared

For each of the five methods named above, the vignette SHALL
report (i) the credible sets returned, (ii) the PIP at the
true causal SNP, and (iii) for the two methods that recover
the correct CS (mvSuSiE with true prior + fixed identity
residual variance, and fSuSiE), an effect-curve overlay plot
of the recovered effect on the truth.

### Requirement: vignette gates mvsusieR availability

The vignette SHALL gate every code chunk that loads or calls
`mvsusieR` with `eval = requireNamespace("mvsusieR", quietly = TRUE)`.
The package SHALL build on environments without mvsusieR
installed; mvsusieR-dependent chunks SHALL be skipped.

#### Scenario: mvsusieR-absent build

When `requireNamespace("mvsusieR", quietly = TRUE)` returns
FALSE, the vignette SHALL render without errors and SHALL
omit the mvsusieR comparison sections; the per-CpG SuSiE,
top-PC SuSiE, and fSuSiE sections SHALL still render.

### Requirement: no port-source mentions

The vignette SHALL NOT mention `fsusieR`, `mvf.susie.alpha`,
`susiF`, `plot_susiF`, or any port-source package by name.
mvSuSiE is named directly because it is a SuSiE-suite peer,
not a port source.

#### Scenario: vignette grep audit

A grep over `vignettes/fsusie_why_functional.Rmd` for
`fsusieR`, `mvf.susie.alpha`, `susiF`, or `plot_susiF` SHALL
return zero matches.

### Requirement: mvsusieR effect-curve plot uses inline base R

The vignette SHALL plot the mvsusieR effect curve via inline
base R `plot()` and `lines()` calls using `coef(fit_mv)`. It
SHALL NOT extend the public `mfsusie_plot()` API to handle
mvsusie fits, and SHALL NOT add a new public plot helper.

#### Scenario: mvsusie effect plot recipe

The vignette SHALL extract the effect curve via
`coef(fit_mv)[lead_snp + 1, ]` and overlay it on the truth
curve in a single base-R `plot() + lines()` block.

### Requirement: vignette reproduces the upstream `Limitation_SuSiE` result for shared methods

The vignette SHALL use the same simulation recipe and seed
as the upstream `Limitation_SuSiE.Rmd` (cosine effect on
positions 41..128, single causal SNP at position 700,
`n = 100` rows of `susieR::N3finemapping$X`, iid Gaussian
noise). For the methods that overlap with upstream
(per-CpG SuSiE and fSuSiE), the vignette's numerical output
SHALL match upstream's at the same fixed seed: the same
credible sets, the same PIP at the true causal SNP, and
(for fSuSiE) the same recovered effect curve at numerical
tolerance `<= 1e-12`.

#### Scenario: per-CpG SuSiE PIP table matches upstream

For the per-CpG SuSiE section, the vignette SHALL produce
the same per-position PIPs that upstream produces under the
same seed. The number of CSes per position (`length(fit_t$sets$cs)`)
SHALL match upstream exactly.

#### Scenario: fSuSiE recovered curve matches upstream within tolerance

For the fSuSiE section, the recovered effect curve
(post `mf_post_smooth(method = "TI")`) and credible band
SHALL match the upstream `susiF()` output under the same
seed at tolerance `<= 1e-12` on the curve and exact match on
the credible-set membership.
