# Fix `fsusie_covariates_adjustment.Rmd` (faithful port from fsusieR)

## Why

The current `vignettes/fsusie_covariates_adjustment.Rmd` is
incorrect: the prose claims that adjusting for covariates
sharpens PIP at the causal SNP, but the actual numerics show
identical results adjusted vs unadjusted. The cause is the
simulation: the covariate confounding is too weak to matter
(`sd = 1` covariates with a heavily-decayed Haar smoothing
helper), so adjustment is effectively a no-op. The narrative
was AI-written to match what *should* happen on a properly-
calibrated example; on this dataset it does not.

The original fsusieR vignette
`~/GIT/fsusieR/vignettes/Adjusting_covariate.Rmd` simulates
covariates with `sd = 2` and per-covariate effect curves drawn
from `simu_IBSS_per_level(7)$sim_func` (an IBSS-prior wavelet
sample). The resulting covariate signal dominates the genetic
signal at each position, so adjustment matters and the
narrative matches the numerics.

This change ports the original vignette faithfully, replaces
the AI-hallucinated narrative, and adds a final OLS section
demonstrating covariate adjustment for scalar traits.

## What changes

### 1. Port `simu_IBSS_per_level` as a public exported function

Add `mf_simu_ibss_per_level()` to a new file
`R/simulation.R`. The function mirrors
`fsusieR::simu_IBSS_per_level` mathematically (sample a
function from the wavelet IBSS prior on a length-`2^lev_res`
grid by drawing each scale's coefficients from a per-scale
ash mixture) and is exported with full roxygen, including the
parameter contract, the underlying generative model, and a
runnable example. Improvements over the upstream version:

- **Argument validation**: `lev_res >= 2`, `length_grid >= 2`,
  `0 <= alpha`, `0 <= prop_decay <= 1`, `pi0` length
  matches `lev_res` if supplied.
- **Roxygen completeness**: every parameter documented; the
  return value enumerated; the IBSS prior reference cited.
- **No hidden RNG calls**: the function draws RNG only via
  `rnorm`, `runif`, `rchisq`, `sample` — same as upstream;
  the user controls the stream with `set.seed()`.
- **Naming**: `mf_simu_ibss_per_level` follows the package
  `mf_*` convention.

### 2. Vignette rewrite

`vignettes/fsusie_covariates_adjustment.Rmd` is rewritten
faithfully against the upstream vignette layout:

- **Keep**: the introduction paragraph (concise, accurate
  motivation for covariate adjustment in functional fine-
  mapping). This is the only section worth keeping.
- **Replace**: every other section.

New section structure:

1. **Setup** — load `mfsusieR`, `susieR`, `wavethresh`. Set
   seed.
2. **Generating the data** — port the original simulation
   verbatim (modulo the cleanup notes below) using
   `mf_simu_ibss_per_level`. `Cov` has `sd = 2`; per-covariate
   effects f1_cov, f2_cov, f3_cov drawn from the IBSS prior;
   genetic signal has the same prior class.
3. **Account for covariates** — call
   `mf_adjust_for_covariates(Y, Cov, method = "wavelet_eb")`,
   plot the first fitted covariate effect against the truth
   `f1_cov` (overlay). The redundant
   `Y_corrected == Y - Cov %*% fitted_func` sanity check
   from the upstream vignette is dropped; the equality holds
   by construction and adds no information.
4. **Recovering the genetic signal** — two-panel scatter
   `Y` vs `target.data` (observed vs genetic-only) and
   `Y_corrected` vs `target.data` (adjusted vs genetic-only).
   The second scatter should hug the diagonal much more
   tightly.
5. **Fine-map** — fit `fsusie(Y_corrected, X, L = 10)` and
   show `mfsusie_plot(fit)`.
6. **Comparison with no adjustment** — fit
   `fsusie(Y, X, L = 10)` and overlay both fits' recovered
   effect curves against the simulated truth `f1`, `f2`. The
   adjusted fit tracks the truth; the unadjusted fit is
   biased toward the covariate-induced confounding.
7. **OLS for scalar traits** (NEW; after the original
   vignette content). When the response is a scalar (`T = 1`)
   or when smoothness across positions is not relevant,
   `mf_adjust_for_covariates(Y_scalar, Cov, X = X,
   method = "ols")` is the closed-form residualization. Show
   the API and explain that the OLS path returns both
   `Y_adjusted` and `X_adjusted`, which is the FWL
   correction — important when covariates are correlated
   with genotype.

### 3. Logic / readability fixes vs the upstream original

- Rename `target.data` -> `genetic_signal` and
  `noisy.data` -> `covariate_signal` for clarity.
- Drop the unused `svd(N3finemapping$X[,500:700])` line.
- Replace `attach(N3finemapping)` with explicit
  `N3finemapping$X[, 1:100]` etc.
- In the comparison plot, match the recovered effect curve
  to its truth curve by lead-SNP match
  (`fit$sets$cs[[i]]` vs `pos1`/`pos2`) rather than
  hard-coding `l = 1` / `l = 2`.

These are naming + robustness fixes; the math is unchanged.

### 4. Coloc FIXME

`vignettes/fsusie_colocalization.Rmd` carries an inline FIXME
("nothing shows up here!! ... read more carefully the
coloc.susie documentation"). Out of scope for this proposal;
addressed under a separate small change after the main fix.

## Impact

- New: `R/simulation.R` (one exported function), updated
  roxygen man page.
- Changed: `vignettes/fsusie_covariates_adjustment.Rmd`
  rewritten; intro paragraph kept verbatim.
- DESCRIPTION: no changes (uses existing dependencies).
- Specs: new
  `inst/openspec/specs/mf-simulation-helpers/spec.md`
  capability for the simulation helper.

## Out of scope

- Porting `fsusieR::EBmvFR` directly. We use
  `mf_adjust_for_covariates(method = "wavelet_eb")`, which is
  the in-package equivalent that fixes documented upstream
  bugs (per `inst/notes/refactor-exceptions.md`).
- Coloc vignette FIXME. Tracked for a follow-up.
- Other `simu_*` helpers from
  `fsusieR/R/Simulations_functions.R`. Only
  `simu_IBSS_per_level` is needed for this vignette.
