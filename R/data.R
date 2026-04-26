# Documentation for the bundled simulated example datasets. Build
# script: `data-raw/make_data.R`. The fixtures use the
# `susieR::N3finemapping` LD scaffold for realistic genotype
# correlation and add explicit known causal effects so vignettes
# can demonstrate fine-mapping, post-smoothing, and the
# `mfsusie_plot()` overlays without depending on any external
# data download.

#' Simulated DNA methylation fine-mapping fixture
#'
#' An `fsusie()`-shaped fixture for a cis-mQTL fine-mapping
#' problem: a single functional response on irregular CpG
#' positions over a ~10-kb region, two causal variants with
#' overlapping but oppositely-signed Gaussian methylation
#' effects.
#'
#' @format A list with components
#' \describe{
#'   \item{`X`}{`n x p` genotype matrix (`p = 200`) sliced from
#'     `susieR::N3finemapping$X` so the LD pattern is real.}
#'   \item{`Y`}{`n x T` methylation matrix on irregular CpG
#'     positions (`T = 64`).}
#'   \item{`pos`}{length-`T` integer vector of CpG positions
#'     (base pairs).}
#'   \item{`causal_snps`}{integer vector of column indices in
#'     `X` flagged as truly causal.}
#'   \item{`causal_betas`}{`length(causal_snps) x T` matrix of
#'     the true per-position effects.}
#'   \item{`description`}{free-text description of the fixture.}
#' }
#'
#' @source Simulated. See `data-raw/make_data.R`.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(dnam_example)
#' fit <- fsusie(dnam_example$Y, dnam_example$X,
#'               pos = dnam_example$pos, L = 3, verbose = FALSE)
#' fit_s <- mf_post_smooth(fit, method = "TI",
#'                         wavelet_filter = 1L,
#'                         wavelet_family = "DaubExPhase")
#' mfsusie_plot(fit_s)
#' }
"dnam_example"

#' Simulated RNA-seq fine-mapping fixture
#'
#' An `fsusie()`-shaped fixture for a cis-eQTL fine-mapping
#' problem: a single functional response on evenly-spaced
#' exon-body positions, one causal variant with a peaked
#' Gaussian effect.
#'
#' @format A list with components
#' \describe{
#'   \item{`X`}{`n x p` genotype matrix (`p = 200`) sliced from
#'     `susieR::N3finemapping$X`.}
#'   \item{`Y`}{`n x T` expression coverage matrix (`T = 128`).}
#'   \item{`pos`}{length-`T` integer vector of exon-body
#'     positions.}
#'   \item{`causal_snps`}{integer vector of column indices in
#'     `X` flagged as truly causal.}
#'   \item{`causal_betas`}{matrix with `length(causal_snps)`
#'     rows of true per-position effects.}
#'   \item{`description`}{free-text description.}
#' }
#'
#' @source Simulated. See `data-raw/make_data.R`.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(rnaseq_example)
#' fit <- fsusie(rnaseq_example$Y, rnaseq_example$X,
#'               pos = rnaseq_example$pos, L = 3, verbose = FALSE)
#' fit_s <- mf_post_smooth(fit, method = "HMM")
#' mfsusie_plot(fit_s)
#' }
"rnaseq_example"

#' Simulated multi-omic fine-mapping fixture
#'
#' An `mfsusie()`-shaped fixture for joint fine-mapping across
#' DNA methylation (functional, T = 64), RNA-seq (functional,
#' T = 32), and two scalar QTLs (eQTL, pQTL). Two causal SNPs
#' shared across all four outcomes with per-outcome shapes and
#' signs.
#'
#' @format A list with components
#' \describe{
#'   \item{`X`}{`n x p` genotype matrix (`p = 150`) sliced from
#'     `susieR::N3finemapping$X`.}
#'   \item{`Y_list`}{named list of four outcomes:
#'     `dnam` (n x 64), `rna` (n x 32), `eqtl` (n x 1),
#'     `pqtl` (n x 1).}
#'   \item{`pos_list`}{named list of position vectors (CpG bp,
#'     exon-body indices, scalar dummy `1L` for the QTLs).}
#'   \item{`causal_snps`}{integer vector of shared causal SNPs.}
#'   \item{`description`}{free-text description.}
#' }
#'
#' @source Simulated. See `data-raw/make_data.R`.
#' @keywords datasets
#' @examples
#' \donttest{
#' data(multiomic_example)
#' fit <- mfsusie(multiomic_example$X, multiomic_example$Y_list,
#'                pos = multiomic_example$pos_list, L = 3,
#'                verbose = FALSE)
#' fit_s <- mf_post_smooth(fit, method = "TI",
#'                         wavelet_filter = 1L,
#'                         wavelet_family = "DaubExPhase")
#' mfsusie_plot(fit_s)
#' }
"multiomic_example"
