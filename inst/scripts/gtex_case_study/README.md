# GTEx-style case-study scripts

This directory ships the case-study scripts originally written
for the GTEx cis-eQTL analyses in the fSuSiE manuscript: AHCYL1,
SCD, HSP90AA1, STMN2, and ABHD17A. Each script reproduces the
genome-browser-overlay figure for one gene by combining an
fSuSiE fit with `Gviz` tracks (ideogram, genome axis, effect
curve with credible band, allele-stratified observed coverage,
gene-region track from `TxDb.Hsapiens.UCSC.hg38.knownGene`).

The scripts are kept here as user-facing recipes — each one is a
self-contained worked example showing how to layer `fsusie()`
output onto a genome browser. They depend on heavyweight
Bioconductor packages that mfsusieR lists as `Suggests` only;
install them with:

    BiocManager::install(c(
      "Gviz", "GenomicFeatures", "AnnotationDbi",
      "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
      "smashr"
    ))

The lighter base-R version of the case study (PIP panel +
per-CS effect curves with credible bands) is the
`vignette("fsusie_gtex_case_study")` worked example, which uses
the bundled `gtex_example` fixture and the package's own
`mfsusie_plot()`.
