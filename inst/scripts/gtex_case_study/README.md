# GTEx-style case-study scripts

This directory ships a self-contained worked example for the
GTEx cis-eQTL genome-browser figures (AHCYL1, SCD, HSP90AA1,
STMN2, ABHD17A) used in published fSuSiE figures. The published
figures use real GTEx individual-level data which cannot be
redistributed with this package, so the script operates on the
bundled simulated `gtex_example` dataset and reproduces the same
multi-track Gviz layout end-to-end.

## Files

- `AHCYL1.R` — the multi-CS Gviz figure on simulated GTEx-style
  data: chromosome axis, fitted effect track with 95% credible
  band, allele-stratified observed difference, and per-genotype
  log1p coverage. Optional TxDb / ideogram block at the bottom
  for users who want to swap in a real `(chrom, locus, gene)`
  triple.

## Bioconductor dependencies

The script uses `Gviz`, `GenomicRanges`, and `IRanges` which are
listed as `Suggests`. Install them via pixi:

    pixi global install --environment r-base \
      bioconductor-gviz bioconductor-genomicranges \
      bioconductor-iranges

or via `BiocManager`:

    BiocManager::install(c("Gviz", "GenomicRanges", "IRanges"))

For the optional ideogram and gene-region tracks add
`TxDb.Hsapiens.UCSC.hg38.knownGene` and `org.Hs.eg.db`.

## Lighter alternative

For the standard PIP + per-CS effect curve figure (no Gviz
dependency stack), see `vignette("fsusie_gtex_case_study")` and
`?mfsusieR::mfsusie_plot`.
