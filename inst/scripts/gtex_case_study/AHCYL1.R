# GTEx AHCYL1-style genome-browser case study (self-contained).
#
# Reproduces the multi-track effect / observed-difference / coverage
# figure layout published for AHCYL1 / SCD / HSP90AA1 / STMN2 /
# ABHD17A. The published figures use real GTEx individual-level
# data which cannot be redistributed; this script operates on the
# bundled simulated `gtex_example` instead and shows the same
# Gviz layout end-to-end.
#
# Run after installing the Bioconductor Suggests:
#   pixi global install --environment r-base bioconductor-gviz \
#     bioconductor-genomicranges bioconductor-iranges \
#     bioconductor-txdb.hsapiens.ucsc.hg38.knowngene \
#     bioconductor-org.hs.eg.db
# or, with BiocManager:
#   BiocManager::install(c("Gviz", "GenomicRanges", "IRanges",
#                          "TxDb.Hsapiens.UCSC.hg38.knownGene",
#                          "org.Hs.eg.db"))
#
# The optional TxDb / org.Hs.eg.db block at the bottom of this
# script overlays the gene-region track on a real chromosome.

library(mfsusieR)
suppressPackageStartupMessages({
  library(Gviz)
  library(GenomicRanges)
  library(IRanges)
})

data(gtex_example)

# ---- Fit fSuSiE -----------------------------------------------------

fit <- fsusie(gtex_example$Y, gtex_example$X,
              pos = gtex_example$pos,
              L = 15, L_greedy = 5,
              verbose = TRUE)
fit_s <- mf_post_smooth(fit, method = "TI",
                        wavelet_filter = 1L,
                        wavelet_family = "DaubExPhase")

# ---- Per-CS Gviz panel ----------------------------------------------

plot_cs <- function(cs) {
  markers <- fit_s$sets$cs[[cs]]
  lead    <- markers[which.max(fit_s$pip[markers])]
  x_lead  <- gtex_example$X[, lead]

  band   <- fit_s$credible_bands[[1L]][[cs]]
  effect <- fit_s$effect_curves[[1L]][[cs]]

  g0  <- which(x_lead <  median(x_lead))   # reference
  g1  <- which(x_lead >= median(x_lead))   # carrier
  mu0 <- colMeans(gtex_example$Y[g0, , drop = FALSE])
  mu1 <- colMeans(gtex_example$Y[g1, , drop = FALSE])

  chrom <- gtex_example$chrom
  pos   <- gtex_example$pos

  effect_track <- DataTrack(
    range = GRanges(chrom, IRanges(pos, pos + 1L)),
    data  = rbind(effect, band[, 1L], band[, 2L]),
    groups = c("effect", "lower", "upper"),
    type = "l", lty = c(1, 2, 2), lwd = c(2, 1, 1),
    col = c("#1f78b4", "#1f78b4", "#1f78b4"),
    name = sprintf("Effect CS %d", cs),
    ylim = c(min(band) * 1.1, max(band) * 1.1),
    background.title = "white", col.axis = "black",
    col.title = "black", fontface.title = 1, legend = FALSE)

  obs_diff_track <- DataTrack(
    range = GRanges(chrom, IRanges(pos, pos + 1L)),
    data  = mu1 - mu0,
    type = "p", pch = 16, cex = 0.5,
    col = "#1f78b4", name = "Observed difference",
    background.title = "white", col.axis = "black",
    col.title = "black", fontface.title = 1)

  count_track <- DataTrack(
    range = GRanges(chrom, IRanges(pos, pos + 1L)),
    data  = rbind(mu0, mu1),
    groups = c(sprintf("ref (n = %d)",     length(g0)),
               sprintf("carrier (n = %d)", length(g1))),
    type = "p", pch = 16, cex = 0.5,
    col = c("navy", "turquoise"),
    name = "Avg. log1p count",
    background.title = "white", col.axis = "black",
    col.title = "black", fontface.title = 1)

  axis_track <- GenomeAxisTrack(col = "black", fontcolor = "black",
                                col.title = "black",
                                background.title = "white")

  # Gene-region track from the simulated transcript / exon
  # annotation shipped with `gtex_example` (two isoforms plus a
  # 5' L1 element), drawn as exon boxes with intron lines.
  gene_track <- GeneRegionTrack(gtex_example$gene_track,
                                chromosome = chrom,
                                name = "", showId = TRUE,
                                geneSymbol = TRUE,
                                col = "salmon", fill = "salmon",
                                background.title = "white",
                                col.axis = "black",
                                col.title = "black")

  plotTracks(
    list(axis_track, effect_track, obs_diff_track, count_track,
         gene_track),
    from = gtex_example$locus[1L], to = gtex_example$locus[2L],
    sizes = c(1, 3, 3, 3, 2),
    main = sprintf("%s cis-eQTL: CS %d (lead SNP rank %d)",
                   gtex_example$gene, cs, lead),
    cex.main = 1.0)
}

# Plot each CS in turn. With three causal SNPs the simulated fit
# resolves three credible sets.
for (cs in seq_along(fit_s$sets$cs)) plot_cs(cs)

# ---- Optional: ideogram + real-TxDb gene-region track --------------
#
# When working with a real chromosome and a real TxDb annotation,
# swap the simulated `gtex_example$gene_track` for a `GeneRegionTrack`
# built off a TxDb package and add an `IdeogramTrack`:
#
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# ideo_track <- IdeogramTrack(genome = "hg38", chromosome = chrom)
# gene_track <- GeneRegionTrack(txdb, genome = "hg38",
#                               chromosome = chrom,
#                               from = gtex_example$locus[1L],
#                               to   = gtex_example$locus[2L],
#                               showId = TRUE, geneSymbol = TRUE,
#                               name = "", col = "salmon",
#                               fill = "salmon",
#                               background.title = "white")
# plotTracks(list(ideo_track, axis_track, effect_track,
#                 obs_diff_track, count_track, gene_track),
#            from = gtex_example$locus[1L],
#            to   = gtex_example$locus[2L],
#            sizes = c(1, 1, 3, 3, 3, 2))
