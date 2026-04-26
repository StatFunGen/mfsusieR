# Build the bundled simulated example datasets in `data/`.
#
# These are SIMULATED fixtures shaped after the kinds of real
# datasets fSuSiE / mfSuSiE were designed for (DNA methylation
# along a CpG-island, RNA-seq coverage along a gene body, plus
# multi-outcome combinations). They are NOT derived from any
# subject-level data; the genotype matrix uses the
# `susieR::N3finemapping` LD matrix as a real-LD scaffold and the
# response is built from explicit known causal variants. Bundling
# the fixtures lets every vignette plot reproduce a known figure
# without depending on any external data download.
#
# Run from the package root with:
#   Rscript data-raw/make_data.R
# The resulting `.rda` files in `data/` are committed.

suppressMessages({
  library(susieR)
})

set.seed(20260426)

# --- DNAm-style: 200 CpGs over a 10-kb region with two cis-mQTLs ---

dnam_example <- local({
  data(N3finemapping)
  n <- nrow(N3finemapping$X)
  # Re-use susieR's real LD scaffold for a 200-SNP cis window.
  X <- N3finemapping$X[, seq(101, 300)]
  p <- ncol(X)

  # Irregular CpG positions across ~10 kb.
  T_m <- 64L
  cpg_pos <- sort(sample(1L:10000L, T_m))

  # Two causal SNPs with smooth-bump methylation effects; both
  # within the same CpG island so the bumps are at similar
  # positions but with opposite sign.
  shape <- function(center, width)
    exp(-((seq_along(cpg_pos) - center)^2) / (2 * width^2))
  beta <- matrix(0, p, T_m)
  beta[37,  ] <-  0.55 * shape(center = 32, width = 6)
  beta[112, ] <- -0.40 * shape(center = 22, width = 4)

  Y <- X %*% beta + matrix(rnorm(n * T_m, sd = 0.20), n)

  list(
    X            = X,
    Y            = Y,
    pos          = cpg_pos,
    causal_snps  = c(37L, 112L),
    causal_betas = beta[c(37L, 112L), ],
    description  = paste0(
      "Simulated DNAm-style cis-mQTL fixture. n = ", n, ", p = ", p,
      " SNPs from the susieR::N3finemapping cis window, T = ", T_m,
      " irregular CpG positions over a 10-kb region. Two causal",
      " variants with overlapping but oppositely-signed Gaussian",
      " methylation effects. Designed to surface fSuSiE's",
      " functional fine-mapping behaviour and post-smoothing",
      " (TI / HMM) credible bands."
    )
  )
})

# --- RNA-seq-style: 128 positions along a gene body, single eQTL ---

rnaseq_example <- local({
  data(N3finemapping)
  X <- N3finemapping$X[, seq(401, 600)]
  n <- nrow(X); p <- ncol(X)
  T_m <- 128L
  exon_pos <- seq_len(T_m)

  # One causal eQTL at SNP 73 with a peaked expression effect
  # concentrated over a 30-position window.
  shape <- exp(-((exon_pos - 70)^2) / (2 * 8^2))
  beta  <- matrix(0, p, T_m); beta[73, ] <- 1.5 * shape

  Y <- X %*% beta + matrix(rnorm(n * T_m, sd = 0.45), n)

  list(
    X            = X,
    Y            = Y,
    pos          = exon_pos,
    causal_snps  = 73L,
    causal_betas = beta[73L, , drop = FALSE],
    description  = paste0(
      "Simulated RNA-seq-style cis-eQTL fixture. n = ", n, ", p = ", p,
      " SNPs from the susieR::N3finemapping cis window, T = ", T_m,
      " evenly-spaced exon-body positions. One causal eQTL at SNP",
      " 73 with a peaked Gaussian effect spanning ~30 positions.",
      " Designed to surface a clean single-effect fit and the",
      " HMM smoother's lfsr overlay."
    )
  )
})

# --- Multi-outcome: DNAm (T=64) + RNA-seq (T=32) + 2 scalar QTLs ---

multiomic_example <- local({
  data(N3finemapping)
  X <- N3finemapping$X[, seq(201, 350)]
  n <- nrow(X); p <- ncol(X)

  # Two causal SNPs shared across all four outcomes.
  causal <- c(42L, 97L)

  T_dnam <- 64L
  cpg_pos <- sort(sample(1L:8000L, T_dnam))
  shape_d <- function(c, w) exp(-((seq_along(cpg_pos) - c)^2) / (2 * w^2))
  b_d <- matrix(0, p, T_dnam)
  b_d[42L, ] <-  0.60 * shape_d(c = 30L, w = 6)
  b_d[97L, ] <- -0.35 * shape_d(c = 22L, w = 5)
  Y_dnam <- X %*% b_d + matrix(rnorm(n * T_dnam, sd = 0.20), n)

  T_rna <- 32L
  exon_pos <- seq_len(T_rna)
  shape_r <- exp(-((exon_pos - 16)^2) / (2 * 4^2))
  b_r <- matrix(0, p, T_rna)
  b_r[42L, ] <-  1.0 * shape_r
  b_r[97L, ] <- -0.5 * shape_r
  Y_rna <- X %*% b_r + matrix(rnorm(n * T_rna, sd = 0.40), n)

  Y_eqtl  <- as.numeric(X[, 42L] *  0.9 + X[, 97L] * -0.5 +
                        rnorm(n, sd = 0.40))
  Y_pqtl  <- as.numeric(X[, 42L] *  0.7 + X[, 97L] * -0.6 +
                        rnorm(n, sd = 0.45))

  list(
    X        = X,
    Y_list   = list(
      dnam = Y_dnam,
      rna  = Y_rna,
      eqtl = matrix(Y_eqtl, ncol = 1L),
      pqtl = matrix(Y_pqtl, ncol = 1L)
    ),
    pos_list = list(
      dnam = cpg_pos,
      rna  = exon_pos,
      eqtl = 1L,
      pqtl = 1L
    ),
    causal_snps = causal,
    description = paste0(
      "Simulated four-outcome cis-QTL fixture. n = ", n, ", p = ", p,
      " SNPs from the susieR::N3finemapping cis window. Outcomes:",
      " DNAm (T = 64, irregular CpGs), RNA-seq (T = 32, exon-body),",
      " plus two scalar QTLs (eQTL, pQTL). Two causal SNPs (42, 97)",
      " shared across all four outcomes; shapes and signs differ",
      " per outcome. Designed to surface mfsusie()'s joint fine-",
      " mapping and the per-outcome panel layout in mfsusie_plot()."
    )
  )
})

# --- Save to data/ ---

# --- GTEx-style: gene-body RNA-seq coverage with multiple CSes ---

gtex_example <- local({
  data(N3finemapping)
  X <- N3finemapping$X[, seq(601, 800)]
  n <- nrow(X); p <- ncol(X)
  T_m <- 256L

  # Place the simulated locus at a realistic gene-body window:
  # chr1, 109,990,000 - 110,022,000 bp, mimicking the AHCYL1 cis
  # window. SNP positions are drawn within the same window.
  chrom    <- "chr1"
  gene     <- "AHCYL1-like"
  locus    <- c(109990000L, 110022000L)
  exon_pos <- round(seq(locus[1L], locus[2L], length.out = T_m))
  set.seed(1)
  snp_pos  <- sort(sample(seq(locus[1L] - 1.5e5, locus[2L] + 1.5e5),
                          p, replace = FALSE))
  set.seed(20260426)

  # Three causal SNPs at different gene-body positions, producing
  # the multi-CS structure characteristic of cis-eQTL fSuSiE fits
  # on real GTEx genes (AHCYL1, SCD, HSP90AA1, etc.). Effects are
  # localized peaks at distinct positions along the gene body.
  shape <- function(c, w) exp(-((seq_len(T_m) - c)^2) / (2 * w^2))
  beta <- matrix(0, p, T_m)
  beta[37,  ] <-  1.20 * shape(c =  60, w = 12)   # 5' peak
  beta[112, ] <- -0.85 * shape(c = 130, w = 10)   # mid-body peak
  beta[173, ] <-  0.95 * shape(c = 200, w =  8)   # 3' peak

  # log1p-coverage observation model: a positive baseline count
  # plus the per-position effect, on the log1p scale.
  mean_cov <- 4 + X %*% beta
  Y <- mean_cov + matrix(rnorm(n * T_m, sd = 0.30), n)

  # Simulated transcript / exon annotation for the gene-region
  # track. Two isoforms ("AHCYL1-like.1" - canonical; ".2" -
  # alternative-start variant skipping exon 1) over the locus,
  # 7 exons of varying widths, plus side L1 element on the 5' end
  # to mirror the published AHCYL1 figure layout.
  exon_starts <- c(109990500L, 109995200L, 109999800L, 110005000L,
                   110010000L, 110015500L, 110020500L)
  exon_ends   <- exon_starts + c(800L, 1100L, 600L, 900L, 1200L,
                                 700L, 1100L)
  gene_track_df <- rbind(
    data.frame(chromosome = chrom, start = exon_starts,
               end = exon_ends, strand = "+",
               feature = "protein_coding",
               gene = gene, exon = paste0("e", seq_along(exon_starts)),
               transcript = paste0(gene, ".1"),
               symbol = gene),
    # alternative isoform skipping exon 1
    data.frame(chromosome = chrom, start = exon_starts[-1L],
               end = exon_ends[-1L], strand = "+",
               feature = "protein_coding",
               gene = gene,
               exon = paste0("e", seq_along(exon_starts[-1L])),
               transcript = paste0(gene, ".2"),
               symbol = gene),
    # 5' L1 element annotation
    data.frame(chromosome = chrom, start = 109987000L,
               end = 109988200L, strand = "+",
               feature = "lincRNA",
               gene = "L1", exon = "e1",
               transcript = "L1.1", symbol = "L1")
  )

  list(
    X            = X,
    Y            = Y,
    pos          = exon_pos,
    snp_pos      = snp_pos,
    chrom        = chrom,
    gene         = gene,
    locus        = locus,
    gene_track   = gene_track_df,
    causal_snps  = c(37L, 112L, 173L),
    causal_betas = beta[c(37L, 112L, 173L), ],
    description  = paste0(
      "Simulated GTEx-style cis-eQTL fixture mimicking the AHCYL1 / ",
      "SCD / HSP90AA1 case studies. Designed for the GTEx ",
      "case-study vignette: the underlying GTEx data are protected ",
      "and cannot be redistributed, so we simulate a fixture that ",
      "reproduces the multi-CS structure published in the figures. ",
      "n = ", n, ", p = ", p, " SNPs over the susieR::N3finemapping ",
      "LD scaffold, T = ", T_m, " gene-body positions over a ~32 kb ",
      "window. Three causal SNPs with localized peak effects at 5', ",
      "mid-body, and 3' positions of the gene."
    )
  )
})

usethis::use_data(dnam_example,      overwrite = TRUE, compress = "xz")
usethis::use_data(rnaseq_example,    overwrite = TRUE, compress = "xz")
usethis::use_data(multiomic_example, overwrite = TRUE, compress = "xz")
usethis::use_data(gtex_example,      overwrite = TRUE, compress = "xz")
