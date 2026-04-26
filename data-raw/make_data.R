# Build the bundled simulated example datasets in `data/`.
#
# Each dataset is simulated. The genotype matrix is sliced from
# `susieR::N3finemapping$X` to retain real LD; the response is
# generated from explicit known causal variants. Datasets cover
# DNA methylation along a CpG island, RNA-seq coverage along a
# gene body, and multi-outcome combinations of the two.
#
# Run from the package root with:
#   Rscript data-raw/make_data.R
# The resulting `.rda` files in `data/` are committed.

suppressMessages({
  library(susieR)
  library(wavethresh)
})

set.seed(20260426)

# --- DNAm: small simulated methylation QTL example ---
# n = 100, p = 12 SNPs, T = 32 CpGs. Three causal SNPs (1, 9, 3)
# act on two CpG clusters: SNPs 1 and 9 affect CpGs 9-16 with
# opposite signs; SNP 3 affects CpGs 25-32. SNP 4 is a high-LD
# near-clone of SNP 3. CpG positions are evenly spaced.

dnam_example <- local({
  set.seed(1)
  n <- 100L
  p <- 12L
  T_m <- 32L

  maf <- 0.05 + 0.45 * runif(p)
  X <- (matrix(runif(n * p), n, p, byrow = TRUE) < rep(maf, each = n)) +
       (matrix(runif(n * p), n, p, byrow = TRUE) < rep(maf, each = n))
  storage.mode(X) <- "double"
  # SNP 4 is a near-clone of SNP 3 (high-LD pair) so per-CpG
  # susie() ambiguates which of {SNP 3, SNP 4} drives the cluster B
  # signal, which is the comparison the methylation vignette uses.
  X[, 4L] <- X[, 3L] + 0.03 * rnorm(n)
  colnames(X) <- paste0("SNP", seq_len(p))

  # Effect matrix: SNPs 1, 9 -> CpG cluster A (positions 9-16),
  # SNP 3 -> CpG cluster B (positions 25-32).
  beta <- matrix(0, p, T_m)
  beta[1L,  9L:16L] <-  2.3
  beta[9L,  9L:16L] <- -2.3
  beta[3L, 25L:32L] <-  2.0

  E <- matrix(3 * rnorm(n * T_m), n, T_m)
  Y <- X %*% beta + E
  Y <- Y - min(Y)                     # nonneg methylation-style scale
  colnames(Y) <- paste0("CpG", seq_len(T_m))

  cpg_pos <- seq_len(T_m)

  # truth_mask is aligned with causal_snps = c(1, 9, 3):
  #   truth_mask[[1]] = CpGs affected by SNP 1 (cluster A)
  #   truth_mask[[2]] = CpGs affected by SNP 9 (cluster A)
  #   truth_mask[[3]] = CpGs affected by SNP 3 (cluster B)
  truth_mask <- list(
    seq.int(9L, 16L),                 # SNP 1
    seq.int(9L, 16L),                 # SNP 9
    seq.int(25L, 32L)                 # SNP 3
  )
  truth_mask <- lapply(truth_mask, function(idx) {
    m <- logical(T_m); m[idx] <- TRUE; m
  })

  list(
    X            = X,
    Y            = Y,
    pos          = cpg_pos,
    causal_snps  = c(1L, 9L, 3L),
    causal_betas = beta[c(1L, 9L, 3L), ],
    truth_mask   = truth_mask,
    description  = paste0(
      "Simulated methylation QTL example. n = ", n,
      " individuals, p = ", p, " SNPs, T = ", T_m,
      " evenly-spaced CpGs. Three causal SNPs (1, 9, 3): SNPs 1",
      " and 9 affect CpGs 9-16 with opposite signs; SNP 3 affects",
      " CpGs 25-32. SNP 4 is a high-LD near-clone of SNP 3.",
      " truth_mask is a length-3 list of length-T boolean",
      " vectors marking the CpGs each CS truly affects, read by",
      " mfsusie_plot_lfsr() when its truth argument is supplied."
    )
  )
})

# --- Coverage QTL: two smooth cis-effects on per-position ---
# Two causal SNPs at positions 25 and 75. Per-position effects
# are draws from `fsusieR::simu_IBSS_per_level(7)$sim_func`, the
# IBSS-prior wavelet sample on a length-128 grid. This matches
# the prior class fSuSiE assumes, so the recovered fitted
# function tracks the truth tightly. The example is
# assay-agnostic: the per-position response could be
# RNA-seq exon-body coverage, ATAC-seq peak coverage,
# WGBS / ChIP-seq read counts on a window, or any signal
# sampled at a power-of-two grid of positions.

coverage_example <- local({
  set.seed(1)
  data(N3finemapping)
  X <- N3finemapping$X[, seq_len(100)]
  n <- nrow(X); p <- ncol(X)
  T_m <- 128L
  exon_pos <- seq_len(T_m)

  rsnr <- 0.5
  pos1 <- 25L
  pos2 <- 75L
  # `simu_IBSS_per_level` is fsusieR's level-7 wavelet-prior
  # sampler; we use it once at data-build time so the bundled
  # dataset does not pull fsusieR at runtime.
  f1   <- fsusieR::simu_IBSS_per_level(7L)$sim_func
  f2   <- fsusieR::simu_IBSS_per_level(7L)$sim_func

  noise <- matrix(rnorm(n * T_m, sd = sqrt(var(f1) / rsnr)), n, T_m)
  Y     <- tcrossprod(X[, pos1], f1) + tcrossprod(X[, pos2], f2) + noise

  beta <- matrix(0, p, T_m)
  beta[pos1, ] <- f1
  beta[pos2, ] <- f2

  list(
    X            = X,
    Y            = Y,
    pos          = exon_pos,
    causal_snps  = c(pos1, pos2),
    causal_betas = beta[c(pos1, pos2), ],
    description  = paste0(
      "Simulated coverage-QTL example. n = ", n, ", p = ", p,
      " SNPs from the susieR::N3finemapping cis window, T = ",
      T_m, " evenly-spaced positions on a power-of-two grid.",
      " Two causal QTLs at SNPs ", pos1, " and ", pos2,
      " with smooth random per-position effects from the",
      " IBSS-prior wavelet sampler. The same shape applies",
      " to RNA-seq, ATAC-seq, WGBS, ChIP-seq, or any per-",
      "position coverage / read-count assay."
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
      "Simulated four-outcome cis-QTL example. n = ", n, ", p = ", p,
      " SNPs from the susieR::N3finemapping cis window. Outcomes:",
      " DNAm (T = 64, irregular CpGs), RNA-seq (T = 32, exon-body),",
      " and two scalar QTLs (eQTL, pQTL). Two causal SNPs (42, 97)",
      " shared across all four outcomes; per-outcome shapes and",
      " signs differ."
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
  peak_idx <- c(60L, 130L, 200L)
  peak_w   <- c(12L,  10L,   8L)
  beta <- matrix(0, p, T_m)
  beta[37,  ] <-  1.20 * shape(c = peak_idx[1L], w = peak_w[1L])
  beta[112, ] <- -0.85 * shape(c = peak_idx[2L], w = peak_w[2L])
  beta[173, ] <-  0.95 * shape(c = peak_idx[3L], w = peak_w[3L])

  # log1p-coverage observation model: a positive baseline count
  # plus the per-position effect, on the log1p scale.
  mean_cov <- 4 + X %*% beta
  Y <- mean_cov + matrix(rnorm(n * T_m, sd = 0.30), n)

  # Simulated transcript / exon annotation aligned with the
  # effect peaks. Real RNA-seq coverage signal comes from
  # transcribed (exonic) sequence, so each non-zero effect
  # region must sit on top of an exon. Three exons are placed
  # at the centers of the three peaks (sized to sit comfortably
  # inside each peak's main support); two short UTR-style exons
  # flank the gene at the 5' and 3' ends as decoration. Each
  # peak roughly covers its exon (the peak's support extends
  # well beyond the exon edges in both directions).
  bp_at <- function(idx) as.integer(round(exon_pos[1L] +
                                          (idx - 1L) * diff(range(exon_pos)) /
                                          (T_m - 1L)))
  peak_bp     <- vapply(peak_idx, bp_at, integer(1L))
  bp_per_idx  <- diff(range(exon_pos)) / (T_m - 1L)
  # Exon half-width = 2 * sigma in bp -- covers the peak's main
  # support (~95% of the Gaussian mass at +/- 2 sigma) so the
  # exon block is roughly as wide as the visible peak.
  half_w      <- as.integer(round(2 * peak_w * bp_per_idx))
  exon_starts <- c(109990500L,
                   peak_bp - half_w,
                   110021000L)
  exon_ends   <- c(109991300L,
                   peak_bp + half_w,
                   110021900L)
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
      "Simulated GTEx-style cis-eQTL example. n = ", n,
      ", p = ", p, " SNPs over the susieR::N3finemapping LD",
      " scaffold, T = ", T_m, " gene-body positions over a ~32 kb",
      " window. Three causal SNPs with localized peak effects at",
      " 5', mid-body, and 3' positions of the gene. The simulation",
      " reproduces the multi-CS structure of the published AHCYL1,",
      " SCD, and HSP90AA1 case studies; the underlying GTEx",
      " individual-level data are protected and cannot be",
      " redistributed."
    )
  )
})

usethis::use_data(dnam_example,      overwrite = TRUE, compress = "xz")
usethis::use_data(coverage_example,  overwrite = TRUE, compress = "xz")
usethis::use_data(multiomic_example, overwrite = TRUE, compress = "xz")
usethis::use_data(gtex_example,      overwrite = TRUE, compress = "xz")
