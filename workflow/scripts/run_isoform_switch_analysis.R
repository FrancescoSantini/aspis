#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(IsoformSwitchAnalyzeR)
  library(readr)
  library(dplyr)
})

# --- CLI parsing ---
option_list <- list(
  make_option("--counts", type = "character", help = "Transcript count matrix (CSV)"),
  make_option("--phenodata", type = "character", help = "Phenodata file (CSV)"),
  make_option("--gtf", type = "character", help = "Merged GTF file"),
  make_option("--outdir", type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
cat("[INFO] Loading transcript counts...\n")
counts <- read_csv(opt$counts)
colnames(counts)[1] <- "isoform_id"

cat("[INFO] Loading phenodata...\n")
pheno <- read.csv(opt$phenodata, row.names = 1)
pheno$sampleID <- rownames(pheno)

# Ensure columns match
samples <- intersect(names(counts)[-1], pheno$sampleID)
counts <- counts[, c("isoform_id", samples)]
pheno <- pheno[pheno$sampleID %in% samples, ]

# --- Build switch list ---
cat("[INFO] Constructing switchAnalyzeRlist...\n")
switchList <- importRdata(
  isoformCountMatrix = counts,
  isoformRepExpression = counts,
  designMatrix = pheno,
  isoformExonAnno = opt$gtf,
  comparisonsToMake = NULL,
  showProgress = TRUE
)

# --- Filter ---
switchListFiltered <- preFilter(
  switchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 1,
  removeSingleIsoformGenes = TRUE
)

# --- Perform switch analysis ---
cat("[INFO] Running isoform switch test (DRIMSeq)...\n")
switchListTested <- isoformSwitchTestDRIMSeq(switchListFiltered)

# --- Extract results ---
cat("[INFO] Extracting switch summary...\n")
summary_path <- file.path(opt$outdir, "isoform_switch_summary.csv")
write.csv(extractSwitchSummary(switchListTested), summary_path, row.names = FALSE)

# --- Optional plots ---
pdf(file.path(opt$outdir, "isoform_switch_qc.pdf"))
switchPlot(switchListTested, gene = switchListTested$isoformFeatures$gene_id[1])
dev.off()

# --- Save list ---
saveRDS(switchListTested, file = file.path(opt$outdir, "switch_list.rds"))
cat("[INFO] Analysis completed and results saved.\n")
