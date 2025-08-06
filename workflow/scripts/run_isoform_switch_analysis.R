#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DRIMSeq)
  library(IsoformSwitchAnalyzeR)
  library(readr)
  library(dplyr)
})

# --- CLI parsing ---
option_list <- list(
  make_option("--counts", type = "character", help = "Transcript count matrix (CSV)"),
  make_option("--phenodata", type = "character", help = "Phenodata file (CSV)"),
  make_option("--gtf", type = "character", help = "Merged GTF file"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--gene-expr", type = "double", default = 1, help = "Gene expression cutoff"),
  make_option("--isoform-expr", type = "double", default = 1, help = "Isoform expression cutoff"),
  make_option("--condition", type = "character", help = "Experimental condition to compare vs control"),
  make_option("--covariate", type = "character", help = "Covariate1 value to subset samples")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
cat("[INFO] Loading transcript counts...\n")
counts <- read_tsv(opt$counts)
colnames(counts)[1] <- "isoform_id"
colnames(counts) <- gsub(".*/", "", colnames(counts))
colnames(counts) <- gsub("\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("_sorted$", "", colnames(counts))

cat("[INFO] Loading phenodata...\n")
pheno <- read.csv(opt$phenodata, row.names = 1)
pheno$sampleID <- rownames(pheno)

# --- Subset by covariate and condition ---
if (!is.null(opt$covariate)) {
  cat(sprintf("[INFO] Subsetting to covariate1 == %s\n", opt$covariate))
  pheno <- pheno[pheno$covariate1 == opt$covariate, ]
}

if (!is.null(opt$condition)) {
  cat(sprintf("[INFO] Subsetting to condition in (control, %s)\n", opt$condition))
  pheno <- pheno[pheno$condition %in% c("control", opt$condition), ]
  pheno$condition <- factor(pheno$condition, levels = c("control", opt$condition))
}

# Drop constant columns (except sampleID and condition)
is_constant <- function(col) length(unique(na.omit(col))) <= 1
to_keep <- !sapply(pheno, is_constant) | names(pheno) %in% c("sampleID", "condition")
pheno <- pheno[, to_keep, drop = FALSE]
pheno <- pheno[, c("sampleID", "condition")]

cat("[INFO] Design matrix after final filtering:\n")
print(str(pheno))

# Match sample columns
samples <- intersect(colnames(counts)[-1], pheno$sampleID)
cat(sprintf("[INFO] Matched %d samples between counts and phenodata:\n", length(samples)))
print(samples)

counts <- counts[, c("isoform_id", samples)]
pheno <- pheno[pheno$sampleID %in% samples, ]

# --- Build switch list ---
cat("[INFO] Constructing switchAnalyzeRlist...\n")
switchList <- importRdata(
  isoformCountMatrix = counts,
  designMatrix = pheno,
  isoformExonAnno = opt$gtf,
  comparisonsToMake = NULL,
  showProgress = TRUE
)

# --- Pre-filter ---
switchListFiltered <- preFilter(
  switchList,
  removeSingleIsoformGenes = TRUE
)

# --- Manual expression filtering ---
cat(sprintf("[INFO] Applying manual expression filters: gene >= %.2f, isoform >= %.2f\n",
            opt$`gene-expr`, opt$`isoform-expr`))

expr_pass <- switchListFiltered$isoformFeatures$gene_expression >= opt$`gene-expr` &
             switchListFiltered$isoformFeatures$isoform_expression >= opt$`isoform-expr`

summary_path <- file.path(opt$outdir, "isoform_switch_summary.csv")
qc_plot_path <- file.path(opt$outdir, "isoform_switch_qc.pdf")

if (sum(expr_pass) > 0) {
  switchListFiltered <- subsetSwitchAnalyzeRlist(switchListFiltered, expr_pass)
} else {
  cat("[WARN] No isoforms pass the expression filter — writing empty result.\n")
  write.csv(data.frame(), summary_path, row.names = FALSE)
  file.create(qc_plot_path)
  file.create(file.path(opt$outdir, "switch_list.rds"))
  quit(save = "no", status = 0)
}

# --- Perform test ---
cat("[INFO] Running isoform switch test (SatuRn)...\n")
switchListTested <- isoformSwitchTestSatuRn(switchListFiltered)

# --- Save results ---
cat("[INFO] Extracting switch summary...\n")
write.csv(extractSwitchSummary(switchListTested), summary_path, row.names = FALSE)

# --- Optional QC plot ---
if (!is.null(switchListTested$isoformFeatures) &&
    nrow(switchListTested$isoformFeatures) > 0 &&
    !all(is.na(switchListTested$isoformFeatures$gene_id))) {

  pdf(qc_plot_path)
  switchPlot(switchListTested, gene = switchListTested$isoformFeatures$gene_id[1])
  dev.off()

} else {
  cat("[INFO] Skipping switchPlot: no valid gene for plotting.\n")
  file.create(qc_plot_path)
}

# --- Save full object ---
saveRDS(switchListTested, file = file.path(opt$outdir, "switch_list.rds"))
cat("[INFO] Analysis completed and results saved.\n")
