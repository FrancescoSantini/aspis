#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

# Argument parsing
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Path to count matrix"),
  make_option(c("-p", "--phenodata"), type = "character", help = "Path to phenodata CSV"),
  make_option(c("-o", "--output"), type = "character", help = "Output for default model results"),
  make_option(c("-x", "--interaction"), type = "character", help = "Output for interaction model results"),
  make_option("--filtered", type = "character", help = "Filtered output for default model"),
  make_option("--interaction_filtered", type = "character", help = "Filtered output for interaction model"),
  make_option("--padj", type = "double", default = 0.1, help = "Adjusted p-value threshold [default: %default]"),
  make_option("--logfc", type = "double", default = 1, help = "Absolute log2 fold change threshold [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Load input
counts <- read.table(opt$counts, header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[, 6:ncol(counts)]
counts <- round(counts)
colnames(counts) <- gsub("\\.bam$", "", basename(colnames(counts)))

coldata <- read.csv(opt$phenodata, row.names = 1)
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(counts)))

# Handle time as numeric if possible
if (all(grepl("^biosampletime\\d+$", coldata$biosample))) {
  cat("[INFO] Using centered numeric time from biosample\n")
  coldata$time_numeric <- as.numeric(gsub("biosampletime", "", coldata$biosample))
  coldata$time_numeric_centered <- scale(coldata$time_numeric, center = TRUE, scale = FALSE)
  time_term <- "time_numeric_centered"
  use_time_numeric <- TRUE
} else {
  cat("[INFO] Using biosample as factor\n")
  time_term <- "biosample"
  use_time_numeric <- FALSE
}

# -------- Default model --------
cat("[INFO] Running default model\n")
if (use_time_numeric) {
  dds1 <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                 design = as.formula("~ time_numeric_centered + condition"))
} else {
  dds1 <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                 design = ~ biosample + condition)
}
dds1 <- DESeq(dds1)
res1 <- lfcShrink(dds1, coef = 2, type = "apeglm")
res1 <- res1[order(res1$padj), ]
write.csv(as.data.frame(res1), file = opt$output)

cat("[INFO] Calculating alternative shrinkage methods for comparison\n")

# Raw log2FC (no shrinkage)
res_raw <- results(dds1)

# Normal shrinkage
res_norm <- lfcShrink(dds1, coef = 2, type = "normal")

# Combine into one table
res_compare <- data.frame(
  baseMean = res1$baseMean,
  padj = res1$padj,
  apeglm_log2FC = res1$log2FoldChange,
  normal_log2FC = res_norm$log2FoldChange,
  raw_log2FC = res_raw$log2FoldChange
)

rownames(res_compare) <- rownames(res1)

# Sort and export
res_compare <- res_compare[order(res_compare$padj), ]
write.csv(res_compare, file = sub("\\.csv$", "_logFCcomparison.csv", opt$output))

cat("[INFO] Wrote logFC comparison to _logFCcomparison.csv\n")

cat(sprintf("[INFO] Using padj threshold = %.3f and log2FC threshold = %.2f\n", opt$padj, opt$logfc))

# Main filtered output using both padj and log2FC (for pipeline consistency)
resSig1 <- subset(res1, padj < opt$padj & abs(log2FoldChange) > opt$logfc)
write.csv(as.data.frame(resSig1), file = opt$filtered)

# Optional: write a padj-only filtered table for exploratory use
resPadjOnly <- subset(res1, padj < opt$padj)
write.csv(as.data.frame(resPadjOnly),
          file = sub("\\.csv$", "_padjonly.csv", opt$filtered))  # adds _padjonly.csv

cat("[DEBUG] Default padj summary:\n")
print(summary(res1$padj))
cat("[DEBUG] Default log2FC summary:\n")
print(summary(res1$log2FoldChange))

# -------- Interaction model --------
cat("[INFO] Running interaction model\n")
if (use_time_numeric) {
  dds2 <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                 design = as.formula("~ time_numeric_centered * condition"))
} else {
  dds2 <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                                 design = ~ biosample * condition)
}
dds2 <- DESeq(dds2)
cat("[DEBUG] Available result names:\n")
print(resultsNames(dds2))

res2 <- lfcShrink(dds2, coef = 2, type = "apeglm")
res2 <- res2[order(res2$padj), ]
write.csv(as.data.frame(res2), file = opt$interaction)

resSig2 <- subset(res2, padj < opt$padj & abs(log2FoldChange) > opt$logfc)
write.csv(as.data.frame(resSig2), file = opt$interaction_filtered)

cat("[DEBUG] Interaction padj summary:\n")
print(summary(res2$padj))
cat("[DEBUG] Interaction log2FC summary:\n")
print(summary(res2$log2FoldChange))
