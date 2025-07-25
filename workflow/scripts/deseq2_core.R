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
  make_option("--filtered", type = "character", help = "Filtered output for default model"),
  make_option("--padj", type = "double", default = 0.1, help = "Adjusted p-value threshold"),
  make_option("--logfc", type = "double", default = 1, help = "Absolute log2 fold change threshold"),
  make_option("--condition", type = "character", help = "Condition to compare against control"),
  make_option("--covariate", type = "character", help = "Covariate1 value to subset samples")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Load and preprocess data
counts <- read.table(opt$counts, header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[, !(colnames(counts) %in% c("Chr", "Start", "End", "Strand", "Length"))]
counts <- round(counts)
colnames(counts) <- gsub("\\.bam$", "", basename(colnames(counts)))

coldata <- read.csv(opt$phenodata, row.names = 1, na.strings = c("NA", ""))
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(counts)))

# Subset to the relevant covariate1
if (!is.null(opt$covariate)) {
  cat(sprintf("[INFO] Subsetting to covariate1 == %s\n", opt$covariate))
  coldata <- coldata[coldata$covariate1 == opt$covariate, ]
  counts <- counts[, rownames(coldata)]
}

# Subset to just control vs selected condition
if (!is.null(opt$condition)) {
  cat(sprintf("[INFO] Comparing condition: control vs %s\n", opt$condition))
  coldata <- coldata[coldata$condition %in% c("control", opt$condition), ]
  counts <- counts[, rownames(coldata)]
  coldata$condition <- factor(coldata$condition, levels = c("control", opt$condition))
}

# Default model (~ condition)
cat("[INFO] Running DESeq2 with design: ~ condition\n")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

cat("[INFO] Available coefficients:\n")
print(resultsNames(dds))

coef_name <- resultsNames(dds)[2]
cat(sprintf("[INFO] Using coef=2 which corresponds to: %s\n", coef_name))

res <- lfcShrink(dds, coef = 2, type = "apeglm")
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file = opt$output)

# Additional outputs
cat("[INFO] Calculating alternative shrinkage methods for comparison\n")
res_raw <- results(dds)
res_norm <- lfcShrink(dds, coef = 2, type = "normal")

res_compare <- data.frame(
  baseMean = res$baseMean,
  padj = res$padj,
  apeglm_log2FC = res$log2FoldChange,
  normal_log2FC = res_norm$log2FoldChange,
  raw_log2FC = res_raw$log2FoldChange
)
rownames(res_compare) <- rownames(res)
res_compare <- res_compare[order(res_compare$padj), ]
write.csv(res_compare, file = sub("\\.csv$", "_logFCcomparison.csv", opt$output))

# Filtered output
cat(sprintf("[INFO] Using padj threshold = %.3f and log2FC threshold = %.2f\n", opt$padj, opt$logfc))
resSig <- subset(res, padj < opt$padj & abs(log2FoldChange) > opt$logfc)
write.csv(as.data.frame(resSig), file = opt$filtered)

resPadjOnly <- subset(res, padj < opt$padj)
write.csv(as.data.frame(resPadjOnly), file = sub("\\.csv$", "_padjonly.csv", opt$filtered))

cat("[DEBUG] padj summary:\n")
print(summary(res$padj))
cat("[DEBUG] log2FC summary:\n")
print(summary(res$log2FoldChange))
