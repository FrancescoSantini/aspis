#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})

# ---- Parse arguments ----
option_list <- list(
  make_option(c("-r", "--results"), type = "character", help = "Full DESeq2 results CSV"),
  make_option(c("-f", "--filtered"), type = "character", help = "Filtered DESeq2 results CSV"),
  make_option(c("-c", "--counts"), type = "character", help = "Path to counts.txt"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to phenodata.csv"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory"),
  make_option(c("--padj"), type = "double", default = 0.1, help = "Adjusted p-value threshold"),
  make_option(c("--logfc"), type = "double", default = 1.0, help = "Absolute log2FC threshold")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Load DE results ----
res <- read.csv(opt$results, row.names = 1)
res_filtered <- read.csv(opt$filtered, row.names = 1)
sig_exists <- nrow(res_filtered) > 0
basename_stub <- sub("\\.csv$", "", basename(opt$results))

plot_padj <- opt$padj
plot_logfc <- opt$logfc

# Filter out rows with NA in padj or log2FC BEFORE plotting
res <- na.omit(res)

# Initialize columns safely
res$Significant <- rep(FALSE, nrow(res))
res$label <- rep("", nrow(res))

# Define significance
res$Significant <- with(res, padj < plot_padj & abs(log2FoldChange) > plot_logfc)

# Label top hits
if (sum(res$Significant) > 0) {
  top_sig <- head(rownames(res[res$Significant, , drop = FALSE]), 10)
  res$label[rownames(res) %in% top_sig] <- top_sig
} else {
  cat("[INFO] Filtered results are empty — will highlight top 10 miRNAs by padj for plotting.\n")
  top10 <- head(order(res$padj), 10)
  res$Significant[top10] <- TRUE
  res$label[top10] <- rownames(res)[top10]
}

# ---- Volcano Plot ----
volcano_title <- if (sum(res$Significant) > 0) {
  sprintf("Volcano Plot (padj < %.2f, |log2FC| > %.2f)", plot_padj, plot_logfc)
} else {
  "Top 10 miRNAs by padj (no hits at current thresholds)"
}

pdf(file.path(opt$outdir, paste0(basename_stub, "_volcano.pdf")))
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(aes(label = label), color = "black", size = 3, max.overlaps = Inf) +
  theme_minimal() +
  labs(title = volcano_title, x = "log2 Fold Change", y = "-log10 Adjusted P-Value")
dev.off()

# ---- Load counts and metadata ----
counts <- read.table(opt$counts, header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[, 6:ncol(counts)]
counts <- round(counts)
colnames(counts) <- gsub("\\.bam$", "", basename(colnames(counts)))

coldata <- read.csv(opt$metadata, row.names = 1)
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(counts)))

# ---- Add time_numeric_centered if biosample is numeric ----
if (all(grepl("^biosampletime\\d+$", coldata$biosample))) {
  cat("[INFO] Using centered numeric time from biosample\n")
  coldata$time_numeric <- as.numeric(gsub("biosampletime", "", coldata$biosample))
  coldata$time_numeric_centered <- scale(coldata$time_numeric, center = TRUE, scale = FALSE)
} else {
  cat("[INFO] Skipping time_numeric_centered: biosample not in expected format\n")
}

# ---- PCA ----
library(limma)

if (!"time_numeric_centered" %in% colnames(coldata)) {
  stop("time_numeric_centered is not available in coldata — check earlier parsing step.")
}
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ time_numeric_centered + condition)
vsd <- vst(dds, blind = FALSE)

# Remove batch effects based on time (if relevant)
assay(vsd) <- removeBatchEffect(assay(vsd), covariates = vsd$time_numeric_centered)

pcaData <- plotPCA(vsd, intgroup = c("condition", "time_numeric"), returnData = TRUE)
pcaData$Sample <- rownames(pcaData)
pcaData$Timepoint <- factor(pcaData$time_numeric,
                            levels = c(24, 48, 72),
                            labels = c("24h", "48h", "72h"))

percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_title <- sprintf("PCA (padj < %.2f, |log2FC| > %.2f) - %s miRNAs",
                     plot_padj, plot_logfc,
                     if (sig_exists) nrow(res_filtered) else "top 30 fallback")

pdf(file.path(opt$outdir, paste0(basename_stub, "_PCA.pdf")))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = Timepoint)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), color = "black", size = 3, max.overlaps = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(pca_title) +
  theme_minimal()
dev.off()


# ---- Heatmap ----

if (all(grepl("^biosampletime\\d+$", coldata$biosample))) {
  cat("[INFO] Parsing numeric time from biosample for heatmap ordering\n")
  coldata$time_numeric <- as.numeric(gsub("biosampletime", "", coldata$biosample))
} else {
  cat("[INFO] biosample format not numeric — skipping time ordering\n")
  coldata$time_numeric <- NA
}

# Filter to significant genes or fallback
res_filtered <- na.omit(res_filtered)
top_genes <- if (sig_exists) {
  rownames(res_filtered)
} else {
  head(rownames(res[order(res$padj), ]), 30)
}

# Build expression matrix
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

# Reorder columns by condition and time if available
if (!all(is.na(coldata$time_numeric))) {
  ordering <- order(coldata$condition, coldata$time_numeric)
  mat <- mat[, ordering]
  annotation <- coldata[ordering, c("condition", "biosample"), drop = FALSE]
} else {
  annotation <- coldata[, c("condition", "biosample"), drop = FALSE]
}

heat_colors <- colorRampPalette(c("green", "black", "red"))(100)

# Title fallback logic
heatmap_title <- if (sig_exists) {
  sprintf("Heatmap (padj < %.2f, |log2FC| > %.2f)", plot_padj, plot_logfc)
} else {
  "Top 30 miRNAs by padj (no hits at current thresholds)"
}

pdf(file.path(opt$outdir, paste0(basename_stub, "_heatmap.pdf")))
pheatmap(mat,
         annotation_col = annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         color = heat_colors,
         fontsize_row = 8,
         fontsize_col = 9,
         main = heatmap_title)
dev.off()

