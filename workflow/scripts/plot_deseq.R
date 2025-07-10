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

# ---- Enhanced Volcano Plot ----

# Define significance categories
res$SigCategory <- "NS"
res$SigCategory[abs(res$log2FoldChange) > plot_logfc] <- "FC"
res$SigCategory[res$padj < plot_padj] <- "FDR"
res$SigCategory[res$padj < plot_padj & abs(res$log2FoldChange) > plot_logfc] <- "FC_FDR"
res$SigCategory <- factor(res$SigCategory, levels = c("NS", "FC", "FDR", "FC_FDR"))

# Initialize label column
res$label <- rep("", nrow(res))

# Add labels to top genes (only FC+FDR)
if (any(res$SigCategory == "FC_FDR")) {
  top_sig <- head(rownames(res[res$SigCategory == "FC_FDR", , drop = FALSE]), 10)
  res$label[rownames(res) %in% top_sig] <- top_sig
} else {
  cat("[INFO] No FC+FDR hits — labeling top 10 by padj\n")
  top10 <- head(order(res$padj), 10)
  res$label[top10] <- rownames(res)[top10]
}

# Define legend title
volcano_title <- if (any(res$SigCategory == "FC_FDR")) {
  sprintf("Volcano Plot (padj < %.2f, |log2FC| > %.2f)", plot_padj, plot_logfc)
} else {
  "Top miRNAs by padj and/or FC (no hits in intersection)"
}

# Plot
pdf(file.path(opt$outdir, paste0(basename_stub, "_volcano.pdf")))
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = SigCategory)) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(
    values = c(NS = "grey30", FC = "forestgreen", FDR = "royalblue", FC_FDR = "red2"),
    labels = c(
      NS = "NS",
      FC = paste("LogFC > |", plot_logfc, "|"),
      FDR = paste("FDR < ", plot_padj),
      FC_FDR = paste("FDR < ", plot_padj, " & LogFC > |", plot_logfc, "|")
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = volcano_title,
    x = expression(log[2]~"fold change"),
    y = expression(-log[10]~adjusted~italic(P))
  ) +
  geom_text_repel(
    data = subset(res, label != ""),
    aes(label = label),
    color = "black", size = 2.5, max.overlaps = Inf
  ) +
  geom_vline(xintercept = c(-plot_logfc, plot_logfc), linetype = "longdash", color = "black", size = 0.4) +
  geom_hline(yintercept = -log10(plot_padj), linetype = "longdash", color = "black", size = 0.4)
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

coldata$Timepoint <- factor(coldata$time_numeric,
                            levels = c(24, 48, 72),
                            labels = c("24h", "48h", "72h"))

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
  annotation <- coldata[ordering, c("condition", "Timepoint"), drop = FALSE]
} else {
  annotation <- coldata[, c("condition", "Timepoint"), drop = FALSE]
}

annotation_colors <- list(
  condition = c(
    "control" = "#1f77b4",    # blue
    "PS_10uM" = "#2ca02c",    # green
    "PS_100uM" = "#d62728"    # red
  ),
  Timepoint = c(
    "24h" = "#a6cee3",  # light blue
    "48h" = "#1f78b4",  # mid blue
    "72h" = "#b2df8a"   # light green
  )
)

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
         annotation_colors = annotation_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         color = heat_colors,
         fontsize_row = 8,
         fontsize_col = 9,
         main = heatmap_title)
dev.off()

