#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(limma)
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

# ---- Load counts and metadata ----
counts <- read.table(opt$counts, header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[, 6:ncol(counts)]
counts <- round(counts)
colnames(counts) <- gsub("\\.bam$", "", basename(colnames(counts)))

coldata <- read.csv(opt$metadata, row.names = 1)
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(counts)))

# ---- Subset to relevant samples for this contrast ----
# Expected format: results/deg/{bioproject}/{condition}_vs_control_{covariate}/
contrast_folder <- basename(dirname(opt$results))
parts <- strsplit(contrast_folder, "_vs_control_")[[1]]

if (length(parts) == 2) {
  condition_of_interest <- parts[1]
  covariate_value <- parts[2]

  cat(sprintf("[INFO] Subsetting to condition = %s and covariate1 = %s\n", condition_of_interest, covariate_value))

  coldata <- coldata[
    coldata$condition %in% c("control", condition_of_interest) &
    coldata$covariate1 == covariate_value, , drop = FALSE
  ]

  counts <- counts[, rownames(coldata)]
} else {
  stop("Failed to parse condition and covariate from result folder name.")
}

cat("[INFO] Plotting PCA and heatmap for samples:\n")
print(rownames(coldata))

# ---- Process covariate1 (e.g., time) as numeric
has_time <- FALSE
if ("covariate1" %in% colnames(coldata)) {
  suppressWarnings({
    coldata$time_numeric <- as.numeric(as.character(coldata$covariate1))
  })
  if (!all(is.na(coldata$time_numeric))) {
    has_time <- TRUE
    coldata$Timepoint <- factor(coldata$time_numeric,
                                levels = sort(unique(coldata$time_numeric)),
                                labels = paste0(sort(unique(coldata$time_numeric)), "h"))
  }
}

# ---- VST transformation
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = FALSE)

# ---- Export VST matrix
write.csv(assay(vsd), file.path(opt$outdir, paste0(basename_stub, "_vst.csv")))

# ---- Write summary line to console
contrast_label <- basename(dirname(opt$results))
cat(sprintf("[SUMMARY] %s: %d significant genes (padj < %.2f & |log2FC| > %.2f)\n",
            contrast_label, nrow(res_filtered), plot_padj, plot_logfc))

# ---- Volcano Plot ----
res <- na.omit(res)
res$SigCategory <- "NS"
res$SigCategory[abs(res$log2FoldChange) > plot_logfc] <- "FC"
res$SigCategory[res$padj < plot_padj] <- "FDR"
res$SigCategory[res$padj < plot_padj & abs(res$log2FoldChange) > plot_logfc] <- "FC_FDR"
res$SigCategory <- factor(res$SigCategory, levels = c("NS", "FC", "FDR", "FC_FDR"))
res$label <- rep("", nrow(res))

if (any(res$SigCategory == "FC_FDR")) {
  top_sig <- head(rownames(res[res$SigCategory == "FC_FDR", , drop = FALSE]), 10)
  res$label[rownames(res) %in% top_sig] <- top_sig
} else {
  cat("[INFO] No FC+FDR hits — labeling top 10 by padj\n")
  top10 <- head(order(res$padj), 10)
  res$label[top10] <- rownames(res)[top10]
}

volcano_title <- if (any(res$SigCategory == "FC_FDR")) {
  sprintf("Volcano Plot (padj < %.2f, |log2FC| > %.2f)", plot_padj, plot_logfc)
} else {
  "Top miRNAs by padj and/or FC (no hits in intersection)"
}

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
    legend.title = element_blank()
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

# ---- PCA ----
if (has_time) {
  pcaData <- plotPCA(vsd, intgroup = c("condition", "time_numeric"), returnData = TRUE)
  pcaData$Sample <- rownames(pcaData)
  pcaData$Timepoint <- factor(pcaData$time_numeric,
                              levels = sort(unique(pcaData$time_numeric)),
                              labels = paste0(sort(unique(pcaData$time_numeric)), "h"))
} else {
  pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
  pcaData$Sample <- rownames(pcaData)
  pcaData$Timepoint <- "NA"
}

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
top_genes <- if (sig_exists) {
  rownames(na.omit(res_filtered))
} else {
  head(rownames(res[order(res$padj), ]), 30)
}
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

if (has_time) {
  ordering <- order(coldata$condition, coldata$time_numeric)
  mat <- mat[, ordering]
  annotation <- coldata[ordering, c("condition", "Timepoint"), drop = FALSE]
} else {
  annotation <- coldata[, "condition", drop = FALSE]
}

# Preferred palettes
preferred_condition_colors <- c("#1f77b4", "#2ca02c", "#d62728", "#ff7f0e")
preferred_time_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99")

condition_levels <- unique(coldata$condition)
condition_colors <- setNames(
  rep(preferred_condition_colors, length.out = length(condition_levels)),
  condition_levels
)

annotation_colors <- list(condition = condition_colors)

if (!all(is.na(coldata$Timepoint))) {
  time_levels <- unique(coldata$Timepoint)
  time_colors <- setNames(
    rep(preferred_time_colors, length.out = length(time_levels)),
    time_levels
  )
  annotation_colors$Timepoint <- time_colors
}

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
         color = colorRampPalette(c("green", "black", "red"))(100),
         fontsize_row = 8,
         fontsize_col = 9,
         main = heatmap_title)
dev.off()
