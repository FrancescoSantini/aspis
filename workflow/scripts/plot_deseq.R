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
basename_stub <- sub("\\.csv$", "", basename(opt$results))

sig_exists <- nrow(res_filtered) > 0
plot_padj <- opt$padj
plot_logfc <- opt$logfc

if (!sig_exists) {
  cat("[INFO] Filtered results are empty — using top 30 unfiltered genes for plotting.\n")
  plot_padj <- 0.2
  plot_logfc <- 0.5
}

# ---- Volcano Plot ----
res$Significant <- with(res, padj < plot_padj & abs(log2FoldChange) > plot_logfc)
res$label <- ""
if (sum(res$Significant, na.rm = TRUE) > 0) {
  top_sig <- head(rownames(res[res$Significant, ]), 10)
  res$label[top_sig] <- top_sig
}
volcano_title <- sprintf("Volcano Plot (padj < %.2f, |log2FC| > %.2f)", plot_padj, plot_logfc)

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

# ---- PCA ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("condition", "biosample"), returnData = TRUE)
pcaData$Sample <- rownames(pcaData)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_title <- sprintf("PCA (padj < %.2f, |log2FC| > %.2f) — %s genes",
                     plot_padj, plot_logfc,
                     if (sig_exists) nrow(res_filtered) else "top 30 fallback")

pdf(file.path(opt$outdir, paste0(basename_stub, "_PCA.pdf")))
ggplot(pcaData, aes(PC1, PC2, color = condition, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(color = "black", size = 3, max.overlaps = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(pca_title) +
  theme_minimal()
dev.off()

# ---- Heatmap ----
top_genes <- if (sig_exists) {
  rownames(res_filtered)
} else {
  head(rownames(res[order(res$padj), ]), 30)
}
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

annotation <- coldata[, c("condition", "biosample"), drop = FALSE]

heatmap_title <- sprintf("Heatmap (padj < %.2f, |log2FC| > %.2f)", plot_padj, plot_logfc)

pdf(file.path(opt$outdir, paste0(basename_stub, "_heatmap.pdf")))
pheatmap(mat,
         annotation_col = annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         main = heatmap_title)
dev.off()
