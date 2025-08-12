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
  make_option(c("-r", "--results"), type = "character"),
  make_option(c("-f", "--filtered"), type = "character"),
  make_option(c("-c", "--counts"), type = "character"),
  make_option(c("-m", "--metadata"), type = "character"),
  make_option(c("-o", "--outdir"), type = "character"),
  make_option(c("--padj"), type = "double", default = 0.1),
  make_option(c("--logfc"), type = "double", default = 1.0),
  make_option(c("--topn"), type = "integer", default = 30),
  make_option(c("--annotate"), type = "character", default = NULL),
  make_option("--label", type = "character", default = "features")
)
opt <- parse_args(OptionParser(option_list = option_list))
plot_padj <- opt$padj
plot_logfc <- opt$logfc

# ---- Load DE results ----
res <- read.csv(opt$results, row.names = 1)
res_filtered <- read.csv(opt$filtered, row.names = 1)

# ---- Clean up rownames ----
rownames(res) <- gsub('^"|"$', '', rownames(res))
rownames(res_filtered) <- gsub('^"|"$', '', rownames(res_filtered))

cat(sprintf("[DEBUG] DESeq result table rows: %d (filtered: %d)\n", nrow(res), nrow(res_filtered)))

# ---- Annotate GeneType (Known / Novel) ----
if (!is.null(opt$annotate)) {
  cat("[DEBUG] Annotating GeneType using transcript metadata file\n")
  meta <- read.csv(opt$annotate)
  gene_id_map <- meta$gene_id
  names(gene_id_map) <- meta$transcript_id

  matched_gene_ids <- gene_id_map[rownames(res)]
  cat("[DEBUG] matched_gene_ids NA count:\n")
  print(sum(is.na(matched_gene_ids)))

  gene_type <- ifelse(grepl("^MSTRG", matched_gene_ids), "Novel", "Known")
  res$GeneType <- factor(gene_type, levels = c("Known", "Novel"))

  res_filtered$GeneType <- res$GeneType[match(rownames(res_filtered), rownames(res))]
} else {
  cat("[DEBUG] Annotating GeneType using rowname pattern (MSTRG)\n")
  gene_type <- ifelse(grepl("^MSTRG", rownames(res)), "Novel", "Known")
  res$GeneType <- factor(gene_type, levels = c("Known", "Novel"))

  res_filtered$GeneType <- res$GeneType[match(rownames(res_filtered), rownames(res))]
}

cat("[DEBUG] GeneType counts in full table:\n")
print(table(res$GeneType, useNA = "ifany"))

cat("[DEBUG] GeneType counts in filtered table:\n")
print(table(res_filtered$GeneType, useNA = "ifany"))

cat("[DEBUG] Unique GeneType values in filtered:\n")
print(unique(res_filtered$GeneType))

# ---- Load counts and metadata ----
counts <- read.table(opt$counts, header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[, !(colnames(counts) %in% c("Chr", "Start", "End", "Strand", "Length"))]
counts <- round(counts)
colnames(counts) <- gsub("_sorted$", "", gsub("\\.bam$", "", basename(colnames(counts))))

coldata <- read.csv(opt$metadata, row.names = 1)
coldata$condition <- factor(coldata$condition)
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(counts)))

# ---- Parse condition + covariate
contrast_folder <- basename(dirname(opt$results))
bioproject <- basename(dirname(dirname(opt$results)))
parts <- strsplit(contrast_folder, "_vs_control_")[[1]]
stopifnot(length(parts) == 2)
condition_of_interest <- parts[1]
covariate_value <- as.numeric(parts[2])
title_prefix <- sprintf("%s - %s_vs_control_%s", bioproject, condition_of_interest, covariate_value)

# ---- Subset coldata/counts
coldata <- coldata[
  coldata$condition %in% c("control", condition_of_interest) &
  as.character(coldata$covariate1) == as.character(covariate_value), , drop = FALSE
]
counts <- counts[, rownames(coldata)]

# ---- Handle timepoint
has_time <- FALSE
if ("covariate1" %in% colnames(coldata)) {
  coldata$time_numeric <- suppressWarnings(as.numeric(as.character(coldata$covariate1)))
  if (!all(is.na(coldata$time_numeric))) {
    has_time <- TRUE
    coldata$Timepoint <- factor(coldata$time_numeric,
      levels = sort(unique(coldata$time_numeric)),
      labels = paste0(sort(unique(coldata$time_numeric)), "h")
    )
  }
}

# ---- VST and PCA
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = FALSE)
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
vsd_top <- vsd[top_var_genes, ]
basename_stub <- sub("\\.csv$", "", basename(opt$results))
write.csv(assay(vsd), file.path(opt$outdir, paste0(basename_stub, "_vst.csv")))

# ---- Setup annotation colors
condition_levels <- unique(coldata$condition)
condition_colors <- setNames(rep(c("#1f77b4", "#2ca02c", "#d62728", "#ff7f0e"), length.out = length(condition_levels)), condition_levels)
annotation_colors <- list(condition = condition_colors)
if (!all(is.na(coldata$Timepoint))) {
  time_levels <- unique(coldata$Timepoint)
  time_colors <- setNames(rep(c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99"), length.out = length(time_levels)), time_levels)
  annotation_colors$Timepoint <- time_colors
}

# ---- Determine subset strategy based on identifier prefixes
id_prefixes <- substr(rownames(res), 1, 5)
has_ensg <- any(grepl("^ENSG", id_prefixes))
has_mstrg <- any(grepl("^MSTRG", id_prefixes))
has_gene_type_annotation <- any(res$GeneType %in% c("Known", "Novel"))

subset_types <- if (has_gene_type_annotation && (has_ensg || has_mstrg)) {
  c("All", "Known", "Novel")
} else {
  message("[INFO] Gene novelty split disabled — using only 'All' mode.")
  c("All")
}

# ---- Volcano plots ----
pdf(file.path(opt$outdir, paste0(basename_stub, "_volcano.pdf")))
for (subset_type in subset_types) {
  subres <- switch(subset_type,
    Known = res[res$GeneType == "Known", ],
    Novel = res[res$GeneType == "Novel", ],
    res
  )

  subres <- na.omit(subres)
  subres$SigCategory <- "NS"
  subres$SigCategory[abs(subres$log2FoldChange) > opt$logfc] <- "FC"
  subres$SigCategory[subres$padj < opt$padj] <- "FDR"
  subres$SigCategory[subres$padj < opt$padj & abs(subres$log2FoldChange) > opt$logfc] <- "FC_FDR"
  subres$SigCategory <- factor(subres$SigCategory, levels = c("NS", "FC", "FDR", "FC_FDR"))

  subres$label <- ""

  # Try to use gene_name from metadata (if available)
  gene_name_map <- NULL
  if (!is.null(opt$annotate)) {
    meta <- read.csv(opt$annotate)
    meta <- meta[!is.na(meta$transcript_id) & !is.na(meta$gene_name), ]
    gene_name_map <- setNames(meta$gene_name, meta$transcript_id)
  }

  if (any(subres$SigCategory == "FC_FDR")) {
    top_fc_fdr <- head(rownames(subres[subres$SigCategory == "FC_FDR", ]), opt$topn)
    if (!is.null(gene_name_map)) {
      subres$label[rownames(subres) %in% top_fc_fdr] <- gene_name_map[rownames(subres)[rownames(subres) %in% top_fc_fdr]]
    } else {
      subres$label[rownames(subres) %in% top_fc_fdr] <- sub(".*\\|", "", top_fc_fdr)
    }
  } else {
    top_padjs <- head(order(subres$padj), opt$topn)
    if (!is.null(gene_name_map)) {
      subres$label[top_padjs] <- gene_name_map[rownames(subres)[top_padjs]]
    } else {
      subres$label[top_padjs] <- rownames(subres)[top_padjs]
    }
  }

  # Fallback: remove empty labels caused by missing names
  subres$label[is.na(subres$label)] <- ""

  plot_title <- sprintf("Volcano Plot - %s subset\n%d features\npadj < %.2f, |log2FC| > %.2f",
                      subset_type, nrow(subres), plot_padj, plot_logfc)

  print(
    ggplot(subres, aes(x = log2FoldChange, y = -log10(padj), color = SigCategory, shape = GeneType)) +
      geom_point(alpha = 0.6, size = 0.8) +
      scale_color_manual(values = c(NS = "grey30", FC = "forestgreen", FDR = "royalblue", FC_FDR = "red2")) +
      scale_shape_manual(values = c(Known = 16, Novel = 17)) +
      theme_bw(base_size = 14) +
      theme(legend.position = "top", legend.title = element_blank(), legend.direction = "vertical") +
      labs(
        title = paste(title_prefix, plot_title, sep = "\n"),
        x = expression(log[2]~"fold change"),
        y = expression(-log[10]~adjusted~italic(P))
      ) +
      geom_text_repel(data = subset(subres, label != ""), aes(label = label),
                      color = "black", size = 2.5, segment.color = "grey60", segment.alpha = 0.5,
                      max.overlaps = Inf) +
      geom_vline(xintercept = c(-opt$logfc, opt$logfc), linetype = "longdash", color = "black", linewidth = 0.4) +
      geom_hline(yintercept = -log10(opt$padj), linetype = "longdash", color = "black", linewidth = 0.4)
  )
}
dev.off()

# ---- PCA plot ----
pcaData <- if (has_time) {
  plotPCA(vsd_top, intgroup = c("condition", "time_numeric"), returnData = TRUE)
} else {
  plotPCA(vsd_top, intgroup = "condition", returnData = TRUE)
}
pcaData$Sample <- rownames(pcaData)
pcaData$Timepoint <- if (has_time) factor(coldata$Timepoint) else "NA"
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_title <- sprintf("PCA\ntop 500 variable %s", opt$label)

pdf(file.path(opt$outdir, paste0(basename_stub, "_PCA.pdf")))
ggplot(pcaData, aes(PC1, PC2, color = condition, shape = Timepoint)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), color = "black", size = 3, max.overlaps = Inf) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle(paste(title_prefix, pca_title, sep = "\n")) +
  theme_minimal()
dev.off()

# ---- Heatmaps ----
pdf(file.path(opt$outdir, paste0(basename_stub, "_heatmap.pdf")))
for (subset_type in c("All", "Known", "Novel")) {
  top_rows <- switch(subset_type,
    Known = res_filtered[res_filtered$GeneType == "Known", ],
    Novel = res_filtered[res_filtered$GeneType == "Novel", ],
    res_filtered
  )

  if (nrow(top_rows) < 2) {
    plot.new()
    title(paste(title_prefix, "Heatmap skipped - too few", subset_type, "features"), cex.main = 1.2)
    next
  }

  top_genes <- rownames(head(top_rows[order(top_rows$padj), ], opt$topn))
  top_genes <- intersect(top_genes, rownames(assay(vsd)))
  if (length(top_genes) < 2) {
    plot.new()
    title(paste(title_prefix, "Heatmap skipped - too few overlapping genes:", subset_type), cex.main = 1.2)
    next
  }

  mat <- assay(vsd)[top_genes, , drop = FALSE]
  mat <- mat - rowMeans(mat)

  # Clean up gene-level rownames: use only gene_name if available
  rownames(mat) <- sub(".*\\|", "", rownames(mat))

  annotation <- if (has_time) {
    ordering <- order(coldata$condition, coldata$time_numeric)
    mat <- mat[, ordering]
    coldata[ordering, c("condition", "Timepoint"), drop = FALSE]
  } else {
    coldata[, "condition", drop = FALSE]
  }

  heatmap_n <- nrow(top_rows)
  heatmap_main <- sprintf("Heatmap (%s - top %d/%d %s)", subset_type, length(top_genes), heatmap_n, opt$label)

  # Try to use gene names as heatmap rownames (if annotate provided)
  if (!is.null(opt$annotate)) {
    meta <- read.csv(opt$annotate)
    meta_named <- setNames(meta$gene_name, meta$transcript_id)
    valid_rows <- rownames(mat) %in% names(meta_named)
    rownames(mat)[valid_rows] <- ifelse(
      !is.na(meta_named[rownames(mat)[valid_rows]]),
      meta_named[rownames(mat)[valid_rows]],
      rownames(mat)[valid_rows]
    )
  }

  pheatmap(mat,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    color = colorRampPalette(c("green", "black", "red"))(100),
    fontsize_row = 8,
    fontsize_col = 9,
    main = paste(title_prefix, heatmap_main, sep = "\n")
  )
}
dev.off()
