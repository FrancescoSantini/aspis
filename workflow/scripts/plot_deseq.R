#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(limma)
  library(dplyr)
})

# ------------------------- CLI -------------------------
option_list <- list(
  make_option(c("-r", "--results"),  type = "character"),
  make_option(c("-f", "--filtered"), type = "character"),
  make_option(c("-c", "--counts"),   type = "character"),
  make_option(c("-m", "--metadata"), type = "character"),
  make_option(c("-o", "--outdir"),   type = "character"),
  make_option(c("--padj"),  type = "double",  default = 0.1),
  make_option(c("--logfc"), type = "double",  default = 1.0),
  make_option(c("--topn"),  type = "integer", default = 30),
  make_option(c("--annotate"), type = "character", default = NULL),
  make_option("--label", type = "character", default = "features")
)
opt <- parse_args(OptionParser(option_list = option_list))
plot_padj  <- opt$padj
plot_logfc <- opt$logfc

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# ------------------------- Load DE results -------------------------
msg("[INFO] Reading DE tables (no rownames yet)")
res <- read.csv(opt$results,  check.names = FALSE)
res_filtered <- read.csv(opt$filtered, check.names = FALSE)

res_ids      <- gsub('^"|"$', '', res[[1]])
filtered_ids <- gsub('^"|"$', '', res_filtered[[1]])

if (anyDuplicated(res_ids)) {
  dups <- unique(res_ids[duplicated(res_ids)])
  stop(sprintf("Duplicate feature IDs in res: %s", paste(head(dups, 10), collapse = ", ")))
}
if (anyDuplicated(filtered_ids)) {
  dups <- unique(filtered_ids[duplicated(filtered_ids)])
  stop(sprintf("Duplicate feature IDs in res_filtered: %s", paste(head(dups, 10), collapse = ", ")))
}

rownames(res) <- res_ids;              res <- res[, -1, drop = FALSE]
rownames(res_filtered) <- filtered_ids; res_filtered <- res_filtered[, -1, drop = FALSE]

msg("[INFO] Rows: res=%d, filtered=%d", nrow(res), nrow(res_filtered))

# ------------------------- Annotation join -------------------------
# Strategy:
# 1) Build keys from rownames like "MSTRG.123|CFTR" or "ENSG...|CFTR"
# 2) Lookup order for GeneType/GeneName:
#    feature_id -> gene_id -> gene_name (right side after '|')
if (!is.null(opt$annotate) && file.exists(opt$annotate)) {
  annot <- read.csv(opt$annotate, check.names = FALSE)
  required_cols <- c("feature_id", "gene_id", "gene_name", "gene_type")
  missing_cols  <- setdiff(required_cols, colnames(annot))
  if (length(missing_cols)) {
    stop("Annotation missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # keys for res
  split_res <- strsplit(rownames(res), "\\|")
  res$feature_id_join <- vapply(split_res, `[`, "", 1)
  res$gene_id_join    <- vapply(split_res, function(x) if (length(x) >= 2) x[2] else NA_character_, "")
  right_name_res      <- sub(".*\\|", "", rownames(res))

  # keys for res_filtered
  split_filt <- strsplit(rownames(res_filtered), "\\|")
  res_filtered$feature_id_join <- vapply(split_filt, `[`, "", 1)
  res_filtered$gene_id_join    <- vapply(split_filt, function(x) if (length(x) >= 2) x[2] else NA_character_, "")
  right_name_filt              <- sub(".*\\|", "", rownames(res_filtered))

  # lookup maps
  gt_by_feature <- setNames(annot$gene_type, annot$feature_id)
  gt_by_gene    <- setNames(annot$gene_type, annot$gene_id)
  gt_by_name    <- setNames(annot$gene_type, annot$gene_name)  # fallback by symbol/name

  gn_by_feature <- setNames(annot$gene_name, annot$feature_id)
  gn_by_gene    <- setNames(annot$gene_name, annot$gene_id)

  # assign GeneType with fallback
  res$GeneType <- dplyr::coalesce(
    gt_by_feature[res$feature_id_join],
    gt_by_gene[res$gene_id_join],
    gt_by_name[right_name_res]
  )
  res_filtered$GeneType <- dplyr::coalesce(
    gt_by_feature[res_filtered$feature_id_join],
    gt_by_gene[res_filtered$gene_id_join],
    gt_by_name[right_name_filt]
  )

  # assign GeneName with fallback (ensures readable labels)
  res$GeneName <- dplyr::coalesce(
    gn_by_feature[res$feature_id_join],
    gn_by_gene[res$gene_id_join],
    right_name_res
  )
  res_filtered$GeneName <- dplyr::coalesce(
    gn_by_feature[res_filtered$feature_id_join],
    gn_by_gene[res_filtered$gene_id_join],
    right_name_filt
  )

  # finalize factors
  res$GeneType <- ifelse(is.na(res$GeneType), "Novel", res$GeneType)
  res_filtered$GeneType <- ifelse(is.na(res_filtered$GeneType), "Novel", res_filtered$GeneType)
  res$GeneType <- factor(res$GeneType, levels = c("Known", "Novel"))
  res_filtered$GeneType <- factor(res_filtered$GeneType, levels = c("Known", "Novel"))

  # concise debug
  gt_full  <- table(res$GeneType, useNA = "no")
  gt_filt  <- table(res_filtered$GeneType, useNA = "no")
  msg("[INFO] GeneType (full):   Known=%d Novel=%d", as.integer(gt_full["Known"]), as.integer(gt_full["Novel"]))
  msg("[INFO] GeneType (filtered): Known=%d Novel=%d", as.integer(gt_filt["Known"]), as.integer(gt_filt["Novel"]))
} else {
  msg("[WARN] --annotate missing; using MSTRG prefix fallback")
  res$GeneType <- factor(ifelse(grepl("^MSTRG", rownames(res)), "Novel", "Known"), levels = c("Known", "Novel"))
  res_filtered$GeneType <- res$GeneType[match(rownames(res_filtered), rownames(res))]
  res$GeneName <- sub(".*\\|", "", rownames(res))
  res_filtered$GeneName <- sub(".*\\|", "", rownames(res_filtered))
}

# ------------------------- Filtering stats -------------------------
n_non_na_padj <- sum(!is.na(res_filtered$padj))
n_padj        <- sum(res_filtered$padj < plot_padj, na.rm = TRUE)
n_fc          <- sum(abs(res_filtered$log2FoldChange) > plot_logfc, na.rm = TRUE)
n_both        <- sum(res_filtered$padj < plot_padj & abs(res_filtered$log2FoldChange) > plot_logfc, na.rm = TRUE)
msg("[INFO] Filter stats (filtered): nonNA padj=%d ; padj<%.3f=%d ; |log2FC|>%.2f=%d ; both=%d",
    n_non_na_padj, plot_padj, n_padj, plot_logfc, n_fc, n_both)

sig_mask <- !is.na(res_filtered$padj) &
            (res_filtered$padj < plot_padj) &
            !is.na(res_filtered$log2FoldChange) &
            (abs(res_filtered$log2FoldChange) > plot_logfc)
sig_res <- res_filtered[sig_mask, , drop = FALSE]
msg("[INFO] Significant rows carried into plots: %d", nrow(sig_res))

# ------------------------- Counts & coldata -------------------------
counts <- read.table(opt$counts, header = TRUE, row.names = 1, check.names = FALSE)
counts <- counts[, !(colnames(counts) %in% c("Chr", "Start", "End", "Strand", "Length")), drop = FALSE]
counts <- round(counts)
colnames(counts) <- gsub("_sorted$", "", gsub("\\.bam$", "", basename(colnames(counts))))

coldata <- read.csv(opt$metadata, row.names = 1)
coldata$condition <- factor(coldata$condition)
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(counts)))

# infer contrast from paths
contrast_folder <- basename(dirname(opt$results))
bioproject      <- basename(dirname(dirname(opt$results)))
parts <- strsplit(contrast_folder, "_vs_control_")[[1]]
stopifnot(length(parts) == 2)
condition_of_interest <- parts[1]
covariate_value       <- as.numeric(parts[2])
title_prefix <- sprintf("%s - %s_vs_control_%s", bioproject, condition_of_interest, covariate_value)

# subset to contrast
coldata <- coldata[
  coldata$condition %in% c("control", condition_of_interest) &
    as.character(coldata$covariate1) == as.character(covariate_value),
  , drop = FALSE
]
counts <- counts[, rownames(coldata), drop = FALSE]

# time/covariate
has_time <- FALSE
if ("covariate1" %in% colnames(coldata)) {
  coldata$time_numeric <- suppressWarnings(as.numeric(as.character(coldata$covariate1)))
  if (!all(is.na(coldata$time_numeric))) {
    has_time <- TRUE
    coldata$Timepoint <- factor(
      coldata$time_numeric,
      levels = sort(unique(coldata$time_numeric)),
      labels = paste0(sort(unique(coldata$time_numeric)), "h")
    )
  }
}

# ------------------------- VST / PCA prep -------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = FALSE)
top_var_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
vsd_top <- vsd[top_var_genes, ]
basename_stub <- sub("\\.csv$", "", basename(opt$results))
write.csv(assay(vsd), file.path(opt$outdir, paste0(basename_stub, "_vst.csv")))

# colors (unchanged)
condition_levels <- unique(coldata$condition)
condition_colors <- setNames(rep(c("#1f77b4", "#2ca02c", "#d62728", "#ff7f0e"),
                                 length.out = length(condition_levels)),
                             condition_levels)
annotation_colors <- list(condition = condition_colors)
if (!all(is.na(coldata$Timepoint))) {
  time_levels <- unique(coldata$Timepoint)
  time_colors <- setNames(rep(c("#a6cee3", "#1f78b4", "#b2df8a", "#fb9a99"),
                              length.out = length(time_levels)),
                          time_levels)
  annotation_colors$Timepoint <- time_colors
}

# novelty splits availability
has_gene_type_annotation <- any(res$GeneType %in% c("Known", "Novel"))
subset_types <- if (has_gene_type_annotation) c("All", "Known", "Novel") else c("All")

# ------------------------- Volcano -------------------------
pdf(file.path(opt$outdir, paste0(basename_stub, "_volcano.pdf")))
for (subset_type in subset_types) {
  subres <- switch(subset_type,
                   Known = res[res$GeneType == "Known", , drop = FALSE],
                   Novel = res[res$GeneType == "Novel", , drop = FALSE],
                   res)

  if (!nrow(subres)) {
    plot.new(); title(paste(title_prefix,
                            sprintf("Volcano Plot - %s subset\n(no features)", subset_type),
                            sep = "\n"))
    next
  }

  # Ensure numeric + drop rows missing padj/log2FC
  if ("padj" %in% colnames(subres)) subres$padj <- suppressWarnings(as.numeric(subres$padj))
  if ("log2FoldChange" %in% colnames(subres)) subres$log2FoldChange <- suppressWarnings(as.numeric(subres$log2FoldChange))
  keep <- !is.na(subres$padj) & !is.na(subres$log2FoldChange)
  subres <- subres[keep, , drop = FALSE]

  if (!nrow(subres)) {
    plot.new(); title(paste(title_prefix,
                            sprintf("Volcano Plot - %s subset\n(no valid rows after padj/log2FC filtering)", subset_type),
                            sep = "\n"))
    next
  }

  # Protect against padj == 0 (Inf on -log10)
  subres$padj_eps <- pmax(subres$padj, .Machine$double.xmin)

  # Significance categories
  subres$SigCategory <- "NS"
  subres$SigCategory[abs(subres$log2FoldChange) > opt$logfc] <- "FC"
  subres$SigCategory[subres$padj < opt$padj] <- "FDR"
  subres$SigCategory[subres$padj < opt$padj & abs(subres$log2FoldChange) > opt$logfc] <- "FC_FDR"
  subres$SigCategory <- factor(subres$SigCategory, levels = c("NS", "FC", "FDR", "FC_FDR"))

  # Labels (prefer FC_FDR; else smallest padj)
  subres$label <- ""
  sel <- if (any(subres$SigCategory == "FC_FDR")) which(subres$SigCategory == "FC_FDR") else order(subres$padj)
  sel <- head(sel, opt$topn)

  prettify <- function(id, gene_name) {
    if (!is.null(gene_name) && !is.na(gene_name) && gene_name != "") return(gene_name)
    if (grepl("\\|", id)) return(sub(".*\\|", "", id))
    id
  }
  subres$label[sel] <- mapply(
    prettify,
    id = rownames(subres)[sel],
    gene_name = if ("GeneName" %in% colnames(subres)) subres$GeneName[sel] else NA
  )
  subres$label[is.na(subres$label)] <- ""

  plot_title <- sprintf("Volcano Plot - %s subset\n%d features\npadj < %.2f, |log2FC| > %.2f",
                        subset_type, nrow(subres), plot_padj, plot_logfc)

  print(
    ggplot(subres, aes(x = log2FoldChange, y = -log10(padj_eps),
                       color = SigCategory, shape = GeneType)) +
      geom_point(alpha = 0.6, size = 0.8) +
      scale_color_manual(values = c(NS = "grey30", FC = "forestgreen",
                                    FDR = "royalblue", FC_FDR = "red2")) +
      scale_shape_manual(values = c(Known = 16, Novel = 17)) +
      theme_bw(base_size = 14) +
      theme(legend.position = "top", legend.title = element_blank(),
            legend.direction = "vertical") +
      labs(
        title = paste(title_prefix, plot_title, sep = "\n"),
        x = expression(log[2]~"fold change"),
        y = expression(-log[10]~adjusted~italic(P))
      ) +
      geom_text_repel(
        data = subset(subres, label != ""),
        aes(label = label, y = -log10(padj_eps)),   # <-- fix
        color = "black", size = 2.5,
        segment.color = "grey60", segment.alpha = 0.5,
        max.overlaps = Inf
      ) +
      geom_vline(xintercept = c(-opt$logfc, opt$logfc),
                 linetype = "longdash", color = "black", linewidth = 0.4) +
      geom_hline(yintercept = -log10(opt$padj),
                 linetype = "longdash", color = "black", linewidth = 0.4)
  )
}
dev.off()

# ------------------------- PCA -------------------------
pcaData <- if (has_time) {
  plotPCA(vsd_top, intgroup = c("condition", "time_numeric"), returnData = TRUE)
} else {
  plotPCA(vsd_top, intgroup = "condition", returnData = TRUE)
}
pcaData$Sample <- rownames(pcaData)
pcaData$Timepoint <- if (has_time) factor(coldata$Timepoint) else "NA"
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_title <- sprintf("PCA\nTop 500 most variable %s", opt$label)

pdf(file.path(opt$outdir, paste0(basename_stub, "_PCA.pdf")))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = Timepoint)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), color = "black", size = 3, max.overlaps = Inf) +
  xlab(sprintf("PC1: %.1f%% variance", percentVar[1])) +
  ylab(sprintf("PC2: %.1f%% variance", percentVar[2])) +
  coord_fixed() +
  ggtitle(paste(title_prefix, pca_title, sep = "\n")) +
  theme_minimal()
dev.off()

# ------------------------- Heatmaps -------------------------
pdf(file.path(opt$outdir, paste0(basename_stub, "_heatmap.pdf")))
for (subset_type in c("All", "Known", "Novel")) {
  top_rows <- switch(subset_type,
                     Known = res_filtered[res_filtered$GeneType == "Known", , drop = FALSE],
                     Novel = res_filtered[res_filtered$GeneType == "Novel", , drop = FALSE],
                     res_filtered)

  if (nrow(top_rows) < 2) {
    plot.new(); title(paste(title_prefix, "Heatmap skipped - too few", subset_type, "features"), cex.main = 1.2)
    next
  }

  # prefer significant rows if available; otherwise fall back to best padj
  cand <- top_rows[!is.na(top_rows$padj) & !is.na(top_rows$log2FoldChange), , drop = FALSE]
  sig_cand <- cand[(cand$padj < plot_padj) & (abs(cand$log2FoldChange) > plot_logfc), , drop = FALSE]
  pick_from <- if (nrow(sig_cand) > 0) sig_cand else cand

  if (nrow(pick_from) < 2) {
    plot.new(); title(paste(title_prefix, "Heatmap skipped - too few usable", subset_type, "features"), cex.main = 1.2)
    next
  }

  pick_from <- pick_from[order(pick_from$padj), , drop = FALSE]
  top_ids <- rownames(head(pick_from, opt$topn))

  # intersection with VST matrix
  top_ids <- intersect(top_ids, rownames(assay(vsd)))
  if (length(top_ids) < 2) {
    plot.new(); title(paste(title_prefix, "Heatmap skipped - too few overlapping features:", subset_type), cex.main = 1.2)
    next
  }


  mat <- assay(vsd)[top_ids, , drop = FALSE]
  mat <- mat - rowMeans(mat)

  # Prefer GeneName; else right side after '|'; else id
  ids   <- rownames(mat)
  right <- sub(".*\\|", "", ids)
  rn <- right
  if ("GeneName" %in% colnames(res)) {
    # build a name map from both res and res_filtered
    map_full <- setNames(res$GeneName, rownames(res))
    map_filt <- if (exists("res_filtered")) setNames(res_filtered$GeneName, rownames(res_filtered)) else NULL
    name_map <- c(map_full, map_filt)
    rn2 <- name_map[ids]
    rn[!is.na(rn2) & rn2 != ""] <- rn2[!is.na(rn2) & rn2 != ""]
  }
  rn[is.na(rn) | rn == ""] <- ids[is.na(rn) | rn == ""]
  rownames(mat) <- rn

  annotation <- if (has_time) {
    ordering <- order(coldata$condition, coldata$time_numeric)
    mat <- mat[, ordering, drop = FALSE]
    coldata[ordering, c("condition", "Timepoint"), drop = FALSE]
  } else {
    coldata[, "condition", drop = FALSE]
  }

  heatmap_n <- nrow(top_rows)
  heatmap_main <- sprintf("Heatmap (%s - top %d/%d %s)", subset_type, length(top_ids), heatmap_n, opt$label)

  pheatmap(mat,
           annotation_col = annotation,
           annotation_colors = annotation_colors,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "ward.D2",
           color = colorRampPalette(c("green", "black", "red"))(100),
           fontsize_row = 8,
           fontsize_col = 9,
           main = paste(title_prefix, heatmap_main, sep = "\n"))
}
dev.off()
