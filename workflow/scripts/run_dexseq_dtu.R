#!/usr/bin/env Rscript
# Native DEXSeq-style transcript-feature usage runner for ASPIS contrast-level tables.

blocked_exit <- 20

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument: %s", key))
    }
    name <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

need <- function(opts, key) {
  value <- opts[[key]]
  if (is.null(value) || identical(value, "")) {
    stop(sprintf("Missing required argument --%s", key))
  }
  value
}

write_empty <- function(path, columns) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns)),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

blocked <- function(message, summary_path, gene_results, feature_results) {
  write_empty(gene_results, c("gene_id", "n_features", "min_pvalue", "min_padj", "status"))
  write_empty(feature_results, c("gene_id", "feature_id", "pvalue", "padj", "status"))
  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
  summary <- data.frame(status = "blocked", reason = message, stringsAsFactors = FALSE)
  write.table(summary, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message(message)
  quit(status = blocked_exit)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
counts_path <- need(opts, "counts")
coldata_path <- need(opts, "coldata")
metadata_path <- need(opts, "metadata")
gene_results <- need(opts, "gene-results")
feature_results <- need(opts, "feature-results")
summary_path <- need(opts, "summary")
condition_col <- ifelse(is.null(opts[["condition-col"]]), "condition", opts[["condition-col"]])
control_label <- ifelse(is.null(opts[["control-label"]]), "control", opts[["control-label"]])
test_label <- ifelse(is.null(opts[["test-label"]]), "treated", opts[["test-label"]])
min_count <- as.integer(ifelse(is.null(opts[["min-count"]]), "10", opts[["min-count"]]))
min_samples <- as.integer(ifelse(is.null(opts[["min-samples"]]), "2", opts[["min-samples"]]))
min_gene_count <- as.integer(ifelse(is.null(opts[["min-gene-count"]]), "10", opts[["min-gene-count"]]))
min_transcripts_per_gene <- as.integer(ifelse(is.null(opts[["min-transcripts-per-gene"]]), "2", opts[["min-transcripts-per-gene"]]))

if (!requireNamespace("DEXSeq", quietly = TRUE)) {
  blocked("R package DEXSeq is not installed in this environment", summary_path, gene_results, feature_results)
}

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
coldata <- read.delim(coldata_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!"sample_id" %in% names(coldata)) {
  blocked("DEXSeq coldata is missing sample_id", summary_path, gene_results, feature_results)
}
if (!condition_col %in% names(coldata)) {
  blocked(sprintf("DEXSeq coldata is missing condition column %s", condition_col), summary_path, gene_results, feature_results)
}
if (ncol(counts) < 3 || nrow(counts) == 0) {
  blocked("DEXSeq count matrix has no transcript rows or too few sample columns", summary_path, gene_results, feature_results)
}

feature_col <- names(counts)[[1]]
metadata_feature_col <- intersect(c("transcript_id", "feature_id", "isoform_id", "id"), names(metadata))
metadata_gene_col <- intersect(c("gene_id", "group_id", "gene"), names(metadata))
if (length(metadata_feature_col) == 0 || length(metadata_gene_col) == 0) {
  blocked("transcript metadata must contain transcript_id/feature_id and gene_id/group_id columns", summary_path, gene_results, feature_results)
}
metadata_feature_col <- metadata_feature_col[[1]]
metadata_gene_col <- metadata_gene_col[[1]]

sample_ids <- coldata$sample_id
missing_samples <- setdiff(sample_ids, names(counts))
if (length(missing_samples) > 0) {
  blocked(sprintf("sample columns missing from count matrix: %s", paste(missing_samples, collapse = ",")), summary_path, gene_results, feature_results)
}

coldata <- coldata[coldata$sample_id %in% names(counts), , drop = FALSE]
coldata[[condition_col]] <- factor(coldata[[condition_col]], levels = c(control_label, test_label))
coldata <- coldata[!is.na(coldata[[condition_col]]), , drop = FALSE]
if (length(unique(coldata[[condition_col]])) != 2) {
  blocked("DEXSeq requires exactly two condition labels after contrast filtering", summary_path, gene_results, feature_results)
}
if (any(table(coldata[[condition_col]]) < min_samples)) {
  blocked("DEXSeq contrast has fewer samples per group than the configured min_samples", summary_path, gene_results, feature_results)
}

merged <- merge(
  counts,
  metadata[, c(metadata_feature_col, metadata_gene_col), drop = FALSE],
  by.x = feature_col,
  by.y = metadata_feature_col,
  all.x = FALSE,
  all.y = FALSE
)
if (nrow(merged) == 0) {
  blocked("no count rows matched transcript metadata", summary_path, gene_results, feature_results)
}

sample_cols <- coldata$sample_id
count_values <- as.data.frame(lapply(merged[, sample_cols, drop = FALSE], as.numeric), check.names = FALSE)
feature_keep <- rowSums(count_values >= min_count, na.rm = TRUE) >= min_samples
merged <- merged[feature_keep, , drop = FALSE]
count_values <- count_values[feature_keep, , drop = FALSE]
if (nrow(merged) == 0) {
  blocked("no transcript features remained after DEXSeq feature count filtering", summary_path, gene_results, feature_results)
}

gene_col <- metadata_gene_col
gene_total <- rowsum(as.matrix(count_values), group = merged[[gene_col]], reorder = FALSE)
expressed_genes <- rownames(gene_total)[rowSums(gene_total >= min_gene_count, na.rm = TRUE) >= min_samples]
keep_gene <- merged[[gene_col]] %in% expressed_genes
merged <- merged[keep_gene, , drop = FALSE]
count_values <- count_values[keep_gene, , drop = FALSE]
if (nrow(merged) == 0) {
  blocked("no genes remained after DEXSeq gene count filtering", summary_path, gene_results, feature_results)
}

gene_sizes <- table(merged[[gene_col]])
multi_genes <- names(gene_sizes)[gene_sizes >= min_transcripts_per_gene]
keep_multi <- merged[[gene_col]] %in% multi_genes
merged <- merged[keep_multi, , drop = FALSE]
count_values <- count_values[keep_multi, , drop = FALSE]
if (nrow(merged) == 0) {
  blocked("no multi-feature genes remained after DEXSeq filtering", summary_path, gene_results, feature_results)
}

sample_data <- data.frame(
  sample = coldata$sample_id,
  condition = coldata[[condition_col]],
  stringsAsFactors = FALSE
)
rownames(sample_data) <- sample_data$sample
count_matrix <- as.matrix(merged[, sample_cols, drop = FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- make.unique(as.character(merged[[feature_col]]))

tryCatch({
  dxd <- DEXSeq::DEXSeqDataSetFromMatrix(
    countData = count_matrix,
    sampleData = sample_data,
    design = ~ sample + exon + condition:exon,
    featureID = as.character(merged[[feature_col]]),
    groupID = as.character(merged[[gene_col]])
  )
  dxd <- DEXSeq::estimateSizeFactors(dxd)
  dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE)
  dxd <- DEXSeq::testForDEU(dxd, reducedModel = ~ sample + exon, fullModel = ~ sample + exon + condition:exon)
  dxd <- DEXSeq::estimateExonFoldChanges(dxd, fitExpToVar = "condition")
  result <- as.data.frame(DEXSeq::DEXSeqResults(dxd), stringsAsFactors = FALSE)
}, error = function(e) {
  stop(e)
})

if (!"groupID" %in% names(result)) {
  blocked("DEXSeq results did not include groupID", summary_path, gene_results, feature_results)
}
if (!"featureID" %in% names(result)) {
  blocked("DEXSeq results did not include featureID", summary_path, gene_results, feature_results)
}
if (!"pvalue" %in% names(result)) {
  blocked("DEXSeq results did not include pvalue", summary_path, gene_results, feature_results)
}
if (!"padj" %in% names(result)) {
  result$padj <- p.adjust(result$pvalue, method = "BH")
}

log2fc_col <- intersect(c("log2fold_test_control", "log2fold_condition_treated_control", "log2fold_condition"), names(result))
log2fc <- if (length(log2fc_col) > 0) result[[log2fc_col[[1]]]] else NA_real_
feature_table <- data.frame(
  gene_id = result$groupID,
  feature_id = result$featureID,
  statistic = if ("stat" %in% names(result)) result$stat else NA_real_,
  log2_fold_change = log2fc,
  pvalue = result$pvalue,
  padj = result$padj,
  event_type = "transcript_feature_usage",
  status = "ok",
  stringsAsFactors = FALSE
)

usage_counts <- as.data.frame(lapply(merged[, sample_cols, drop = FALSE], as.numeric), check.names = FALSE)
gene_totals <- usage_counts
for (sample in sample_cols) {
  gene_totals[[sample]] <- ave(usage_counts[[sample]], merged[[gene_col]], FUN = sum)
}
usage_props <- usage_counts
for (sample in sample_cols) {
  usage_props[[sample]] <- ifelse(gene_totals[[sample]] > 0, usage_counts[[sample]] / gene_totals[[sample]], NA_real_)
}
control_samples <- coldata$sample_id[coldata[[condition_col]] == control_label]
test_samples <- coldata$sample_id[coldata[[condition_col]] == test_label]
mean_or_na <- function(frame, cols) {
  if (length(cols) == 0) {
    return(rep(NA_real_, nrow(frame)))
  }
  rowMeans(frame[, cols, drop = FALSE], na.rm = TRUE)
}
usage_table <- data.frame(
  gene_id = merged[[gene_col]],
  feature_id = merged[[feature_col]],
  mean_usage_control = mean_or_na(usage_props, control_samples),
  mean_usage_test = mean_or_na(usage_props, test_samples),
  delta_usage = mean_or_na(usage_props, test_samples) - mean_or_na(usage_props, control_samples),
  stringsAsFactors = FALSE
)
feature_table <- merge(feature_table, usage_table, by = c("gene_id", "feature_id"), all.x = TRUE)

dir.create(dirname(feature_results), recursive = TRUE, showWarnings = FALSE)
write.table(feature_table, file = feature_results, sep = "\t", quote = FALSE, row.names = FALSE)

gene_summary <- aggregate(
  cbind(pvalue, padj) ~ gene_id,
  data = feature_table,
  FUN = function(x) suppressWarnings(min(x, na.rm = TRUE))
)
names(gene_summary) <- c("gene_id", "min_pvalue", "min_padj")
feature_counts <- aggregate(feature_id ~ gene_id, data = feature_table, FUN = length)
names(feature_counts) <- c("gene_id", "n_features")
gene_summary <- merge(gene_summary, feature_counts, by = "gene_id", all.x = TRUE)
gene_summary$status <- "ok"
write.table(gene_summary, file = gene_results, sep = "\t", quote = FALSE, row.names = FALSE)

summary <- data.frame(
  status = "ok",
  reason = "DEXSeq run on transcript features grouped by gene; exon-bin DEXSeq requires separate exon-count inputs",
  n_input_transcripts = nrow(counts),
  n_tested_genes = length(unique(feature_table$gene_id)),
  n_usage_transcripts = nrow(feature_table),
  control_label = control_label,
  test_label = test_label,
  stringsAsFactors = FALSE
)
write.table(summary, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)
