#!/usr/bin/env Rscript
# Native DRIMSeq DTU runner for ASPIS contrast-level transcript count tables.

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

blocked <- function(message, summary_path, gene_results, transcript_results) {
  write_empty(gene_results, c("gene_id", "pvalue", "padj", "status"))
  write_empty(transcript_results, c("gene_id", "feature_id", "pvalue", "padj", "status"))
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
transcript_results <- need(opts, "transcript-results")
summary_path <- need(opts, "summary")
condition_col <- ifelse(is.null(opts[["condition-col"]]), "condition", opts[["condition-col"]])
control_label <- ifelse(is.null(opts[["control-label"]]), "control", opts[["control-label"]])
test_label <- ifelse(is.null(opts[["test-label"]]), "treated", opts[["test-label"]])
min_count <- as.integer(ifelse(is.null(opts[["min-count"]]), "10", opts[["min-count"]]))
min_samples <- as.integer(ifelse(is.null(opts[["min-samples"]]), "2", opts[["min-samples"]]))
min_proportion <- as.numeric(ifelse(is.null(opts[["min-proportion"]]), "0.05", opts[["min-proportion"]]))
min_gene_count <- as.integer(ifelse(is.null(opts[["min-gene-count"]]), "10", opts[["min-gene-count"]]))
min_transcripts_per_gene <- as.integer(ifelse(is.null(opts[["min-transcripts-per-gene"]]), "2", opts[["min-transcripts-per-gene"]]))

if (!requireNamespace("DRIMSeq", quietly = TRUE)) {
  blocked("R package DRIMSeq is not installed in this environment", summary_path, gene_results, transcript_results)
}

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
coldata <- read.delim(coldata_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!"sample_id" %in% names(coldata)) {
  blocked("DRIMSeq coldata is missing sample_id", summary_path, gene_results, transcript_results)
}
if (!condition_col %in% names(coldata)) {
  blocked(sprintf("DRIMSeq coldata is missing condition column %s", condition_col), summary_path, gene_results, transcript_results)
}
if (ncol(counts) < 3 || nrow(counts) == 0) {
  blocked("DRIMSeq count matrix has no transcript rows or too few sample columns", summary_path, gene_results, transcript_results)
}

transcript_col <- names(counts)[[1]]
metadata_transcript_col <- intersect(c("transcript_id", "feature_id", "isoform_id", "id"), names(metadata))
metadata_gene_col <- intersect(c("gene_id", "group_id", "gene"), names(metadata))
if (length(metadata_transcript_col) == 0 || length(metadata_gene_col) == 0) {
  blocked("transcript metadata must contain transcript_id/feature_id and gene_id/group_id columns", summary_path, gene_results, transcript_results)
}
metadata_transcript_col <- metadata_transcript_col[[1]]
metadata_gene_col <- metadata_gene_col[[1]]

sample_ids <- coldata$sample_id
missing_samples <- setdiff(sample_ids, names(counts))
if (length(missing_samples) > 0) {
  blocked(sprintf("sample columns missing from count matrix: %s", paste(missing_samples, collapse = ",")), summary_path, gene_results, transcript_results)
}

coldata <- coldata[coldata$sample_id %in% names(counts), , drop = FALSE]
coldata[[condition_col]] <- factor(coldata[[condition_col]], levels = c(control_label, test_label))
coldata <- coldata[!is.na(coldata[[condition_col]]), , drop = FALSE]
if (length(unique(coldata[[condition_col]])) != 2) {
  blocked("DRIMSeq requires exactly two condition labels after contrast filtering", summary_path, gene_results, transcript_results)
}
if (any(table(coldata[[condition_col]]) < min_samples)) {
  blocked("DRIMSeq contrast has fewer samples per group than the configured min_samples", summary_path, gene_results, transcript_results)
}

merged <- merge(
  counts,
  metadata[, c(metadata_transcript_col, metadata_gene_col), drop = FALSE],
  by.x = transcript_col,
  by.y = metadata_transcript_col,
  all.x = FALSE,
  all.y = FALSE
)
if (nrow(merged) == 0) {
  blocked("no count rows matched transcript metadata", summary_path, gene_results, transcript_results)
}

sample_cols <- coldata$sample_id
count_values <- as.data.frame(lapply(merged[, sample_cols, drop = FALSE], as.numeric), check.names = FALSE)
feature_keep <- rowSums(count_values >= min_count, na.rm = TRUE) >= min_samples
merged <- merged[feature_keep, , drop = FALSE]
count_values <- count_values[feature_keep, , drop = FALSE]
if (nrow(merged) == 0) {
  blocked("no transcripts remained after DRIMSeq feature count filtering", summary_path, gene_results, transcript_results)
}

gene_col <- metadata_gene_col
feature_col <- transcript_col
gene_total <- rowsum(as.matrix(count_values), group = merged[[gene_col]], reorder = FALSE)
expressed_genes <- rownames(gene_total)[rowSums(gene_total >= min_gene_count, na.rm = TRUE) >= min_samples]
merged <- merged[merged[[gene_col]] %in% expressed_genes, , drop = FALSE]
if (nrow(merged) == 0) {
  blocked("no genes remained after DRIMSeq gene count filtering", summary_path, gene_results, transcript_results)
}

gene_sizes <- table(merged[[gene_col]])
multi_genes <- names(gene_sizes)[gene_sizes >= min_transcripts_per_gene]
merged <- merged[merged[[gene_col]] %in% multi_genes, , drop = FALSE]
if (nrow(merged) == 0) {
  blocked("no multi-isoform genes remained after DRIMSeq filtering", summary_path, gene_results, transcript_results)
}

dm_counts <- data.frame(
  gene_id = merged[[gene_col]],
  feature_id = merged[[feature_col]],
  merged[, sample_cols, drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
dm_samples <- data.frame(
  sample_id = coldata$sample_id,
  condition = coldata[[condition_col]],
  stringsAsFactors = FALSE
)

tryCatch({
  d <- DRIMSeq::dmDSdata(counts = dm_counts, samples = dm_samples)
  d <- DRIMSeq::dmFilter(
    d,
    min_samps_feature_expr = min_samples,
    min_feature_expr = min_count,
    min_samps_feature_prop = min_samples,
    min_feature_prop = min_proportion,
    min_samps_gene_expr = min_samples,
    min_gene_expr = min_gene_count
  )
  if (nrow(DRIMSeq::counts(d)) == 0) {
    blocked("DRIMSeq internal filtering removed all features", summary_path, gene_results, transcript_results)
  }
  design <- model.matrix(~ condition, data = DRIMSeq::samples(d))
  d <- DRIMSeq::dmPrecision(d, design = design)
  d <- DRIMSeq::dmFit(d, design = design)
  d <- DRIMSeq::dmTest(d, coef = 2)
  result <- as.data.frame(DRIMSeq::results(d), stringsAsFactors = FALSE)
}, error = function(e) {
  stop(e)
})

if (!"gene_id" %in% names(result)) {
  blocked("DRIMSeq results did not include gene_id", summary_path, gene_results, transcript_results)
}
if (!"feature_id" %in% names(result) && "feature" %in% names(result)) {
  names(result)[names(result) == "feature"] <- "feature_id"
}
if (!"pvalue" %in% names(result) && "p.value" %in% names(result)) {
  names(result)[names(result) == "p.value"] <- "pvalue"
}
if (!"padj" %in% names(result) && "adj_pvalue" %in% names(result)) {
  names(result)[names(result) == "adj_pvalue"] <- "padj"
}
if (!"padj" %in% names(result) && "adj.p.value" %in% names(result)) {
  names(result)[names(result) == "adj.p.value"] <- "padj"
}
if (!"padj" %in% names(result) && "pvalue" %in% names(result)) {
  result$padj <- ave(result$pvalue, result$gene_id, FUN = function(x) p.adjust(x, method = "BH"))
}
result$status <- "ok"
dir.create(dirname(gene_results), recursive = TRUE, showWarnings = FALSE)
write.table(result, file = gene_results, sep = "\t", quote = FALSE, row.names = FALSE)

usage_counts <- as.data.frame(lapply(merged[, sample_cols, drop = FALSE], as.numeric), check.names = FALSE)
rownames(usage_counts) <- seq_len(nrow(usage_counts))
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
transcript_usage <- data.frame(
  gene_id = merged[[gene_col]],
  feature_id = merged[[feature_col]],
  mean_usage_control = mean_or_na(usage_props, control_samples),
  mean_usage_test = mean_or_na(usage_props, test_samples),
  delta_usage = mean_or_na(usage_props, test_samples) - mean_or_na(usage_props, control_samples),
  mean_count_control = mean_or_na(usage_counts, control_samples),
  mean_count_test = mean_or_na(usage_counts, test_samples),
  status = "ok",
  stringsAsFactors = FALSE
)
write.table(transcript_usage, file = transcript_results, sep = "\t", quote = FALSE, row.names = FALSE)
summary <- data.frame(
  status = "ok",
  reason = "",
  n_input_transcripts = nrow(counts),
  n_tested_genes = nrow(result),
  n_usage_transcripts = nrow(transcript_usage),
  control_label = control_label,
  test_label = test_label,
  stringsAsFactors = FALSE
)
write.table(summary, file = summary_path, sep = "\t", quote = FALSE, row.names = FALSE)