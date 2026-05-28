#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(optparse)
})

option_list <- list(
  make_option("--counts", type = "character", help = "Contrast count matrix TSV"),
  make_option("--coldata", type = "character", help = "Contrast sample metadata TSV"),
  make_option("--metadata", type = "character", default = "", help = "Feature metadata TSV"),
  make_option(
    "--feature-id-column",
    type = "character",
    dest = "feature_id_column",
    default = "Geneid",
    help = "Feature identifier column in counts and metadata TSVs"
  ),
  make_option("--results", type = "character", help = "Full DESeq2 result TSV"),
  make_option("--filtered", type = "character", help = "Significant result TSV"),
  make_option(
    "--normalized-counts",
    type = "character",
    dest = "normalized_counts",
    help = "Normalized counts TSV"
  ),
  make_option("--summary", type = "character", help = "Summary TSV"),
  make_option("--condition-col", type = "character", dest = "condition_col", default = "condition"),
  make_option("--control-label", type = "character", dest = "control_label", default = "control"),
  make_option("--test-label", type = "character", dest = "test_label", help = "Condition to compare to control"),
  make_option("--design-formula", type = "character", dest = "design_formula", default = ""),
  make_option("--padj", type = "double", default = 0.1),
  make_option("--log2fc", type = "double", default = 1.0),
  make_option("--min-count", type = "integer", dest = "min_count", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("counts", "coldata", "results", "filtered", "normalized_counts", "summary", "test_label")
missing <- required[vapply(required, function(name) is.null(opt[[name]]) || opt[[name]] == "", logical(1))]
if (length(missing)) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

counts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"

coldata <- read.table(opt$coldata, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
coldata <- coldata[colnames(counts), , drop = FALSE]
if (!all(rownames(coldata) == colnames(counts))) {
  stop("coldata rownames do not match count matrix columns")
}
if (!(opt$condition_col %in% colnames(coldata))) {
  stop("Condition column not found in coldata: ", opt$condition_col)
}

keep <- rowSums(counts) >= opt$min_count
counts_filtered <- counts[keep, , drop = FALSE]
if (nrow(counts_filtered) == 0) {
  stop("No features retained after min-count filtering")
}

coldata[[opt$condition_col]] <- factor(
  coldata[[opt$condition_col]],
  levels = c(opt$control_label, opt$test_label)
)
if (any(is.na(coldata[[opt$condition_col]]))) {
  stop("Unexpected condition labels after contrast subsetting")
}

design_text <- trimws(opt$design_formula)
if (design_text == "") {
  design_text <- paste("~", opt$condition_col)
}
design_formula <- as.formula(design_text)
design_variables <- all.vars(design_formula)
missing_design_variables <- setdiff(design_variables, colnames(coldata))
if (length(missing_design_variables)) {
  stop("Design formula references missing coldata column(s): ", paste(missing_design_variables, collapse = ", "))
}
if (!(opt$condition_col %in% design_variables)) {
  stop("Design formula must include condition column: ", opt$condition_col)
}
for (column in design_variables) {
  if (is.character(coldata[[column]])) {
    coldata[[column]] <- factor(coldata[[column]])
  }
}
dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = coldata, design = design_formula)
dds <- tryCatch(
  DESeq(dds),
  error = function(err) {
    message <- conditionMessage(err)
    if (!grepl("all gene-wise dispersion estimates are within 2 orders of magnitude", message, fixed = TRUE)) {
      stop(err)
    }
    dds <- estimateSizeFactors(dds)
    dds <- estimateDispersionsGeneEst(dds)
    dispersions(dds) <- mcols(dds)$dispGeneEst
    nbinomWaldTest(dds)
  }
)

res <- results(dds, contrast = c(opt$condition_col, opt$test_label, opt$control_label))
res <- as.data.frame(res)
feature_id_column <- opt$feature_id_column
res[[feature_id_column]] <- rownames(res)
res <- res[, c(feature_id_column, setdiff(colnames(res), feature_id_column)), drop = FALSE]
res <- res[order(res$padj, na.last = TRUE), , drop = FALSE]

if (opt$metadata != "" && file.exists(opt$metadata)) {
  metadata <- read.table(opt$metadata, header = TRUE, sep = "\t", check.names = FALSE)
  if (feature_id_column %in% colnames(metadata)) {
    res <- merge(metadata, res, by = feature_id_column, all.y = TRUE, sort = FALSE)
  }
}

dir.create(dirname(opt$results), recursive = TRUE, showWarnings = FALSE)
write.table(res, opt$results, sep = "\t", quote = FALSE, row.names = FALSE)

significant <- res[
  !is.na(res$padj) &
    res$padj < opt$padj &
    !is.na(res$log2FoldChange) &
    abs(res$log2FoldChange) >= opt$log2fc,
  ,
  drop = FALSE
]
write.table(significant, opt$filtered, sep = "\t", quote = FALSE, row.names = FALSE)

normalized <- as.data.frame(counts(dds, normalized = TRUE))
normalized[[feature_id_column]] <- rownames(normalized)
normalized <- normalized[, c(feature_id_column, setdiff(colnames(normalized), feature_id_column)), drop = FALSE]
write.table(normalized, opt$normalized_counts, sep = "\t", quote = FALSE, row.names = FALSE)

summary_row <- data.frame(
  status = "ok",
  condition_col = opt$condition_col,
  control_label = opt$control_label,
  test_label = opt$test_label,
  design_formula = design_text,
  feature_id_column = feature_id_column,
  n_samples = ncol(counts),
  n_features_input = nrow(counts),
  n_features_tested = nrow(counts_filtered),
  n_genes_input = nrow(counts),
  n_genes_tested = nrow(counts_filtered),
  n_significant = nrow(significant),
  padj_threshold = opt$padj,
  log2fc_threshold = opt$log2fc,
  min_count = opt$min_count
)
write.table(summary_row, opt$summary, sep = "\t", quote = FALSE, row.names = FALSE)
