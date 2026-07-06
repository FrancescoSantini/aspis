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
  make_option(
    "--shrunken-results",
    type = "character",
    dest = "shrunken_results",
    default = "",
    help = "Optional TSV containing DESeq2 results with shrinkage metadata"
  ),
  make_option(
    "--transformed-counts",
    type = "character",
    dest = "transformed_counts",
    default = "",
    help = "Optional variance-stabilized count matrix TSV"
  ),
  make_option("--summary", type = "character", help = "Summary TSV"),
  make_option("--condition-col", type = "character", dest = "condition_col", default = "condition"),
  make_option("--control-label", type = "character", dest = "control_label", default = "control"),
  make_option("--test-label", type = "character", dest = "test_label", help = "Condition to compare to control"),
  make_option("--design-formula", type = "character", dest = "design_formula", default = ""),
  make_option(
    "--lfc-shrinkage",
    type = "character",
    dest = "lfc_shrinkage",
    default = "none",
    help = "LFC shrinkage method: none, normal, apeglm, ashr, or auto"
  ),
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

sanitize_coef_component <- function(value) {
  make.names(value)
}

condition_coef_name <- function(dds, condition_col, test_label, control_label) {
  names <- resultsNames(dds)
  candidates <- unique(c(
    paste0(condition_col, "_", test_label, "_vs_", control_label),
    paste0(sanitize_coef_component(condition_col), "_", sanitize_coef_component(test_label), "_vs_", sanitize_coef_component(control_label)),
    paste0(condition_col, test_label),
    paste0(sanitize_coef_component(condition_col), sanitize_coef_component(test_label))
  ))
  matched <- candidates[candidates %in% names]
  if (length(matched)) {
    return(matched[[1]])
  }
  pattern <- paste0("^", gsub("([\\W])", "\\\\\\1", sanitize_coef_component(condition_col)), ".*", gsub("([\\W])", "\\\\\\1", sanitize_coef_component(test_label)))
  regex_match <- grep(pattern, names, value = TRUE)
  regex_match <- regex_match[grepl("vs", regex_match)]
  if (length(regex_match)) {
    return(regex_match[[1]])
  }
  ""
}

add_feature_id <- function(result, feature_id_column) {
  result <- as.data.frame(result)
  result[[feature_id_column]] <- rownames(result)
  result[, c(feature_id_column, setdiff(colnames(result), feature_id_column)), drop = FALSE]
}

merge_metadata <- function(result, metadata_path, feature_id_column) {
  if (metadata_path != "" && file.exists(metadata_path)) {
    metadata <- read.table(metadata_path, header = TRUE, sep = "\t", check.names = FALSE)
    if (feature_id_column %in% colnames(metadata)) {
      result <- merge(metadata, result, by = feature_id_column, all.y = TRUE, sort = FALSE)
    }
  }
  result
}

clean_display_value <- function(value) {
  cleaned <- trimws(as.character(value))
  cleaned[is.na(cleaned)] <- ""
  missing <- tolower(cleaned) %in% c("", "na", "n/a", "none", "null", "nan", ".")
  cleaned[missing] <- ""
  wrapped <- grepl("^(NA|N/A|None|NULL|NaN)\\s*\\([^)]+\\)$", cleaned, ignore.case = TRUE)
  cleaned[wrapped] <- sub("^(NA|N/A|None|NULL|NaN)\\s*\\(([^)]+)\\)$", "\\2", cleaned[wrapped], ignore.case = TRUE)
  cleaned
}

first_clean_column <- function(result, columns) {
  for (column in columns) {
    if (column %in% colnames(result)) {
      cleaned <- clean_display_value(result[[column]])
      if (any(cleaned != "")) {
        return(cleaned)
      }
    }
  }
  rep("", nrow(result))
}

gene_display_label <- function(gene_id, gene_name) {
  gene_id <- clean_display_value(gene_id)
  gene_name <- clean_display_value(gene_name)
  ifelse(
    gene_name != "" & gene_id != "" & gene_name != gene_id,
    paste0(gene_name, " (", gene_id, ")"),
    ifelse(gene_name != "", gene_name, gene_id)
  )
}

transcript_display_label <- function(transcript_id, gene_id, gene_name) {
  transcript_id <- clean_display_value(transcript_id)
  gene_label <- gene_display_label(gene_id, gene_name)
  ifelse(
    gene_label != "" & transcript_id != "",
    paste0(gene_label, " | ", transcript_id),
    ifelse(transcript_id != "", transcript_id, gene_label)
  )
}

add_display_columns <- function(result, feature_id_column) {
  feature_id <- if (feature_id_column %in% colnames(result)) {
    clean_display_value(result[[feature_id_column]])
  } else {
    rep("", nrow(result))
  }
  gene_id <- first_clean_column(result, c("gene_id", "Geneid", "gene"))
  if (feature_id_column == "Geneid") {
    gene_id <- ifelse(gene_id != "", gene_id, feature_id)
  }
  gene_name <- first_clean_column(result, c("gene_name", "GeneName", "gene_symbol", "symbol"))
  gene_display <- first_clean_column(result, c("gene_display"))
  gene_display <- ifelse(gene_display != "", gene_display, gene_display_label(gene_id, gene_name))

  transcript_id <- first_clean_column(result, c("transcript_id", "isoform_id"))
  if (feature_id_column %in% c("transcript_id", "isoform_id")) {
    transcript_id <- ifelse(transcript_id != "", transcript_id, feature_id)
  }
  transcript_display <- first_clean_column(result, c("transcript_display"))
  transcript_display <- ifelse(
    transcript_display != "",
    transcript_display,
    ifelse(transcript_id != "", transcript_display_label(transcript_id, gene_id, gene_name), "")
  )
  existing_feature_display <- first_clean_column(result, c("feature_display"))
  computed_feature_display <- ifelse(
    transcript_display != "",
    transcript_display,
    ifelse(gene_display != "", gene_display, feature_id)
  )

  result$gene_display <- gene_display
  result$feature_display <- ifelse(
    existing_feature_display != "",
    existing_feature_display,
    computed_feature_display
  )
  preferred <- c(feature_id_column, "feature_display", "gene_display")
  result[, c(preferred[preferred %in% colnames(result)], setdiff(colnames(result), preferred)), drop = FALSE]
}

write_feature_matrix <- function(matrix, feature_id_column, path) {
  output <- as.data.frame(matrix, check.names = FALSE)
  output[[feature_id_column]] <- rownames(output)
  output <- output[, c(feature_id_column, setdiff(colnames(output), feature_id_column)), drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

run_shrinkage <- function(dds, contrast_vector, coef_name, requested_method) {
  requested <- tolower(trimws(requested_method))
  if (requested == "" || requested == "false") {
    requested <- "none"
  }
  if (requested == "none") {
    return(list(result = NULL, method = "none", reason = "not requested"))
  }
  methods <- switch(
    requested,
    auto = c("apeglm", "normal"),
    apeglm = "apeglm",
    normal = "normal",
    ashr = "ashr",
    stop("Unsupported --lfc-shrinkage method: ", requested_method)
  )
  reasons <- character()
  for (method in methods) {
    if (method == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)) {
      reasons <- c(reasons, "apeglm package is not installed")
      next
    }
    if (method == "ashr" && !requireNamespace("ashr", quietly = TRUE)) {
      reasons <- c(reasons, "ashr package is not installed")
      next
    }
    result <- tryCatch(
      {
        if (method %in% c("apeglm", "ashr") && coef_name != "") {
          lfcShrink(dds, coef = coef_name, type = method)
        } else if (method == "normal") {
          lfcShrink(dds, contrast = contrast_vector, type = method)
        } else {
          stop("coefficient name is unavailable for ", method)
        }
      },
      error = function(err) {
        reasons <<- c(reasons, paste(method, conditionMessage(err), sep = ": "))
        NULL
      }
    )
    if (!is.null(result)) {
      return(list(result = result, method = method, reason = ""))
    }
  }
  list(result = NULL, method = "none", reason = paste(reasons, collapse = "; "))
}

make_transformed_counts <- function(dds, feature_id_column, path) {
  if (path == "") {
    return(list(method = "", reason = "not requested"))
  }
  reason <- ""
  transformed <- tryCatch(
    {
      assay(varianceStabilizingTransformation(dds, blind = FALSE))
    },
    error = function(err) {
      reason <<- conditionMessage(err)
      normalized <- counts(dds, normalized = TRUE)
      log2(normalized + 1)
    }
  )
  method <- if (reason == "") "varianceStabilizingTransformation" else "log2_normalized_fallback"
  write_feature_matrix(transformed, feature_id_column, path)
  list(method = method, reason = reason)
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

contrast_vector <- c(opt$condition_col, opt$test_label, opt$control_label)
raw_res <- results(dds, contrast = contrast_vector)
feature_id_column <- opt$feature_id_column
coef_name <- condition_coef_name(dds, opt$condition_col, opt$test_label, opt$control_label)
shrinkage <- run_shrinkage(dds, contrast_vector, coef_name, opt$lfc_shrinkage)

raw_table <- add_feature_id(raw_res, feature_id_column)
res <- raw_table
res$raw_log2FoldChange <- raw_table$log2FoldChange
if (!is.null(shrinkage$result)) {
  shrunken_table <- add_feature_id(shrinkage$result, feature_id_column)
  shrunken_lfc <- shrunken_table$log2FoldChange[match(res[[feature_id_column]], shrunken_table[[feature_id_column]])]
  res$shrunken_log2FoldChange <- shrunken_lfc
  res$log2FoldChange <- shrunken_lfc
} else {
  res$shrunken_log2FoldChange <- NA_real_
}
res$lfc_shrinkage_method <- shrinkage$method
res <- res[order(res$padj, na.last = TRUE), , drop = FALSE]

res <- merge_metadata(res, opt$metadata, feature_id_column)
res <- add_display_columns(res, feature_id_column)

dir.create(dirname(opt$results), recursive = TRUE, showWarnings = FALSE)
write.table(res, opt$results, sep = "\t", quote = FALSE, row.names = FALSE)

if (opt$shrunken_results != "") {
  write.table(res, opt$shrunken_results, sep = "\t", quote = FALSE, row.names = FALSE)
}

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
write_feature_matrix(normalized, feature_id_column, opt$normalized_counts)
transform_info <- make_transformed_counts(dds, feature_id_column, opt$transformed_counts)

summary_row <- data.frame(
  status = "ok",
  condition_col = opt$condition_col,
  control_label = opt$control_label,
  test_label = opt$test_label,
  design_formula = design_text,
  contrast = paste(contrast_vector, collapse = ","),
  coefficient = coef_name,
  feature_id_column = feature_id_column,
  n_samples = ncol(counts),
  n_features_input = nrow(counts),
  n_features_tested = nrow(counts_filtered),
  n_genes_input = nrow(counts),
  n_genes_tested = nrow(counts_filtered),
  n_significant = nrow(significant),
  lfc_shrinkage_requested = opt$lfc_shrinkage,
  lfc_shrinkage_method = shrinkage$method,
  lfc_shrinkage_reason = shrinkage$reason,
  transformed_counts_method = transform_info$method,
  transformed_counts_reason = transform_info$reason,
  padj_threshold = opt$padj,
  log2fc_threshold = opt$log2fc,
  min_count = opt$min_count
)
write.table(summary_row, opt$summary, sep = "\t", quote = FALSE, row.names = FALSE)
