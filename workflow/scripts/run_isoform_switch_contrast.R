#!/usr/bin/env Rscript

parse_args <- function(argv) {
  out <- list(
    gene_expr = "1",
    isoform_expr = "1",
    padj = "0.1",
    dif = "0.1",
    max_genes = "30",
    genome_object = ""
  )
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (i == length(argv)) {
      stop("Missing value for argument: ", key)
    }
    name <- chartr("-", "_", substring(key, 3))
    out[[name]] <- argv[[i + 1]]
    i <- i + 2
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || value == "") {
    stop("Missing required argument: --", chartr("_", "-", name))
  }
  value
}

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input does not exist: ", path)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(data, path) {
  ensure_parent(path)
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

blank_pdf <- function(path, title, message) {
  ensure_parent(path)
  grDevices::pdf(path)
  plot.new()
  title(main = title)
  text(0.5, 0.5, message)
  grDevices::dev.off()
}

write_empty_outputs <- function(args, reason) {
  empty <- data.frame(reason = reason, check.names = FALSE)
  write_tsv(empty, args[["results"]])
  write_tsv(empty, args[["summary"]])
  write_tsv(empty, args[["consequences"]])
  write_tsv(empty, args[["detailed"]])
  writeLines(reason, args[["expression_summary"]])
  blank_pdf(args[["qc_pdf"]], args[["contrast_id"]], reason)
  blank_pdf(args[["dif_distribution_pdf"]], args[["contrast_id"]], reason)
  ensure_parent(args[["switch_rds"]])
  saveRDS(list(reason = reason), file = args[["switch_rds"]])
  ensure_parent(args[["nt_fasta"]])
  writeLines(character(), args[["nt_fasta"]])
  ensure_parent(args[["aa_fasta"]])
  writeLines(character(), args[["aa_fasta"]])
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required R package is not installed: ", package)
  }
}

resolve_genome_object <- function(spec) {
  if (is.null(spec) || spec == "") {
    return(NULL)
  }
  parts <- strsplit(spec, "::", fixed = TRUE)[[1]]
  if (length(parts) != 2 || parts[[1]] == "" || parts[[2]] == "") {
    stop("--genome-object must use package::object syntax")
  }
  require_package(parts[[1]])
  get(parts[[2]], envir = asNamespace(parts[[1]]))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required <- c(
    "counts", "design", "gtf", "results", "summary", "qc_pdf", "switch_rds",
    "consequences", "detailed", "dif_distribution_pdf", "nt_fasta", "aa_fasta",
    "expression_summary", "contrast_id", "control_label", "test_label"
  )
  for (name in required) {
    required_arg(args, name)
  }

  require_package("IsoformSwitchAnalyzeR")

  gene_expr <- as.numeric(args[["gene_expr"]])
  isoform_expr <- as.numeric(args[["isoform_expr"]])
  padj_cutoff <- as.numeric(args[["padj"]])
  dif_cutoff <- as.numeric(args[["dif"]])
  max_genes <- as.integer(args[["max_genes"]])
  if (is.na(max_genes) || max_genes < 1) {
    stop("--max-genes must be a positive integer")
  }

  counts <- read_tsv(args[["counts"]])
  design <- read_tsv(args[["design"]])
  if (ncol(counts) < 3) {
    write_empty_outputs(args, "Isoform-switch analysis requires at least two samples")
    return(invisible(0))
  }
  if (!all(c("sampleID", "condition") %in% colnames(design))) {
    stop("Design table must contain sampleID and condition columns")
  }

  colnames(counts)[[1]] <- "isoform_id"
  samples <- intersect(colnames(counts)[-1], design$sampleID)
  if (length(samples) < 2) {
    write_empty_outputs(args, "No overlapping count/design samples for isoform-switch analysis")
    return(invisible(0))
  }
  counts <- counts[, c("isoform_id", samples), drop = FALSE]
  for (sample in samples) {
    counts[[sample]] <- suppressWarnings(as.numeric(counts[[sample]]))
  }
  design <- design[match(samples, design$sampleID), , drop = FALSE]
  design$condition <- factor(design$condition, levels = c(args[["control_label"]], args[["test_label"]]))
  design <- design[, c("sampleID", "condition"), drop = FALSE]

  isa <- asNamespace("IsoformSwitchAnalyzeR")
  switch_list <- isa$importRdata(
    isoformCountMatrix = counts,
    designMatrix = design,
    isoformExonAnno = args[["gtf"]],
    comparisonsToMake = NULL,
    ignoreAfterBar = FALSE,
    showProgress = TRUE
  )

  count_matrix <- as.data.frame(counts, check.names = FALSE)
  rownames(count_matrix) <- count_matrix$isoform_id
  count_matrix <- count_matrix[, samples, drop = FALSE]
  lib_sizes <- colSums(count_matrix, na.rm = TRUE)
  lib_sizes[lib_sizes == 0] <- NA
  cpm <- t(t(as.matrix(count_matrix)) / lib_sizes * 1e6)
  isoform_expression <- rowMeans(cpm, na.rm = TRUE)
  features <- as.data.frame(switch_list$isoformFeatures)
  features$isoform_expression <- isoform_expression[match(features$isoform_id, names(isoform_expression))]
  gene_expression <- tapply(features$isoform_expression, features$gene_id, sum, na.rm = TRUE)
  features$gene_expression <- gene_expression[match(features$gene_id, names(gene_expression))]
  switch_list$isoformFeatures <- features

  switch_list <- isa$preFilter(switch_list, removeSingleIsoformGenes = TRUE)
  if (nrow(switch_list$isoformFeatures) == 0) {
    write_empty_outputs(args, "No multi-isoform genes remain after pre-filtering")
    return(invisible(0))
  }
  expression_pass <- with(
    switch_list$isoformFeatures,
    !is.na(gene_expression) & !is.na(isoform_expression) &
      gene_expression >= gene_expr & isoform_expression >= isoform_expr
  )
  if (!any(expression_pass)) {
    write_empty_outputs(args, "No isoforms pass configured expression filters")
    return(invisible(0))
  }
  switch_list <- isa$subsetSwitchAnalyzeRlist(switch_list, expression_pass)

  switch_list <- isa$isoformSwitchTestDEXSeq(switch_list)

  genome <- resolve_genome_object(args[["genome_object"]])
  generated_nt <- ""
  generated_aa <- ""
  if (!is.null(genome)) {
    sequence_prefix <- file.path(dirname(args[["nt_fasta"]]), "isoformSwitchAnalyzeR")
    switch_list <- isa$analyzeORF(switch_list, orfMethod = "longest", genomeObject = genome)
    switch_list <- isa$extractSequence(switch_list, genomeObject = genome, outputPrefix = sequence_prefix)
    generated_nt <- paste0(sequence_prefix, "_isoform_nt.fasta")
    generated_aa <- paste0(sequence_prefix, "_isoform_AA.fasta")
    switch_list <- isa$analyzeIntronRetention(switch_list)
    switch_list <- isa$analyzeSwitchConsequences(
      switch_list,
      consequencesToAnalyze = c("intron_retention", "NMD_status", "ORF_seq_similarity")
    )
  }

  analysis <- as.data.frame(switch_list$isoformSwitchAnalysis)
  write_tsv(analysis, args[["results"]])

  summary_df <- tryCatch(
    as.data.frame(isa$extractSwitchSummary(switch_list)),
    error = function(err) data.frame(reason = conditionMessage(err), check.names = FALSE)
  )
  write_tsv(summary_df, args[["summary"]])

  features <- as.data.frame(switch_list$isoformFeatures)
  feature_cols <- intersect(c("isoform_id", "gene_id", "gene_name"), colnames(features))
  if (length(feature_cols) && "isoform_id" %in% colnames(analysis)) {
    detailed <- merge(analysis, features[, feature_cols, drop = FALSE], by = "isoform_id", all.x = TRUE)
  } else {
    detailed <- analysis
  }
  write_tsv(detailed, args[["detailed"]])

  consequences <- if (is.data.frame(switch_list$switchConsequence)) {
    as.data.frame(switch_list$switchConsequence)
  } else {
    data.frame()
  }
  write_tsv(consequences, args[["consequences"]])

  expression_summary <- capture.output({
    cat("contrast_id\t", args[["contrast_id"]], "\n", sep = "")
    cat("n_isoforms\t", nrow(features), "\n", sep = "")
    if ("gene_id" %in% colnames(features)) {
      cat("n_genes\t", length(unique(features$gene_id)), "\n", sep = "")
    }
    cat("gene_expr_cutoff\t", gene_expr, "\n", sep = "")
    cat("isoform_expr_cutoff\t", isoform_expr, "\n", sep = "")
  })
  ensure_parent(args[["expression_summary"]])
  writeLines(expression_summary, args[["expression_summary"]])

  if ("dIF" %in% colnames(analysis) && nrow(analysis) > 0) {
    ensure_parent(args[["dif_distribution_pdf"]])
    grDevices::pdf(args[["dif_distribution_pdf"]])
    hist(analysis$dIF, breaks = 30, col = "steelblue", border = "white", main = "Distribution of dIF values", xlab = "dIF")
    grDevices::dev.off()
  } else {
    blank_pdf(args[["dif_distribution_pdf"]], args[["contrast_id"]], "No dIF values available")
  }

  if (all(c("gene_id", "dIF", "padj") %in% colnames(detailed))) {
    candidates <- detailed[!is.na(detailed$padj) & detailed$padj < padj_cutoff & abs(detailed$dIF) >= dif_cutoff, , drop = FALSE]
    candidates <- candidates[order(abs(candidates$dIF), decreasing = TRUE), , drop = FALSE]
    candidates <- candidates[!duplicated(candidates$gene_id), , drop = FALSE]
    candidates <- head(candidates, max_genes)
  } else {
    candidates <- data.frame()
  }
  ensure_parent(args[["qc_pdf"]])
  grDevices::pdf(args[["qc_pdf"]], width = 8, height = 6)
  if (nrow(candidates) == 0) {
    plot.new()
    title("No significant isoform switches")
  } else {
    labels <- if ("gene_name" %in% colnames(candidates)) candidates$gene_name else candidates$gene_id
    barplot(candidates$dIF, names.arg = labels, las = 2, col = "steelblue", main = "Top isoform switches", ylab = "dIF")
    abline(h = c(-dif_cutoff, dif_cutoff), lty = 2, col = "grey40")
  }
  grDevices::dev.off()

  ensure_parent(args[["nt_fasta"]])
  ensure_parent(args[["aa_fasta"]])
  if (generated_nt != "" && file.exists(generated_nt)) {
    file.rename(generated_nt, args[["nt_fasta"]])
  } else {
    writeLines(character(), args[["nt_fasta"]])
  }
  if (generated_aa != "" && file.exists(generated_aa)) {
    file.rename(generated_aa, args[["aa_fasta"]])
  } else {
    writeLines(character(), args[["aa_fasta"]])
  }

  ensure_parent(args[["switch_rds"]])
  saveRDS(switch_list, file = args[["switch_rds"]])
}

main()
