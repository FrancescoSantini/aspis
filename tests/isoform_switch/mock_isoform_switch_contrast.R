#!/usr/bin/env Rscript

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (i == length(argv)) {
      stop("Missing value for argument: ", key)
    }
    out[[chartr("-", "_", substring(key, 3))]] <- argv[[i + 1]]
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

write_tsv <- function(data, path) {
  ensure_parent(path)
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

blank_pdf <- function(path, title_text) {
  ensure_parent(path)
  grDevices::pdf(path)
  plot.new()
  title(main = title_text)
  text(0.5, 0.5, "isoform-switch ready-contract smoke")
  grDevices::dev.off()
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c(
  "counts", "design", "gtf", "results", "summary", "qc_pdf", "switch_rds",
  "consequences", "detailed", "dif_distribution_pdf", "nt_fasta", "aa_fasta",
  "expression_summary", "contrast_id", "control_label", "test_label"
)
for (name in required) {
  required_arg(args, name)
}

counts <- read.delim(args[["counts"]], check.names = FALSE, stringsAsFactors = FALSE)
design <- read.delim(args[["design"]], check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(colnames(counts)[[1]], "isoform_id")) {
  stop("counts first column must be isoform_id")
}
if (!all(c("sampleID", "condition") %in% colnames(design))) {
  stop("design must contain sampleID and condition")
}
if (!file.exists(args[["gtf"]])) {
  stop("annotation GTF does not exist: ", args[["gtf"]])
}

sample_count <- ncol(counts) - 1
isoform_count <- nrow(counts)
result <- data.frame(
  isoform_id = counts[["isoform_id"]],
  gene_id = sub("_[0-9]+$", "", counts[["isoform_id"]]),
  dIF = seq(-0.25, 0.25, length.out = isoform_count),
  padj = seq(0.01, 0.2, length.out = isoform_count),
  check.names = FALSE
)
write_tsv(result, args[["results"]])
write_tsv(
  data.frame(
    status = "ok",
    contrast_id = args[["contrast_id"]],
    n_isoforms = isoform_count,
    n_samples = sample_count,
    control_label = args[["control_label"]],
    test_label = args[["test_label"]],
    check.names = FALSE
  ),
  args[["summary"]]
)
write_tsv(result, args[["detailed"]])
write_tsv(data.frame(contrast_id = args[["contrast_id"]], consequence = "mock", check.names = FALSE), args[["consequences"]])
writeLines(c("mock_status\tok", paste0("n_samples\t", sample_count)), args[["expression_summary"]])
blank_pdf(args[["qc_pdf"]], args[["contrast_id"]])
blank_pdf(args[["dif_distribution_pdf"]], paste(args[["contrast_id"]], "dIF"))
ensure_parent(args[["switch_rds"]])
saveRDS(list(counts = counts, design = design), file = args[["switch_rds"]])
ensure_parent(args[["nt_fasta"]])
writeLines(c(">mock_nt", "ATGGCC"), args[["nt_fasta"]])
ensure_parent(args[["aa_fasta"]])
writeLines(c(">mock_aa", "MA"), args[["aa_fasta"]])
cat("mock isoform-switch contrast completed\n")
