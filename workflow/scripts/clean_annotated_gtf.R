#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(dplyr)
})

option_list <- list(
  make_option("--input", type = "character", help = "Path to raw annotated.gtf"),
  make_option("--output", type = "character", help = "Path to cleaned annotated.gtf")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Load GTF
gtf <- rtracklayer::import(opt$input)
gtf_df <- as.data.frame(gtf)

# --- Patch gene_id and gene_name
fixed_gene_id <- ifelse(!is.na(gtf_df$ref_gene_id), gtf_df$ref_gene_id, gtf_df$gene_id)
fixed_gene_name <- ifelse(!is.na(gtf_df$ref_gene_name), gtf_df$ref_gene_name, gtf_df$gene_name)

mcols(gtf)$gene_id <- fixed_gene_id
mcols(gtf)$gene_name <- fixed_gene_name

# --- Export
rtracklayer::export(gtf, opt$output, format = "gtf")

cat(sprintf("[INFO] Cleaned GTF written to: %s\n", opt$output))
