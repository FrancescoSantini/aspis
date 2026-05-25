#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg) == 0) {
  stop("Cannot locate wrapper script path from commandArgs()")
}
wrapper_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
feature_script <- file.path(dirname(wrapper_path), "run_deseq2_feature.R")
source(feature_script, chdir = TRUE)
