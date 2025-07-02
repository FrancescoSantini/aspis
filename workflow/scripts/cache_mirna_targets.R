#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(multiMiR)
  library(dplyr)
  library(readr)
})

# --- CLI Arguments ---
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CSV file with miRNAs, log2FoldChange, padj"),
  make_option(c("-o", "--output"), type = "character", help = "Output .rds file"),
  make_option(c("-l", "--log"), type = "character", default = NULL, help = "Optional log file")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$log)) sink(opt$log, split = TRUE)

# --- Load input ---
if (!file.exists(opt$input)) stop("Input file not found: ", opt$input)
df <- read_csv(opt$input, show_col_types = FALSE)

required <- c("miRNA", "log2FoldChange", "padj")
if (!all(required %in% names(df))) {
  stop("Input must contain columns: ", paste(required, collapse = ", "))
}

cat(sprintf("[INFO] Loaded %d miRNAs\n", nrow(df)))
print(head(df, 5))

results <- list()
failed <- character(0)

# --- Query each miRNA one by one ---
for (i in seq_len(nrow(df))) {
  row <- df[i, ]
  mirna <- row$miRNA
  lfc <- row$log2FoldChange
  padj <- row$padj

  cat(sprintf("[INFO] [%3d/%3d] Querying %s | log2FC=%.2f | padj=%.3g\n", i, nrow(df), mirna, lfc, padj))

  res <- tryCatch(
    get_multimir(mirna = mirna, table = "validated", summary = TRUE, legacy.out = FALSE),
    error = function(e) {
      warning(sprintf("[WARN] Query failed for %s: %s", mirna, e$message))
      return(NULL)
    }
  )

  if (!is.null(res) && "data" %in% slotNames(res)) {
    dat <- res@data
    if (!is.null(dat) && nrow(dat) > 0) {
      dat <- mutate(dat,
                    input_miRNA = mirna,
                    log2FoldChange = lfc,
                    padj = padj)
      results[[mirna]] <- dat
    } else {
      cat(sprintf("[DEBUG] No valid data for %s\n", mirna))
      failed <- c(failed, mirna)
    }
  } else {
    cat(sprintf("[DEBUG] Invalid object returned for %s\n", mirna))
    failed <- c(failed, mirna)
  }
}

# --- Save output ---
if (length(results) == 0) {
  warning("[WARN] No valid results. Saving empty RDS.")
  saveRDS(data.frame(), opt$output)
} else {
  combined <- bind_rows(results)
  cat(sprintf("[INFO] Retrieved %d interactions across %d miRNAs.\n", nrow(combined), length(results)))
  saveRDS(combined, opt$output)
  cat(sprintf("[INFO] Saved to: %s\n", opt$output))
}

# --- Save failures ---
if (length(failed) > 0) {
  failfile <- sub("\\.rds$", "_failed_mirnas.txt", opt$output)
  writeLines(failed, failfile)
  cat(sprintf("[INFO] %d miRNAs failed. Logged to %s\n", length(failed), failfile))
}
