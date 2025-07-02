#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(miRBaseConverter)
  library(dplyr)
  library(readr)
  library(tibble)
})

# --- Parse CLI ---
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CSV input file"),
  make_option(c("-o", "--output"), type = "character", help = "Output CSV file"),
  make_option(c("-n", "--topn"), type = "integer", default = 30, help = "Top N by padj and log2FC")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Read input ---
if (!file.exists(opt$input)) stop("Input file not found: ", opt$input)
df <- read_csv(opt$input, show_col_types = FALSE)

cat("[DEBUG] Headers of loaded CSV:\n")
print(colnames(df))

cat("[DEBUG] First rows of data:\n")
print(head(df))

# Rename unnamed rowname column if present
if ("...1" %in% colnames(df)) {
  cat("[DEBUG] Renaming ...1 to precursor\n")
  colnames(df)[which(colnames(df) == "...1")] <- "precursor"
}

cat("[DEBUG] Column names after normalization:\n")
print(colnames(df))

# --- Filter and select top N ---
df <- df %>% filter(!is.na(padj))
cat(sprintf("[DEBUG] Total miRNAs with padj: %d\n", nrow(df)))

# Top by padj and FC
top_padj <- df %>% arrange(padj) %>% head(opt$topn)
top_fc   <- df %>% arrange(desc(abs(log2FoldChange))) %>% head(opt$topn)
selected <- bind_rows(top_padj, top_fc) %>% distinct(precursor, .keep_all = TRUE)

cat(sprintf("[DEBUG] Selected top %d unique miRNAs\n", nrow(selected)))
cat("[DEBUG] FULL selected dataframe:\n")
print(selected)

# --- Map precursors to mature ---
cat("[DEBUG] Mapping precursor names to mature...\n")
map <- miRNA_PrecursorToMature(selected$precursor)

cat("[DEBUG] FULL precursor-to-mature mapping:\n")
print(map)

if (nrow(map) == 0) {
  warning("[WARN] No mature miRNAs mapped.")
  write_csv(tibble(miRNA = character(), precursor = character(), log2FoldChange = numeric(), padj = numeric()), opt$output)
  quit(save = "no", status = 0)
}

mature_list <- unique(na.omit(c(map$Mature1, map$Mature2)))
cat("[DEBUG] Mature miRNA candidates:\n")
print(mature_list)

# --- Build output
mature_df <- bind_rows(lapply(mature_list, function(m) {
  cat(sprintf("[DEBUG] Processing mature: %s\n", m))
  prec_matches <- map$OriginalName[which(map$Mature1 == m | map$Mature2 == m)]
  cat("[DEBUG] Matching precursors:\n")
  print(prec_matches)

  prec <- prec_matches[prec_matches %in% selected$precursor][1]
  cat(sprintf("[DEBUG] Selected precursor match: %s\n", prec))

  if (is.na(prec) || is.null(prec)) {
    cat(sprintf("[WARN] Skipping mature miRNA %s — no precursor match in selected.\n", m))
    return(NULL)
  }

  row <- selected %>% filter(precursor == prec)
  if (nrow(row) == 0) {
    cat(sprintf("[WARN] Precursor %s not found in selected rows. Skipping %s\n", prec, m))
    return(NULL)
  }

  row <- row %>% slice(1)
  tibble(
    miRNA = m,
    precursor = prec,
    log2FoldChange = row$log2FoldChange,
    padj = row$padj
  )
}))

if (nrow(mature_df) == 0) {
  warning("[WARN] No data to write. Output will be empty.")
} else {
  cat(sprintf("[DEBUG] Final mature miRNAs: %d\n", nrow(mature_df)))
}

write_csv(mature_df, opt$output)
cat(sprintf("[INFO] Output written: %s\n", opt$output))
