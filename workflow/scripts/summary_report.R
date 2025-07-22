#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

# --- CLI ---
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Filtered DESeq2 results CSV"),
  make_option(c("-o", "--output"), type = "character", help = "Output HTML file"),
  make_option(c("-l", "--label"), type = "character", default = "feature", help = "Feature label (e.g. miRNA, gene, transcript) for display")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Paths ---
bioproject_dir <- dirname(opt$input)
get_path <- function(name) file.path(bioproject_dir, paste0("deseq_results_", name))

paths <- list(
  volcano = get_path("volcano.pdf"),
  pca = get_path("PCA.pdf"),
  heatmap = get_path("heatmap.pdf"),
  go = get_path("GO_dotplot.pdf"),
  kegg = get_path("KEGG_dotplot.pdf"),
  reactome = get_path("Reactome_dotplot.pdf")
)

# --- Output ---
sink(opt$output)
cat("<html><head><title>DESeq2 Summary Report</title></head><body style='font-family:sans-serif;'>\n")
cat(sprintf("<h1>DESeq2 Summary Report: %s</h1>\n", basename(bioproject_dir)))

# --- Top DE features ---
cat(sprintf("<h2>Top Differentially Expressed %ss</h2>\n", opt$label))
df <- read_csv(opt$input, show_col_types = FALSE)

if (nrow(df) == 0 && opt$label == "miRNA") {
  fallback_path <- file.path(bioproject_dir, "all_mirnas.csv")
  if (file.exists(fallback_path)) {
    cat("<p><i>No filtered DE miRNAs found — falling back to all_mirnas.csv</i></p>\n")
    df <- read_csv(fallback_path, show_col_types = FALSE)
  } else {
    cat("<p><i>No filtered DE miRNAs found — no table shown (no fallback available).</i></p>\n")
    df <- data.frame()
  }
} else if (nrow(df) == 0) {
  cat(sprintf("<p><i>No filtered DE %ss found — no table shown.</i></p>\n", opt$label))
  df <- data.frame()
}


if ("...1" %in% names(df)) names(df)[names(df) == "...1"] <- "id" else names(df)[1] <- "id"  # fallback if header exists but not named

df <- df %>% arrange(padj) %>% select(id, log2FoldChange, padj) %>% head(20)


cat(sprintf("<table border=1 cellpadding=5 cellspacing=0>\n<tr><th>%s</th><th>log2FC</th><th>padj</th></tr>\n", opt$label))
for (i in seq_len(nrow(df))) {
  row <- df[i, ]
  cat(sprintf("<tr><td>%s</td><td>%.3f</td><td>%.2e</td></tr>\n", row$id, row$log2FoldChange, row$padj))
}

cat("</table>\n")

# --- Embed Plots ---
embed_plot <- function(label, path) {
  if (file.exists(path)) {
    cat(sprintf("<h2>%s</h2>\n", label))
    png_path <- sub("\\.pdf$", ".png", path)
    system(sprintf("convert -density 150 '%s' -quality 90 '%s'", path, png_path), ignore.stderr = TRUE)
    if (file.exists(png_path)) {
      cat(sprintf("<img src='%s' width='800px'><br>\n", basename(png_path)))
    } else {
      cat(sprintf("<p><i>Could not render %s</i></p>\n", label))
    }
  } else {
    cat(sprintf("<p><i>%s not found.</i></p>\n", label))
  }
}

embed_plot("Principal Component Analysis (PCA)", paths$pca)
embed_plot("Volcano Plot", paths$volcano)
embed_plot("Heatmap", paths$heatmap)
embed_plot("GO Enrichment", paths$go)
embed_plot("KEGG Enrichment", paths$kegg)
embed_plot("Reactome Enrichment", paths$reactome)

cat("</body></html>\n")
sink()
