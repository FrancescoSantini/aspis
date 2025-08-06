#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# --- CLI options ---
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CSV with filtered DESeq2 results (gene or transcript IDs as rownames)"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory"),
  make_option(c("--padj-threshold"), type = "double", default = 0.1, help = "Adjusted p-value threshold [default: %default]"),
  make_option(c("--logfc-threshold"), type = "double", default = 0, help = "Minimum absolute log2FC threshold [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- Read DE results ---
de <- read_csv(opt$input, show_col_types = FALSE)
id_col <- if ("...1" %in% names(de)) "...1" else names(de)[1]
de[[id_col]] <- as.character(de[[id_col]])

# Filter by padj and log2FC
sig <- de %>%
  filter(!is.na(padj), padj < opt$`padj-threshold`, abs(log2FoldChange) > opt$`logfc-threshold`)

if (nrow(sig) == 0) {
  cat("No significant entries after filtering. Exiting.\n")
  quit(save = "no", status = 0)
}

# Extract Ensembl or gene symbols
gene_ids <- sig[[id_col]]
gene_ids <- gsub("\\|.*", "", gene_ids)

if (all(grepl("^ENST", gene_ids))) {
  from_key <- "ENSEMBLTRANS"
} else if (all(grepl("^ENSG", gene_ids))) {
  from_key <- "ENSEMBL"
} else if (all(grepl("^MSTRG", gene_ids))) {
  cat("[INFO] Skipping enrichment: only MSTRG IDs detected.\n")
  quit(save = "no", status = 0)
} else {
  from_key <- "ENSEMBL"  # fallback
}

entrez <- bitr(gene_ids, fromType = from_key, toType = "ENTREZID", OrgDb = org.Hs.eg.db)

entrez_ids <- na.omit(entrez$ENTREZID)

if (length(entrez_ids) == 0) {
  cat("No Entrez IDs mapped. Exiting.\n")
  quit(save = "no", status = 0)
}

# Output files
go_file <- file.path(opt$outdir, "GO_enrichment.csv")
go_plot <- file.path(opt$outdir, "GO_dotplot.pdf")
kegg_file <- file.path(opt$outdir, "KEGG_enrichment.csv")
kegg_plot <- file.path(opt$outdir, "KEGG_dotplot.pdf")
reactome_file <- file.path(opt$outdir, "Reactome_enrichment.csv")
reactome_plot <- file.path(opt$outdir, "Reactome_dotplot.pdf")
summary_file <- file.path(opt$outdir, "all_enrichments_summary.csv")

# --- GO ---
go <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = opt$`padj-threshold`,
  readable = TRUE
)

if (!is.null(go) && nrow(go) > 0) {
  write.csv(as.data.frame(go), file = go_file)
  pdf(go_plot); print(dotplot(go, showCategory = 20) + ggtitle("GO Biological Process")); dev.off()
} else {
  file.create(go_file)
  pdf(go_plot); plot.new(); title("No GO terms enriched"); dev.off()
}

# --- KEGG ---
kegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = opt$`padj-threshold`
)

if (!is.null(kegg) && nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg), file = kegg_file)
  pdf(kegg_plot); print(dotplot(kegg, showCategory = 20) + ggtitle("KEGG Pathways")); dev.off()
} else {
  file.create(kegg_file)
  pdf(kegg_plot); plot.new(); title("No KEGG pathways enriched"); dev.off()
}

# --- Reactome ---
reactome <- enrichPathway(
  gene = entrez_ids,
  organism = "human",
  pAdjustMethod = "BH",
  qvalueCutoff = opt$`padj-threshold`
)

if (!is.null(reactome) && nrow(reactome) > 0) {
  write.csv(as.data.frame(reactome), file = reactome_file)
  pdf(reactome_plot); print(dotplot(reactome, showCategory = 20) + ggtitle("Reactome Pathways")); dev.off()
} else {
  file.create(reactome_file)
  pdf(reactome_plot); plot.new(); title("No Reactome terms enriched"); dev.off()
}

# --- Combine summaries
all <- list()
if (!is.null(go) && nrow(go) > 0) all$GO <- as.data.frame(go)
if (!is.null(kegg) && nrow(kegg) > 0) all$KEGG <- as.data.frame(kegg)
if (!is.null(reactome) && nrow(reactome) > 0) all$Reactome <- as.data.frame(reactome)

if (length(all) > 0) {
  combined <- bind_rows(lapply(names(all), function(name) {
    df <- all[[name]]
    df$Source <- name
    df
  }))
  write.csv(combined, summary_file, row.names = FALSE)
} else {
  file.create(summary_file)
}

cat("[INFO] Enrichment analysis complete.\n")
