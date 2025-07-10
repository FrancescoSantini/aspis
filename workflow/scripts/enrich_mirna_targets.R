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
  make_option(c("-i", "--input"), type = "character", help = "CSV with mature miRNA list and stats"),
  make_option(c("-t", "--targets"), type = "character", help = "RDS from multiMiR"),
  make_option(c("-o", "--outdir"), type = "character", help = "Output directory"),
  make_option(c("--padj-threshold"), type = "double", default = 0.1, help = "Adjusted p-value threshold [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Read miRNA list ---
res <- read_csv(opt$input, show_col_types = FALSE)
if (!"miRNA" %in% names(res)) stop("Missing 'miRNA' column in input.")
mirnas <- unique(na.omit(res$miRNA))
if (length(mirnas) == 0) {
  cat("[INFO] No miRNAs to process — skipping enrichment.\n")
  quit(save = "no", status = 0)
}

# --- Load validated targets ---
cat("[INFO] Loading cached validated targets from RDS...\n")
if (!file.exists(opt$targets)) stop("Cached RDS file not found: ", opt$targets)
all_targets <- readRDS(opt$targets)

required_cols <- c("mature_mirna_id", "target_symbol", "target_entrez", "database")
if (!all(required_cols %in% names(all_targets))) {
  stop("RDS missing required columns: ", paste(required_cols, collapse = ", "))
}

targets <- all_targets[all_targets$mature_mirna_id %in% mirnas, ]
if (nrow(targets) == 0) {
  cat("[INFO] No validated targets matched — skipping enrichment.\n")
  quit(save = "no", status = 0)
}

# --- Output paths ---
go_file       <- file.path(opt$outdir, "deseq_results_GO_enrichment.csv")
go_plot       <- file.path(opt$outdir, "deseq_results_GO_dotplot.pdf")
kegg_file     <- file.path(opt$outdir, "deseq_results_KEGG_enrichment.csv")
kegg_plot     <- file.path(opt$outdir, "deseq_results_KEGG_dotplot.pdf")
reactome_file <- file.path(opt$outdir, "deseq_results_Reactome_enrichment.csv")
reactome_plot <- file.path(opt$outdir, "deseq_results_Reactome_dotplot.pdf")
summary_file  <- file.path(opt$outdir, "deseq_results_all_enrichments_summary.csv")
log_file      <- file.path(opt$outdir, "enrichment_counts.txt")
entrez_file   <- file.path(opt$outdir, "deseq_results_target_genes.csv")
mapping_file  <- file.path(opt$outdir, "deseq_results_mirna_targets.csv")

# --- Save targets ---
entrez_ids <- na.omit(unique(targets$target_entrez))
target_pairs <- unique(targets[, c("mature_mirna_id", "target_symbol", "target_entrez", "database")])
write_csv(data.frame(EntrezID = entrez_ids), entrez_file)
write_csv(target_pairs, mapping_file)
cat(sprintf("[INFO] %d validated targets mapped to Entrez IDs.\n", length(entrez_ids)))

# --- GO enrichment ---
if (length(entrez_ids) > 0) {
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
    pdf(go_plot); print(dotplot(go, showCategory = 20) + ggtitle("GO Biological Process Enrichment")); dev.off()
  } else {
    file.create(go_file); pdf(go_plot); plot.new(); title("No GO enrichment results"); dev.off()
  }
}

# --- KEGG enrichment ---
valid_kegg_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
entrez_ids_kegg <- intersect(entrez_ids, valid_kegg_ids)

if (length(entrez_ids_kegg) > 0) {
  kegg <- enrichKEGG(
    gene = entrez_ids_kegg,
    organism = 'hsa',
    pAdjustMethod = "BH",
    qvalueCutoff = opt$`padj-threshold`
  )
  if (!is.null(kegg) && nrow(kegg) > 0) {
    write.csv(as.data.frame(kegg), file = kegg_file)
    pdf(kegg_plot); print(dotplot(kegg, showCategory = 20) + ggtitle("KEGG Pathway Enrichment")); dev.off()
  } else {
    file.create(kegg_file); pdf(kegg_plot); plot.new(); title("No KEGG enrichment results"); dev.off()
  }
}

# --- Reactome enrichment ---
if (length(entrez_ids) > 0) {
  reactome <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pAdjustMethod = "BH",
    qvalueCutoff = opt$`padj-threshold`
  )
  if (!is.null(reactome) && nrow(reactome) > 0) {
    write.csv(as.data.frame(reactome), file = reactome_file)
    pdf(reactome_plot); print(dotplot(reactome, showCategory = 20) + ggtitle("Reactome Pathway Enrichment")); dev.off()
  } else {
    file.create(reactome_file); pdf(reactome_plot); plot.new(); title("No Reactome enrichment results"); dev.off()
  }
}

# --- Combined summary ---
all_enrichments <- list()
add_enrichment <- function(obj, label) {
  if (!is.null(obj) && nrow(obj) > 0) {
    df <- as.data.frame(obj)
    df$Source <- label
    all_enrichments[[label]] <<- df
  }
}
if (exists("go")) add_enrichment(go, "GO")
if (exists("kegg")) add_enrichment(kegg, "KEGG")
if (exists("reactome")) add_enrichment(reactome, "Reactome")

if (length(all_enrichments) > 0) {
  all_cols <- unique(unlist(lapply(all_enrichments, names)))
  all_enrichments <- lapply(all_enrichments, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[all_cols]
  })
  combined <- do.call(rbind, all_enrichments)
  write.csv(combined, file = summary_file, row.names = FALSE)
} else {
  file.create(summary_file)
}

# --- Summary log ---
sink(log_file)
cat("Enrichment summary (validated targets only):\n")
cat(sprintf("miRNAs:             %d\n", length(mirnas)))
cat(sprintf("Validated targets:  %d\n", nrow(target_pairs)))
cat(sprintf("Unique Entrez IDs:  %d\n", length(entrez_ids)))
cat(sprintf("GO terms:           %d\n", if (exists("go") && !is.null(go)) nrow(go) else 0))
cat(sprintf("KEGG pathways:      %d\n", if (exists("kegg") && !is.null(kegg)) nrow(kegg) else 0))
cat(sprintf("Reactome terms:     %d\n", if (exists("reactome") && !is.null(reactome)) nrow(reactome) else 0))
sink()
