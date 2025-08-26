#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
})

# ---------------- CLI ----------------
option_list <- list(
  make_option(c("-i","--input"),  type="character", help="CSV from extract_mirna_names.R (miRNA, precursor, log2FoldChange, padj)"),
  make_option(c("-t","--targets"),type="character", help="RDS with validated targets (from cache_mirna_targets.R)"),
  make_option(c("-o","--outdir"), type="character", help="Output directory"),
  make_option("--padj-threshold",  type="double", default=0.02, help="miRNA padj cutoff for significance [default: %default]"),
  make_option("--logfc-threshold", type="double", default=0.2,  help="|miRNA log2FC| cutoff for significance [default: %default]"),
  make_option("--min-gs-size",     type="integer", default=5,   help="Minimum gene set size for ORA [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), !is.null(opt$targets), !is.null(opt$outdir))
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

# ------------- helpers -------------
`%||%` <- function(a,b) if (is.null(a)) b else a

write_empty_ora_csv <- function(path) {
  writeLines('"ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count"', path)
}

dotplot_or_blank <- function(res, title, subtitle = NULL) {
  # wrap long subtitles to avoid running off the page
  if (!is.null(subtitle)) {
    subtitle <- paste(strwrap(subtitle, width = 80), collapse = "\n")
  }

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    p <- ggplot() + theme_void() + ggtitle(paste0(title, " (no terms)"))
    if (!is.null(subtitle)) p <- p + labs(subtitle = subtitle)
    return(p)
  }

  dotplot(res, showCategory = 20) +
    scale_y_discrete(labels = function(x) stringr::str_trunc(stringr::str_wrap(x, 40), 80)) +
    theme(
      axis.text.y   = element_text(size = 9),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.margin   = margin(5.5, 20, 5.5, 5.5)  # a little right padding for long labels
    ) +
    ggtitle(title) +
    (if (!is.null(subtitle)) labs(subtitle = subtitle) else NULL)
}

nnz <- function(x) sum(!is.na(x) & x != "")
clean_ids <- function(x) unique(as.character(na.omit(x)))

ok_set <- function(set, base, minGS=5) nnz(set) >= minGS && nnz(base) >= 50

run_GO <- function(genes, universe, ont, minGS)
  suppressMessages(tryCatch(
    enrichGO(gene=genes, universe=universe, OrgDb=org.Hs.eg.db,
             keyType="ENTREZID", ont=ont, pAdjustMethod="BH",
             pvalueCutoff=0.05, qvalueCutoff=0.2, minGSSize=minGS, readable=FALSE),
    error=function(e) NULL))

run_KEGG <- function(genes, universe, minGS)
  suppressMessages(tryCatch(
    enrichKEGG(gene=genes, universe=universe, organism="hsa",
               pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, minGSSize=minGS),
    error=function(e) NULL))

run_Reactome <- function(genes, universe, minGS)
  suppressMessages(tryCatch(
    enrichPathway(gene=genes, universe=universe, organism="human",
                  pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                  minGSSize=minGS, readable=TRUE),
    error=function(e) NULL))

# ------------- I/O -------------
res <- suppressMessages(read_csv(opt$input, show_col_types = FALSE))
if (!all(c("miRNA","log2FoldChange","padj") %in% names(res)))
  stop("Input CSV missing required columns")

targets <- readRDS(opt$targets)
if (!all(c("mature_mirna_id","target_symbol","target_entrez","database") %in% names(targets)))
  stop("Targets RDS missing required columns")

# ------------- filter miRNAs -------------
sig_tbl <- res %>% filter(!is.na(padj) & padj <= opt$`padj-threshold`,
                          !is.na(log2FoldChange) & abs(log2FoldChange) >= opt$`logfc-threshold`)

# if none, still emit empty files
if (nrow(sig_tbl) == 0) {
  message("[INFO] No significant miRNAs found — writing empty outputs.")
  # mapping, entrez, counts
  write.csv(data.frame(), file.path(opt$outdir, "deseq_results_mirna_targets.csv"), row.names=FALSE)
  write.csv(data.frame(), file.path(opt$outdir, "deseq_results_target_genes.csv"), row.names=FALSE)
  writeLines("Enrichment summary: 0 miRNAs, 0 targets", file.path(opt$outdir, "enrichment_counts.txt"))
  # ORA empty
  for (ont in c("BP","CC","MF")) for (subset in c("ALL","UP","DOWN"))
    write_empty_ora_csv(file.path(opt$outdir, sprintf("ORA_GO_%s_%s_enrichment.csv", ont, subset)))
  for (subset in c("ALL","UP","DOWN")) {
    write_empty_ora_csv(file.path(opt$outdir, sprintf("ORA_KEGG_%s_enrichment.csv", subset)))
    write_empty_ora_csv(file.path(opt$outdir, sprintf("ORA_Reactome_%s_enrichment.csv", subset)))
  }
  pdf(file.path(opt$outdir, "ORA_GO_dotplot_ALL_UP_DOWN_BPCCMF.pdf")); plot.new(); title("No GO enrichment"); dev.off()
  pdf(file.path(opt$outdir, "ORA_KEGG_dotplot.pdf")); plot.new(); title("No KEGG enrichment"); dev.off()
  pdf(file.path(opt$outdir, "ORA_Reactome_dotplot.pdf")); plot.new(); title("No Reactome enrichment"); dev.off()
  quit(save="no", status=0)
}

# ------------- build target sets -------------
sig_names <- unique(sig_tbl$miRNA)
tar_sig <- targets %>% filter(mature_mirna_id %in% sig_names) %>%
  mutate(target_entrez = as.character(target_entrez))

# export mapping + entrez list
write.csv(tar_sig[,c("mature_mirna_id","target_symbol","target_entrez","database")],
          file.path(opt$outdir,"deseq_results_mirna_targets.csv"), row.names=FALSE)
write.csv(data.frame(EntrezID=unique(tar_sig$target_entrez)),
          file.path(opt$outdir,"deseq_results_target_genes.csv"), row.names=FALSE)

# summary counts
sink(file.path(opt$outdir,"enrichment_counts.txt"))
cat("Significant miRNAs:", nrow(sig_tbl), "\n")
cat("Validated targets (rows):", nrow(tar_sig), "\n")
cat("Unique Entrez IDs:", length(unique(tar_sig$target_entrez)), "\n")
sink()

# define ALL/UP/DOWN
mi_up   <- sig_tbl$miRNA[sig_tbl$log2FoldChange > 0]
mi_down <- sig_tbl$miRNA[sig_tbl$log2FoldChange < 0]

ent_ALL  <- clean_ids(tar_sig$target_entrez)
ent_UP   <- clean_ids(tar_sig$target_entrez[tar_sig$mature_mirna_id %in% mi_up])
ent_DOWN <- clean_ids(tar_sig$target_entrez[tar_sig$mature_mirna_id %in% mi_down])

ent_UNIV <- clean_ids(keys(org.Hs.eg.db, keytype="ENTREZID") %||% character(0))

minGS <- opt$`min-gs-size`
run_or_save <- function(fun, genes, univ, path, label) {
  if (ok_set(genes, univ, minGS)) {
    obj <- fun(genes, univ, minGS)
    if (!is.null(obj) && nrow(as.data.frame(obj)) > 0) {
      write.csv(as.data.frame(obj), path, row.names=FALSE); return(obj)
    }
  }
  write_empty_ora_csv(path); return(NULL)
}

# --- GO ORA ---
go_list <- list()
for (ont in c("BP","CC","MF")) {
  go_list[[paste0(ont,"_ALL")]]  <- run_or_save(function(g,u,ms) run_GO(g,u,ont,ms), ent_ALL, ent_UNIV,
                                                file.path(opt$outdir,sprintf("ORA_GO_%s_ALL_enrichment.csv",ont)), ont)
  go_list[[paste0(ont,"_UP")]]   <- run_or_save(function(g,u,ms) run_GO(g,u,ont,ms), ent_UP,  ent_UNIV,
                                                file.path(opt$outdir,sprintf("ORA_GO_%s_UP_enrichment.csv",ont)), ont)
  go_list[[paste0(ont,"_DOWN")]] <- run_or_save(function(g,u,ms) run_GO(g,u,ont,ms), ent_DOWN,ent_UNIV,
                                                file.path(opt$outdir,sprintf("ORA_GO_%s_DOWN_enrichment.csv",ont)), ont)
}

# --- KEGG ORA ---
kegg_ALL  <- run_or_save(run_KEGG, ent_ALL,  ent_UNIV, file.path(opt$outdir,"ORA_KEGG_ALL_enrichment.csv"),"KEGG")
kegg_UP   <- run_or_save(run_KEGG, ent_UP,   ent_UNIV, file.path(opt$outdir,"ORA_KEGG_UP_enrichment.csv"),"KEGG")
kegg_DOWN <- run_or_save(run_KEGG, ent_DOWN, ent_UNIV, file.path(opt$outdir,"ORA_KEGG_DOWN_enrichment.csv"),"KEGG")

# --- Reactome ORA ---
react_ALL  <- run_or_save(run_Reactome, ent_ALL,  ent_UNIV, file.path(opt$outdir,"ORA_Reactome_ALL_enrichment.csv"),"Reactome")
react_UP   <- run_or_save(run_Reactome, ent_UP,   ent_UNIV, file.path(opt$outdir,"ORA_Reactome_UP_enrichment.csv"),"Reactome")
react_DOWN <- run_or_save(run_Reactome, ent_DOWN, ent_UNIV, file.path(opt$outdir,"ORA_Reactome_DOWN_enrichment.csv"),"Reactome")

# --- PDFs ---
pdf(file.path(opt$outdir,"ORA_GO_dotplot_ALL_UP_DOWN_BPCCMF.pdf"))
for (ont in c("BP","CC","MF")) for (subset in c("ALL","UP","DOWN")) {
  obj <- go_list[[paste0(ont,"_",subset)]]
  print(dotplot_or_blank(obj, paste0("GO ORA [",ont,"] (",subset,")"),
                         sprintf("STRICT: p.adjust≤%.2g; minGS=%d | set=%d, univ≈%d",
                                 0.05, minGS,
                                 if (subset=="ALL") nnz(ent_ALL) else if (subset=="UP") nnz(ent_UP) else nnz(ent_DOWN),
                                 nnz(ent_UNIV))))
}
dev.off()

pdf(file.path(opt$outdir,"ORA_KEGG_dotplot.pdf"))
for (subset in c("ALL","UP","DOWN")) {
  obj <- switch(subset,ALL=kegg_ALL,UP=kegg_UP,DOWN=kegg_DOWN)
  size <- switch(subset,ALL=nnz(ent_ALL),UP=nnz(ent_UP),DOWN=nnz(ent_DOWN))
  print(dotplot_or_blank(obj,paste0("KEGG ORA (",subset,")"),
                         sprintf("STRICT: p.adjust≤%.2g; minGS=%d | set=%d, univ≈%d",
                                 0.05,minGS,size,nnz(ent_UNIV))))
}
dev.off()

pdf(file.path(opt$outdir,"ORA_Reactome_dotplot.pdf"))
for (subset in c("ALL","UP","DOWN")) {
  obj <- switch(subset,ALL=react_ALL,UP=react_UP,DOWN=react_DOWN)
  size <- switch(subset,ALL=nnz(ent_ALL),UP=nnz(ent_UP),DOWN=nnz(ent_DOWN))
  print(dotplot_or_blank(obj,paste0("Reactome ORA (",subset,")"),
                         sprintf("STRICT: p.adjust≤%.2g; minGS=%d | set=%d, univ≈%d",
                                 0.05,minGS,size,nnz(ent_UNIV))))
}
dev.off()

message("[INFO] Enrichment complete.")
