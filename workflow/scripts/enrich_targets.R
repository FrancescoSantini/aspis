#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
})

# ---------------- CLI ----------------
option_list <- list(
  make_option(c("-i","--input"),  type="character", help="Filtered DE CSV"),
  make_option(c("-r","--raw"),    type="character", help="Raw DE CSV"),
  make_option(c("-o","--outdir"), type="character", help="Output directory"),
  make_option("--gene-meta",      type="character", default=NA, help="Path to gene_metadata.csv"),
  make_option("--tx-meta",        type="character", default=NA, help="Path to transcript_metadata.csv"),
  make_option("--padj-threshold",  type="double",   default=0.05),
  make_option("--logfc-threshold", type="double",   default=0.5),
  make_option("--do-gsea", action="store_true", default=FALSE),
  make_option("--kegg_org",  type="character", default="hsa"),
  make_option("--organism",  type="character", default="human"),
  make_option("--min-gs-size",   type="integer", default=5),
  make_option("--gsea-q-cutoff", type="double",  default=0.25),
  make_option("--pc-list", type="character", default=NA,
              help="Optional path to pc_gene_ids.txt (one ENSG per line). Applied AFTER mapping to ENSG.")
)
opt <- parse_args(OptionParser(option_list=option_list))
stopifnot(!is.null(opt$input), !is.null(opt$raw), !is.null(opt$outdir))
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)

# ---------------- Helpers ----------------
`%||%` <- function(a,b) if (is.null(a)) b else a
strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))
looks_ensg <- function(x) grepl("^ENSG", x %||% "", ignore.case = FALSE)
looks_enst <- function(x) grepl("^ENST", x %||% "", ignore.case = FALSE)

wrap_sub <- function(x, width = 80) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  paste(strwrap(x, width = width), collapse = "\n")
}

discover_gene_meta <- function(outdir, cli_path){
  if (!is.na(cli_path) && nzchar(cli_path) && file.exists(cli_path)) return(cli_path)
  f <- file.path(dirname(outdir), "gene_metadata.csv")
  if (file.exists(f)) f else NA_character_
}
discover_tx_meta <- function(outdir, cli_path){
  if (!is.na(cli_path) && nzchar(cli_path) && file.exists(cli_path)) return(cli_path)
  f <- file.path(dirname(outdir), "transcript_metadata.csv")
  if (file.exists(f)) f else NA_character_
}

GENE_META <- discover_gene_meta(opt$outdir, opt$`gene-meta`)
TX_META   <- discover_tx_meta(opt$outdir,   opt$`tx-meta`)
writeLines(sprintf("gene_meta=%s\ntx_meta=%s",
                   ifelse(is.na(GENE_META),"NA",GENE_META),
                   ifelse(is.na(TX_META),"NA",TX_META)),
           file.path(opt$outdir, "debug_meta.txt"))

# ------- Metadata helpers -------
sym_from_gene_meta <- function(ensg_vec){
  if (is.na(GENE_META)) return(rep(NA_character_, length(ensg_vec)))
  gm <- suppressWarnings(tryCatch(fread(GENE_META, showProgress=FALSE), error=function(e) NULL))
  if (is.null(gm)) return(rep(NA_character_, length(ensg_vec)))
  need <- c("feature_id","gene_name")
  if (!all(need %in% names(gm))) return(rep(NA_character_, length(ensg_vec)))
  key <- data.table(feature_id = strip_version(as.character(gm$feature_id)),
                    gene_name  = as.character(gm$gene_name))
  setkey(key, feature_id)
  m <- key[.(strip_version(ensg_vec))]
  m$gene_name
}

tx2gene_from_meta <- function(txs){
  if (is.na(TX_META)) return(rep(NA_character_, length(txs)))
  dt <- suppressWarnings(tryCatch(fread(TX_META, showProgress=FALSE), error=function(e) NULL))
  if (is.null(dt)) return(rep(NA_character_, length(txs)))
  tcol <- intersect(c("transcript_id","feature_id","tx_id","tx"), names(dt))
  gcol <- intersect(c("gene_id","parent_gene","gene"), names(dt))
  if (!length(tcol) || !length(gcol)) return(rep(NA_character_, length(txs)))
  dt[, tx   := strip_version(get(tcol[1]))]
  dt[, gene := strip_version(get(gcol[1]))]
  setkey(dt, tx)
  dt[.(strip_version(txs)), gene]
}

# ====== (1) Robust ID -> ENSG mapping ======
to_gene_ensembl <- function(ids){
  ids <- strip_version(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  if (!length(ids)) return(character(0))

  keep_ensg <- ids[looks_ensg(ids)]
  rest      <- ids[!looks_ensg(ids)]

  ensg_from_enst <- character(0)
  if (length(rest)) {
    txs <- unique(rest[looks_enst(rest)])
    if (length(txs)) {
      m1 <- suppressWarnings(tryCatch(
        AnnotationDbi::select(org.Hs.eg.db, keys=txs, keytype="ENSEMBLTRANS", columns="ENSEMBL"),
        error=function(e) NULL))
      if (!is.null(m1) && nrow(m1)) ensg_from_enst <- unique(strip_version(m1$ENSEMBL))
      need <- setdiff(strip_version(txs), unique(strip_version(m1$ENSEMBLTRANS %||% character(0))))
      if (length(need)) ensg_from_enst <- unique(c(ensg_from_enst, strip_version(tx2gene_from_meta(need))))
    }
  }

  sym_mask <- !looks_ensg(rest) & !looks_enst(rest)
  ensg_from_sym <- character(0)
  if (any(sym_mask)) {
    uniq <- unique(rest[sym_mask])
    m <- suppressWarnings(tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys=uniq, keytype="SYMBOL", columns="ENSEMBL"),
      error=function(e) NULL))
    if (!is.null(m) && nrow(m)) ensg_from_sym <- unique(strip_version(m$ENSEMBL))
  }

  out <- unique(c(keep_ensg, ensg_from_enst, ensg_from_sym))
  out[looks_ensg(out)]
}

# ====== (2) ENSG -> ENTREZ with SYMBOL rescue ======
ensg_to_entrez <- function(ensg){
  ensg <- unique(strip_version(ensg))
  ensg <- ensg[looks_ensg(ensg)]
  if (!length(ensg)) return(data.frame(ENSEMBL=character(), ENTREZID=character()))

  m1 <- suppressWarnings(tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys=ensg, keytype="ENSEMBL", columns="ENTREZID"),
    error=function(e) NULL))
  m1 <- if (!is.null(m1) && nrow(m1)) unique(m1[,c("ENSEMBL","ENTREZID")]) else
        data.frame(ENSEMBL=character(),ENTREZID=character())

  need <- setdiff(ensg, unique(m1$ENSEMBL))
  if (length(need)) {
    sy <- sym_from_gene_meta(need)
    keep <- !is.na(sy) & sy!=""
    if (any(keep)) {
      m2 <- suppressWarnings(tryCatch(
        AnnotationDbi::select(org.Hs.eg.db, keys=unique(sy[keep]),
                              keytype="SYMBOL", columns=c("ENTREZID","SYMBOL")),
        error=function(e) NULL))
      if (!is.null(m2) && nrow(m2)) {
        dt <- data.table(ENSEMBL=need, SYMBOL=sy)
        m2 <- data.table(m2)
        j  <- m2[dt, on=.(SYMBOL), nomatch=0L][, .(ENSEMBL, ENTREZID)]
        if (nrow(j)) m1 <- unique(rbind(m1, as.data.frame(j)))
      }
    }
  }
  unique(m1[!is.na(m1$ENTREZID) & m1$ENTREZID!="", ])
}

# For GSEA vector mapping
map_ids_vector_to_ensg <- function(ids){
  ids <- strip_version(ids)
  out <- rep(NA_character_, length(ids))
  m <- looks_ensg(ids); out[m] <- ids[m]
  m <- looks_enst(ids)
  if (any(m)) {
    uniq <- unique(ids[m])
    m1 <- suppressWarnings(tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys=uniq, keytype="ENSEMBLTRANS", columns="ENSEMBL"),
      error=function(e) NULL))
    if (!is.null(m1) && nrow(m1)) {
      m1$ENSEMBLTRANS <- strip_version(m1$ENSEMBLTRANS)
      m1$ENSEMBL      <- strip_version(m1$ENSEMBL)
      setDT(m1); setkey(m1, ENSEMBLTRANS)
      out[m] <- m1[.(ids[m]), ENSEMBL]
    }
    need <- m & is.na(out)
    if (any(need)) out[need] <- strip_version(tx2gene_from_meta(ids[need]))
  }
  m <- is.na(out)
  if (any(m)) {
    uniq <- unique(ids[m])
    m2 <- suppressWarnings(tryCatch(
      AnnotationDbi::select(org.Hs.eg.db, keys=uniq, keytype="SYMBOL", columns="ENSEMBL"),
      error=function(e) NULL))
    if (!is.null(m2) && nrow(m2)) {
      m2$SYMBOL  <- as.character(m2$SYMBOL)
      m2$ENSEMBL <- strip_version(m2$ENSEMBL)
      setDT(m2); setkey(m2, SYMBOL)
      out[m] <- m2[.(ids[m]), ENSEMBL]
    }
  }
  out[!looks_ensg(out)] <- NA_character_
  out
}

# I/O helpers
save_enrich_csv <- function(res, path, kind=c("ORA","GSEA")){
  kind <- match.arg(kind)
  if (is.null(res) || nrow(as.data.frame(res))==0) {
    if (kind=="ORA")
      writeLines('"ID","Description","GeneRatio","BgRatio","RichFactor","FoldEnrichment","zScore","pvalue","p.adjust","qvalue","geneID","Count"', path)
    else
      writeLines('"ID","Description","p.adjust","qvalue"', path)
  } else {
    data.table::fwrite(as.data.frame(res), path)
  }
}

dotplot_wrapped <- function(res, title, subtitle = NULL){
  if (!is.null(subtitle)) subtitle <- wrap_sub(subtitle, 80)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    p <- ggplot() + theme_void() + ggtitle(paste0(title, " (no terms)"))
    if (!is.null(subtitle)) p <- p + labs(subtitle = subtitle)
    return(p)
  } else {
    dotplot(res, showCategory = 20) +
      scale_y_discrete(labels = function(x) str_trunc(str_wrap(x, 40), 80)) +
      theme(axis.text.y   = element_text(size = 9),
            plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(size = 9)) +
      ggtitle(title) +
      (if (!is.null(subtitle)) labs(subtitle = subtitle) else NULL)
  }
}

# ===== ORA plotting helpers =====
pdf_triplet <- function(res_list, pdffile, title_prefix="", subtitle=NULL){
  pdf(pdffile, width=8, height=6)
  for (nm in c("ALL","UP","DOWN")) {
    print(dotplot_wrapped(res_list[[nm]], paste0(title_prefix," (",nm,")"), subtitle))
  }
  dev.off()
}
pdf_go_9pages <- function(res_nested, pdffile, title_prefix="", subtitle=NULL){
  pdf(pdffile, width=8, height=6)
  for (ont in c("BP","CC","MF")) {
    for (nm in c("ALL","UP","DOWN")) {
      res <- res_nested[[ont]][[nm]]
      print(dotplot_wrapped(res, paste0(title_prefix," [", ont, "] (", nm, ")"), subtitle))
    }
  }
  dev.off()
}

# ---------------- Read input ----------------
filt <- fread(opt$input, showProgress=FALSE)
raw  <- fread(opt$raw,   showProgress=FALSE)
first_col_name <- names(filt)[1L]
writeLines(sprintf("First column name in filtered CSV: '%s'", first_col_name),
           file.path(opt$outdir, "debug_first_colname.txt"))

ids_filt <- strip_version(filt[[1L]])
ids_raw  <- strip_version(raw[[1L]])

if (!all(c("padj","log2FoldChange") %in% names(filt)))
  stop("Filtered DE must contain columns: padj, log2FoldChange")

padj <- as.numeric(filt$padj)
lfc  <- as.numeric(filt$log2FoldChange)
sig  <- !is.na(padj) & padj <= opt$`padj-threshold`

ids_all0  <- ids_filt[sig]
ids_up0   <- ids_filt[sig & lfc >=  opt$`logfc-threshold`]
ids_down0 <- ids_filt[sig & lfc <= -opt$`logfc-threshold`]

# universe = all tested (non-NA padj) from RAW, fallback to filtered
univ0 <- if ("padj" %in% names(raw)) ids_raw[!is.na(raw$padj)] else ids_filt[!is.na(filt$padj)]

writeLines(sprintf("Rows in sets (ALL/UP/DOWN/UNIV raw ids): %d / %d / %d / %d",
                   length(ids_all0), length(ids_up0), length(ids_down0), length(univ0)),
           file.path(opt$outdir,"debug_summary.txt"))

# ---------------- Map to ENSG ----------------
ids_all_ENSG  <- to_gene_ensembl(ids_all0)
ids_up_ENSG   <- to_gene_ensembl(ids_up0)
ids_down_ENSG <- to_gene_ensembl(ids_down0)
univ_ENSG     <- to_gene_ensembl(univ0)

# ---------------- Optional PC filter ----------------
pc <- character(0)
pc_applied <- FALSE
if (!is.na(opt$`pc-list`) && nzchar(opt$`pc-list`) && file.exists(opt$`pc-list`)) {
  pc <- unique(strip_version(readLines(opt$`pc-list`)))
  ids_all_ENSG  <- intersect(ids_all_ENSG,  pc)
  ids_up_ENSG   <- intersect(ids_up_ENSG,   pc)
  ids_down_ENSG <- intersect(ids_down_ENSG, pc)
  univ_ENSG     <- intersect(univ_ENSG,     pc)
  pc_applied <- TRUE
}

# ---------------- ENSG -> ENTREZ ----------------
map_ALL  <- ensg_to_entrez(ids_all_ENSG)
map_UP   <- ensg_to_entrez(ids_up_ENSG)
map_DOWN <- ensg_to_entrez(ids_down_ENSG)
map_UNIV <- ensg_to_entrez(univ_ENSG)

fwrite(map_ALL,  file.path(opt$outdir,"debug_map_ALL.tsv"),  sep="\t")
fwrite(map_UP,   file.path(opt$outdir,"debug_map_UP.tsv"),   sep="\t")
fwrite(map_DOWN, file.path(opt$outdir,"debug_map_DOWN.tsv"), sep="\t")
fwrite(map_UNIV, file.path(opt$outdir,"debug_map_BASE.tsv"), sep="\t")

# Summaries
dbg1 <- data.frame(
  Stage    = c("POST_MAP_UNIV","POST_MAP_ALL","POST_MAP_UP","POST_MAP_DOWN"),
  N_ENSG   = c(length(univ_ENSG), length(ids_all_ENSG), length(ids_up_ENSG), length(ids_down_ENSG)),
  N_ENTREZ = c(length(unique(map_UNIV$ENTREZID)),
               length(unique(map_ALL$ENTREZID)),
               length(unique(map_UP$ENTREZID)),
               length(unique(map_DOWN$ENTREZID))),
  PC_filter = if (pc_applied) "applied" else "not_applied"
)
fwrite(dbg1, file.path(opt$outdir, "debug_pc_filter_summary.tsv"), sep="\t")

pc_total <- length(pc)
dbg2 <- data.frame(
  stage     = c("ALL","UP","DOWN","UNIV"),
  before    = c(length(ids_all0), length(ids_up0), length(ids_down0), length(univ0)),
  after_ora = c(length(unique(map_ALL$ENTREZID)),
                length(unique(map_UP$ENTREZID)),
                length(unique(map_DOWN$ENTREZID)),
                length(unique(map_UNIV$ENTREZID))),
  pc_total  = pc_total,
  mode      = "ora"
)
fwrite(dbg2, file.path(opt$outdir,"debug_pc_filter.tsv"), sep="\t")

# ---------------- ORA (ENTREZ) ----------------
minGS <- opt$`min-gs-size`

nnz <- function(x) sum(!is.na(x) & x != "")
clean_ids <- function(x) {
  x <- unique(as.character(x))
  x[!is.na(x) & x != ""]
}

ent_all  <- clean_ids(map_ALL$ENTREZID)
ent_up   <- clean_ids(map_UP$ENTREZID)
ent_down <- clean_ids(map_DOWN$ENTREZID)
ent_univ <- clean_ids(map_UNIV$ENTREZID)

# enforce set ⊆ universe
ent_all  <- intersect(ent_all,  ent_univ)
ent_up   <- intersect(ent_up,   ent_univ)
ent_down <- intersect(ent_down, ent_univ)

ok_set <- function(set, base) nnz(set) >= minGS && nnz(base) >= 50

skip_msgs <- c()
do_or_skip <- function(set, base, fun, label, ...) {
  if (!ok_set(set, base)) {
    skip_msgs <<- c(skip_msgs, sprintf("%s skipped: set or universe too small (set=%d, univ=%d)",
                                       label, nnz(set), nnz(base)))
    return(NULL)
  }
  suppressMessages(fun(...))
}

# STRICT (for CSV + primary plots when present)
go_one_strict <- function(genes, universe, ont, label) {
  do_or_skip(genes, universe, enrichGO, label,
             gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = ont, universe = universe, pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             minGSSize = minGS, readable = FALSE)
}
kegg_strict <- function(genes, universe, label)
  do_or_skip(genes, universe, enrichKEGG, label,
             gene=genes, universe=universe, organism=opt$kegg_org,
             pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
             minGSSize=minGS)
reactome_strict <- function(genes, universe, label)
  do_or_skip(genes, universe, enrichPathway, label,
             gene=genes, universe=universe, organism=opt$organism,
             pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
             minGSSize=minGS, readable=TRUE)

onts <- c("BP","CC","MF")
res_go <- setNames(vector("list", length(onts)), onts)
for (ont in onts) {
  res_go[[ont]] <- list(
    ALL  = go_one_strict(ent_all,  ent_univ, ont, paste0("GO ORA ",ont," ALL")),
    UP   = go_one_strict(ent_up,   ent_univ, ont, paste0("GO ORA ",ont," UP")),
    DOWN = go_one_strict(ent_down, ent_univ, ont, paste0("GO ORA ",ont," DOWN"))
  )
}
ekg_ALL  <- kegg_strict(ent_all,  ent_univ, "KEGG ORA ALL")
ekg_UP   <- kegg_strict(ent_up,   ent_univ, "KEGG ORA UP")
ekg_DOWN <- kegg_strict(ent_down, ent_univ, "KEGG ORA DOWN")
erp_ALL  <- reactome_strict(ent_all,  ent_univ, "Reactome ORA ALL")
erp_UP   <- reactome_strict(ent_up,   ent_univ, "Reactome ORA UP")
erp_DOWN <- reactome_strict(ent_down, ent_univ, "Reactome ORA DOWN")

if (length(skip_msgs)) writeLines(skip_msgs, file.path(opt$outdir,"debug_skips.txt"))

# ---- Write ORA CSVs (STRICT ONLY) ----
for (ont in onts) {
  save_enrich_csv(res_go[[ont]]$ALL,  file.path(opt$outdir, sprintf("ORA_GO_%s_ALL_enrichment.csv",  ont)))
  save_enrich_csv(res_go[[ont]]$UP,   file.path(opt$outdir, sprintf("ORA_GO_%s_UP_enrichment.csv",   ont)))
  save_enrich_csv(res_go[[ont]]$DOWN, file.path(opt$outdir, sprintf("ORA_GO_%s_DOWN_enrichment.csv", ont)))
}
save_enrich_csv(ekg_ALL,  file.path(opt$outdir,"ORA_KEGG_ALL_enrichment.csv"))
save_enrich_csv(ekg_UP,   file.path(opt$outdir,"ORA_KEGG_UP_enrichment.csv"))
save_enrich_csv(ekg_DOWN, file.path(opt$outdir,"ORA_KEGG_DOWN_enrichment.csv"))
save_enrich_csv(erp_ALL,  file.path(opt$outdir,"ORA_Reactome_ALL_enrichment.csv"))
save_enrich_csv(erp_UP,   file.path(opt$outdir,"ORA_Reactome_UP_enrichment.csv"))
save_enrich_csv(erp_DOWN, file.path(opt$outdir,"ORA_Reactome_DOWN_enrichment.csv"))

# ---- ORA fallback (plotting only) ----
go_fallback <- function(genes, universe, ont) {
  do_or_skip(genes, universe, enrichGO, paste0("GO ORA ",ont," fallback"),
             gene=genes, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
             ont=ont, universe=universe, pAdjustMethod="BH",
             pvalueCutoff=1, qvalueCutoff=1,
             minGSSize=minGS, maxGSSize=20000, readable=FALSE)
}
kegg_fallback <- function(genes, universe)
  do_or_skip(genes, universe, enrichKEGG, "KEGG ORA fallback",
             gene=genes, universe=universe, organism=opt$kegg_org,
             pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
             minGSSize=minGS)
reactome_fallback <- function(genes, universe)
  do_or_skip(genes, universe, enrichPathway, "Reactome ORA fallback",
             gene=genes, universe=universe, organism=opt$organism,
             pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1,
             minGSSize=minGS, readable=TRUE)

# Utility: build per-panel subtitle that shows ACTUAL thresholds + sizes
panel_subtitle <- function(mode=c("STRICT","FALLBACK"), set_size, univ_size, minGS, extra=NULL){
  mode <- match.arg(mode)
  if (mode=="STRICT") {
    base <- sprintf("STRICT: p.adjust≤%.2g, q≤%.2g, minGS=%d | set=%d, univ=%d",
                    0.05, 0.20, minGS, set_size, univ_size)
  } else {
    base <- sprintf("FALLBACK (viz only): p≤1, q≤1, maxGS=20000, minGS=%d | set=%d, univ=%d",
                    minGS, set_size, univ_size)
  }
  if (!is.null(extra) && nzchar(extra)) base <- paste0(base, " | ", extra)
  base
}

# track which panels used fallback
fb_log <- data.frame(collection=character(), subset=character(), mode=character(),
                     set_size=integer(), universe=integer(), stringsAsFactors=FALSE)

# GO plot objects + subtitles
get_set_vec <- function(which) switch(which, ALL=ent_all, UP=ent_up, DOWN=ent_down)
plot_go <- list(BP=list(), CC=list(), MF=list())
sub_go  <- list(BP=list(), CC=list(), MF=list())

for (ont in onts) {
  for (nm in c("ALL","UP","DOWN")) {
    strict_obj <- res_go[[ont]][[nm]]
    genes <- get_set_vec(nm)
    if (!is.null(strict_obj) && nrow(as.data.frame(strict_obj)) > 0) {
      plot_go[[ont]][[nm]] <- strict_obj
      sub_go[[ont]][[nm]]  <- panel_subtitle("STRICT", set_size=nnz(genes), univ_size=nnz(ent_univ), minGS=minGS)
      fb_log <- rbind(fb_log, data.frame(collection=paste0("GO_",ont), subset=nm,
                                         mode="STRICT", set_size=nnz(genes), universe=nnz(ent_univ)))
    } else {
      fb_obj <- go_fallback(genes, ent_univ, ont)
      plot_go[[ont]][[nm]] <- fb_obj
      sub_go[[ont]][[nm]]  <- panel_subtitle("FALLBACK", set_size=nnz(genes), univ_size=nnz(ent_univ), minGS=minGS)
      fb_log <- rbind(fb_log, data.frame(collection=paste0("GO_",ont), subset=nm,
                                         mode="FALLBACK", set_size=nnz(genes), universe=nnz(ent_univ)))
    }
  }
}

# KEGG/Reactome plot objects + subtitles
use_or_fallback <- function(obj_strict, fb_fun, set_vec, coll_label){
  if (!is.null(obj_strict) && nrow(as.data.frame(obj_strict))>0) {
    list(obj=obj_strict,
         sub=panel_subtitle("STRICT", set_size=nnz(set_vec), univ_size=nnz(ent_univ), minGS=minGS),
         mode="STRICT", coll=coll_label)
  } else {
    fb <- fb_fun()
    list(obj=fb,
         sub=panel_subtitle("FALLBACK", set_size=nnz(set_vec), univ_size=nnz(ent_univ), minGS=minGS),
         mode="FALLBACK", coll=coll_label)
  }
}

plot_kegg <- list(
  ALL  = use_or_fallback(ekg_ALL,  function() kegg_fallback(ent_all,  ent_univ), ent_all,  "KEGG"),
  UP   = use_or_fallback(ekg_UP,   function() kegg_fallback(ent_up,   ent_univ), ent_up,   "KEGG"),
  DOWN = use_or_fallback(ekg_DOWN, function() kegg_fallback(ent_down, ent_univ), ent_down, "KEGG")
)
plot_react <- list(
  ALL  = use_or_fallback(erp_ALL,  function() reactome_fallback(ent_all,  ent_univ), ent_all,  "Reactome"),
  UP   = use_or_fallback(erp_UP,   function() reactome_fallback(ent_up,   ent_univ), ent_up,   "Reactome"),
  DOWN = use_or_fallback(erp_DOWN, function() reactome_fallback(ent_down, ent_univ), ent_down, "Reactome")
)

for (nm in c("ALL","UP","DOWN")) {
  fb_log <- rbind(fb_log,
                  data.frame(collection="KEGG", subset=nm, mode=plot_kegg[[nm]]$mode,
                             set_size=nnz(get_set_vec(nm)), universe=nnz(ent_univ)))
  fb_log <- rbind(fb_log,
                  data.frame(collection="Reactome", subset=nm, mode=plot_react[[nm]]$mode,
                             set_size=nnz(get_set_vec(nm)), universe=nnz(ent_univ)))
}

# Write fallback/strict usage log
fwrite(fb_log, file.path(opt$outdir,"debug_fallback_panels.tsv"), sep="\t")

# ---- ORA PDFs (titles + per-panel subtitles) ----
pdf(file.path(opt$outdir,"ORA_GO_dotplot_ALL_UP_DOWN_BPCCMF.pdf"), width=8, height=6)
for (ont in c("BP","CC","MF")) {
  for (nm in c("ALL","UP","DOWN")) {
    obj  <- plot_go[[ont]][[nm]]
    subt <- sub_go[[ont]][[nm]]
    print(dotplot_wrapped(obj, paste0("GO ORA [",ont,"] (",nm,")"), subt))
  }
}
dev.off()

pdf(file.path(opt$outdir,"ORA_KEGG_dotplot.pdf"), width=8, height=6)
for (nm in c("ALL","UP","DOWN")) {
  print(dotplot_wrapped(plot_kegg[[nm]]$obj, paste0("KEGG ORA (",nm,")"), plot_kegg[[nm]]$sub))
}
dev.off()

pdf(file.path(opt$outdir,"ORA_Reactome_dotplot.pdf"), width=8, height=6)
for (nm in c("ALL","UP","DOWN")) {
  print(dotplot_wrapped(plot_react[[nm]]$obj, paste0("Reactome ORA (",nm,")"), plot_react[[nm]]$sub))
}
dev.off()

# ---------------- GSEA (optional) ----------------
nterms <- function(x) if (is.null(x)) 0 else nrow(as.data.frame(x))

if (opt$`do-gsea`) {
  p <- if ("pvalue" %in% names(raw)) raw$pvalue else raw$padj
  p[is.na(p)] <- 1
  p <- pmin(p, 1e-300)
  lfc <- as.numeric(raw$log2FoldChange); lfc[is.na(lfc)] <- 0
  scores <- sign(lfc) * (-log10(p))
  names(scores) <- strip_version(raw[[1L]])

  dt_r <- data.table(id = names(scores), score = as.numeric(scores))
  dt_r[, ensg := map_ids_vector_to_ensg(id)]
  dt_r <- dt_r[!is.na(ensg) & grepl("^ENSG", ensg)]
  if (length(pc)) dt_r <- dt_r[ensg %in% pc]

  ranks_go <- if (nrow(dt_r)) {
    gdt <- dt_r[, .(score = max(score, na.rm = TRUE)), by = ensg][order(-score)]
    r <- gdt$score; names(r) <- gdt$ensg; r
  } else numeric()

  ranks_er <- if (length(ranks_go)) {
    m <- ensg_to_entrez(names(ranks_go))
    if (nrow(m)) {
      mm <- data.table(m)
      mm[, score := ranks_go[ENSEMBL]]
      mm <- mm[!is.na(score) & !is.na(ENTREZID) & ENTREZID != ""]
      erdt <- mm[, .(score = max(score, na.rm = TRUE)), by = ENTREZID][order(-score)]
      r <- erdt$score; names(r) <- erdt$ENTREZID; r
    } else numeric()
  } else numeric()

  onts <- c("BP","CC","MF")
  g_go <- setNames(vector("list", length(onts)), onts)
  for (ont in onts) {
    g_go[[ont]] <- if (length(ranks_go)) tryCatch(
      gseGO(geneList       = sort(ranks_go, decreasing = TRUE),
            OrgDb          = org.Hs.eg.db,
            keyType        = "ENSEMBL",
            ont            = ont,
            minGSSize      = opt$`min-gs-size`,
            pAdjustMethod  = "BH",
            pvalueCutoff   = opt$`gsea-q-cutoff`,
            verbose        = FALSE),
      error = function(e) NULL) else NULL
    save_enrich_csv(g_go[[ont]], file.path(opt$outdir, sprintf("GSEA_GO_%s_ALL_enrichment.csv", ont)), "GSEA")
  }

  g_keg <- if (length(ranks_er)) tryCatch(
    gseKEGG(geneList     = sort(ranks_er, decreasing = TRUE),
            organism     = opt$kegg_org,
            minGSSize    = opt$`min-gs-size`,
            pAdjustMethod= "BH",
            pvalueCutoff = opt$`gsea-q-cutoff`,
            verbose      = FALSE),
    error = function(e) NULL) else NULL

  g_rea <- if (length(ranks_er)) tryCatch(
    gsePathway(geneList   = sort(ranks_er, decreasing = TRUE),
               organism   = opt$organism,
               minGSSize  = opt$`min-gs-size`,
               pAdjustMethod = "BH",
               pvalueCutoff  = opt$`gsea-q-cutoff`,
               verbose    = FALSE),
    error = function(e) NULL) else NULL

  gsea_sub <- sprintf("GSEA STRICT: q≤%.2g, minGS=%d", opt$`gsea-q-cutoff`, opt$`min-gs-size`)

  pdf(file.path(opt$outdir, "GSEA_GO_dotplot_BPCMF.pdf"), width = 8, height = 6)
  for (ont in onts) print(dotplot_wrapped(g_go[[ont]], paste0("GO GSEA (", ont, ")"), gsea_sub))
  dev.off()

  pdf(file.path(opt$outdir, "GSEA_KEGG_dotplot.pdf"), width = 8, height = 6)
  print(dotplot_wrapped(g_keg, "KEGG GSEA (ALL)", gsea_sub))
  dev.off()

  pdf(file.path(opt$outdir, "GSEA_Reactome_dotplot.pdf"), width = 8, height = 6)
  print(dotplot_wrapped(g_rea, "Reactome GSEA (ALL)", gsea_sub))
  dev.off()

  save_enrich_csv(g_keg, file.path(opt$outdir, "GSEA_KEGG_ALL_enrichment.csv"), "GSEA")
  save_enrich_csv(g_rea, file.path(opt$outdir, "GSEA_Reactome_ALL_enrichment.csv"), "GSEA")

  writeLines(sprintf(
    "[GSEA] GO (BP/CC/MF) = %d / %d / %d  |  KEGG=%d  Reactome=%d (q<=%.2f)",
    nterms(g_go$BP), nterms(g_go$CC), nterms(g_go$MF),
    nterms(g_keg), nterms(g_rea), opt$`gsea-q-cutoff`),
    file.path(opt$outdir, "debug_gsea_summary.txt"))
}
