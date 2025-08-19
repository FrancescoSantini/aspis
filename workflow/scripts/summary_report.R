#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(data.table)
})

# ----------------- CLI -----------------
option_list <- list(
  make_option(c("-o","--output"), type="character"),
  make_option("--title",  type="character", default="Summary Report"),
  make_option("--topn-gene", type="integer", default=30),
  make_option("--topn-tx",   type="integer", default=30),

  # optional metadata paths (recommended)
  make_option("--gene-meta", type="character", default=NA),
  make_option("--tx-meta",   type="character", default=NA),

  # Gene DE + plots
  make_option("--gene-filtered",     type="character"),
  make_option("--gene-pca-pdf",      type="character"),
  make_option("--gene-volcano-pdf",  type="character"),
  make_option("--gene-heatmap-pdf",  type="character"),
  # Gene ORA
  make_option("--gene-go-pdf",       type="character"),
  make_option("--gene-kegg-pdf",     type="character"),
  make_option("--gene-reactome-pdf", type="character"),
  # Gene GSEA (optional)
  make_option("--gene-gsea-go-pdf",       type="character", default=NA),
  make_option("--gene-gsea-kegg-pdf",     type="character", default=NA),
  make_option("--gene-gsea-reactome-pdf", type="character", default=NA),

  # Tx DE + plots
  make_option("--tx-filtered",     type="character"),
  make_option("--tx-pca-pdf",      type="character"),
  make_option("--tx-volcano-pdf",  type="character"),
  make_option("--tx-heatmap-pdf",  type="character"),
  # Tx ORA
  make_option("--tx-go-pdf",       type="character"),
  make_option("--tx-kegg-pdf",     type="character"),
  make_option("--tx-reactome-pdf", type="character"),
  # Tx GSEA (optional)
  make_option("--tx-gsea-go-pdf",       type="character", default=NA),
  make_option("--tx-gsea-kegg-pdf",     type="character", default=NA),
  make_option("--tx-gsea-reactome-pdf", type="character", default=NA),

  # Isoform switch
  make_option("--switch-qc-pdf", type="character")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ----------------- Basics -----------------
outdir <- dirname(opt$output)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

esc <- function(x) {
  x <- gsub("&","&amp;", x, fixed=TRUE)
  x <- gsub("<","&lt;", x, fixed=TRUE)
  x <- gsub(">","&gt;", x, fixed=TRUE)
  x <- gsub('"',"&quot;",x, fixed=TRUE)
  x
}
has_file <- function(p) is.character(p) && !is.na(p) && nzchar(p) && file.exists(p)
fail <- function(msg, also_html=TRUE, code=1){
  if (also_html) {
    cat("<!doctype html><meta charset='utf-8'><pre style='white-space:pre-wrap'>", esc(msg), "</pre>",
        file = opt$output)
  }
  quit(save="no", status=code)
}

# ----------------- PDF conversion -----------------
which_pdftoppm <- Sys.which("pdftoppm")
which_convert  <- Sys.which("convert")
if (which_pdftoppm == "" && which_convert == "") {
  fail("Neither 'pdftoppm' (poppler-utils) nor 'convert' (ImageMagick) found in PATH.\nInstall poppler-utils (recommended).")
}

# Convert a pdf to PNGs in subdir, return basenames
pdf_to_pngs <- function(pdf, subdir_stub, dpi = 220L) {
  if (!has_file(pdf)) return(character(0))
  subdir <- file.path(outdir, subdir_stub)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

  if (which_pdftoppm != "") {
    prefix <- file.path(subdir, subdir_stub)
    args <- c("-png", "-r", as.character(dpi), pdf, prefix)
    suppressWarnings(system2(which_pdftoppm, args, stdout = TRUE, stderr = TRUE))
    g <- Sys.glob(paste0(prefix, "-*.png"))
    if (length(g)) {
      ord <- order(as.integer(gsub(".*-(\\d+)\\.png$","\\1", g)))
      return(basename(normalizePath(g[ord], mustWork = FALSE)))
    }
  } else {
    outpat <- file.path(subdir, paste0(subdir_stub, "_%03d.png"))
    args <- c("-density", as.character(dpi), pdf, "-quality", "92", outpat)
    suppressWarnings(system2(which_convert, args, stdout = TRUE, stderr = TRUE))
    g <- Sys.glob(file.path(subdir, paste0(subdir_stub, "_*.png")))
    if (length(g)) {
      ord <- order(as.integer(gsub(".*_(\\d+)\\.png$","\\1", g)))
      return(basename(normalizePath(g[ord], mustWork = FALSE)))
    }
  }
  character(0)
}

# ----------------- Metadata mapping (length-safe) -----------------
read_meta <- function(p) {
  if (!has_file(p)) return(NULL)
  df <- suppressMessages(tryCatch(readr::read_csv(p, show_col_types = FALSE), error=function(e) NULL))
  if (is.null(df)) return(NULL)
  if (!"feature_id" %in% names(df)) df$feature_id <- NA_character_
  if (!"gene_id"    %in% names(df)) df$gene_id    <- NA_character_
  if (!"transcript_id" %in% names(df)) df$transcript_id <- NA_character_
  if (!"gene_name"  %in% names(df)) df$gene_name  <- NA_character_
  df <- df %>% mutate(across(c(feature_id,gene_id,transcript_id,gene_name), as.character))
  df
}

build_display_names <- function(ids, meta, mode=c("gene","transcript")) {
  mode <- match.arg(mode)
  disp <- ifelse(grepl("\\|", ids), sub(".*\\|", "", ids), NA_character_)
  if (!is.null(meta)) {
    if (mode == "gene") {
      if (all(c("feature_id","gene_name") %in% names(meta))) {
        m <- meta %>% filter(!is.na(feature_id), !is.na(gene_name), nzchar(feature_id), nzchar(gene_name)) %>%
          distinct(feature_id, .keep_all = TRUE)
        idx <- which(is.na(disp) | !nzchar(disp))
        if (length(idx)) {
          mm <- match(ids[idx], m$feature_id)
          take <- which(!is.na(mm)); if (length(take)) disp[idx[take]] <- m$gene_name[mm[take]]
        }
      }
      if (all(c("gene_id","gene_name") %in% names(meta))) {
        m <- meta %>% filter(!is.na(gene_id), !is.na(gene_name), nzchar(gene_id), nzchar(gene_name)) %>%
          distinct(gene_id, .keep_all = TRUE)
        idx <- which(is.na(disp) | !nzchar(disp))
        if (length(idx)) {
          mm <- match(ids[idx], m$gene_id)
          take <- which(!is.na(mm)); if (length(take)) disp[idx[take]] <- m$gene_name[mm[take]]
        }
      }
    } else {
      if (all(c("transcript_id","gene_name") %in% names(meta))) {
        m <- meta %>% filter(!is.na(transcript_id), !is.na(gene_name), nzchar(transcript_id), nzchar(gene_name)) %>%
          distinct(transcript_id, .keep_all = TRUE)
        idx <- which(is.na(disp) | !nzchar(disp))
        if (length(idx)) {
          mm <- match(ids[idx], m$transcript_id)
          take <- which(!is.na(mm)); if (length(take)) disp[idx[take]] <- m$gene_name[mm[take]]
        }
      }
      if (all(c("feature_id","gene_name") %in% names(meta))) {
        m <- meta %>% filter(!is.na(feature_id), !is.na(gene_name), nzchar(feature_id), nzchar(gene_name)) %>%
          distinct(feature_id, .keep_all = TRUE)
        idx <- which(is.na(disp) | !nzchar(disp))
        if (length(idx)) {
          mm <- match(ids[idx], m$feature_id)
          take <- which(!is.na(mm)); if (length(take)) disp[idx[take]] <- m$gene_name[mm[take]]
        }
      }
    }
  }
  disp[is.na(disp) | !nzchar(disp)] <- ids[is.na(disp) | !nzchar(disp)]
  disp
}

read_top_table <- function(csv_path, meta_path, mode=c("gene","transcript"), topn=30) {
  mode <- match.arg(mode)
  if (!has_file(csv_path)) return(NULL)
  df <- suppressMessages(readr::read_csv(csv_path, show_col_types = FALSE))
  if (!nrow(df)) return(NULL)
  if ("...1" %in% names(df)) names(df)[names(df) == "...1"] <- "id"
  if (!"id" %in% names(df)) names(df)[1] <- "id"
  keep <- c("id", intersect(c("log2FoldChange","padj"), names(df)))
  df <- df[, keep, drop=FALSE]
  df$log2FoldChange <- suppressWarnings(as.numeric(df$log2FoldChange))
  df$padj <- suppressWarnings(as.numeric(df$padj))
  meta <- read_meta(meta_path)
  df$display_name <- build_display_names(df$id, meta, mode = mode)
  df %>% arrange(is.na(padj), padj, desc(abs(log2FoldChange))) %>% slice_head(n = topn)
}

# ----------------- HTML skeleton -----------------
con <- file(opt$output, open="wt"); on.exit(close(con), add=TRUE)
cat("<!doctype html><html><head><meta charset='utf-8'>\n", file=con)
cat("<title>", esc(opt$title), "</title>\n", file=con, sep="")
cat("<style>
:root{
  --gap: 18px;
}

/* Use (almost) full screen width, but keep a sane hard cap */
body{
  font-family:system-ui,-apple-system,Segoe UI,Roboto,Ubuntu,Arial,sans-serif;
  margin:24px auto;
  width:calc(100vw - 48px);
  max-width:2200px;   /* gives you more room on big 16:9 displays */
}

/* Headings */
h1{margin:0 0 12px}
h2{margin:28px 0 10px;border-bottom:1px solid #eee;padding-bottom:6px}
h3{margin:18px 0 8px}

/* Tables centered and tidy */
.tablewrap{display:flex;justify-content:center}
table{border-collapse:collapse;margin:12px 0 24px;width:100%;max-width:1500px}
th,td{border:1px solid #ddd;padding:8px 10px;font-size:0.95rem}
th{background:#f7f7f7;text-align:left}

/* Tiles */
.grid{display:grid;gap:var(--gap);margin:8px 0 26px;justify-content:stretch}
.grid.grid3{grid-template-columns:repeat(3, 1fr)}   /* 3 equal columns fill the width */
.grid.grid2{grid-template-columns:repeat(2, 1fr)}   /* used for isoform tiles */

.rowpair{display:flex;flex-wrap:wrap;justify-content:center;gap:var(--gap);margin:8px 0 26px}
/* Make each of the 2 tiles in a pair exactly the width of a 3-col tile */
.rowpair .tile{flex:0 0 calc((100% - 2*var(--gap)) / 3)}

.tile{background:#fff;border:1px solid #ddd;border-radius:10px;padding:8px}
.tile img{width:100%;height:auto;display:block;border-radius:6px}

/* Hero images ~30% smaller than full width and centered */
.hero{width:50%;margin:10px auto 22px}
.hero img{width:100%;height:auto;border:1px solid #ddd;border-radius:8px}

/* Utility */
.muted{color:#888}
.centered{display:flex;justify-content:center}
</style>\n", file=con)
cat("</head><body>\n", file=con)
cat("<h1>", esc(opt$title), "</h1>\n", file=con, sep="")
cat("<div class='muted'>Generated: ", esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "</div>\n", file=con, sep="")

# helper to render a 3‑col grid of images
grid3 <- function(title, base, vec) {
  cat(sprintf("<h3>%s</h3>\n", esc(title)), file=con)
  if (!length(vec)) { cat("<div class='muted'>Not found.</div>", file=con); return() }
  cat("<div class='grid grid3'>\n", file=con)
  for (b in vec) cat("<div class='tile'><img src='", base, "/", esc(b), "'/></div>\n", file=con, sep="")
  cat("</div>\n", file=con)
}

# a centered pair, same tile width as 3‑col rows
pair_centered <- function(title, base_left, left_vec, base_right, right_vec) {
  cat(sprintf("<h3>%s</h3>\n", esc(title)), file=con)
  if (!length(left_vec) && !length(right_vec)) { cat("<div class='muted'>Not found.</div>", file=con); return() }
  cat("<div class='rowpair'>\n", file=con)
  if (length(left_vec))  cat("<div class='tile'><img src='",  base_left,  "/", esc(left_vec[1]),  "'/></div>\n",  file=con, sep="")
  if (length(right_vec)) cat("<div class='tile'><img src='", base_right, "/", esc(right_vec[1]), "'/></div>\n", file=con, sep="")
  cat("</div>\n", file=con)
}

# ----------------- Gene section -----------------
cat("<h2>Gene-level (DESeq2 / ORA / GSEA)</h2>\n", file=con)
gene_meta_guess <- if (has_file(opt$`gene-meta`)) opt$`gene-meta` else {
  if (has_file(opt$`gene-filtered`)) {
    f <- file.path(dirname(dirname(opt$`gene-filtered`)), "gene_metadata.csv")
    if (has_file(f)) f else NA_character_
  } else NA_character_
}
gtab <- read_top_table(opt$`gene-filtered`, gene_meta_guess, mode="gene", topn = opt$`topn-gene`)
if (!is.null(gtab) && nrow(gtab)) {
  cat("<div class='tablewrap'><table><tr><th>Gene</th><th>id</th><th>log2FC</th><th>padj</th></tr>\n", file=con)
  apply(gtab, 1, function(r) {
    cat(sprintf("<tr><td>%s</td><td><code>%s</code></td><td>%.3f</td><td>%s</td></tr>\n",
                esc(as.character(r[["display_name"]])),
                esc(as.character(r[["id"]])),
                as.numeric(r[["log2FoldChange"]]),
                ifelse(is.na(as.numeric(r[["padj"]])), "NA",
                       formatC(as.numeric(r[["padj"]]), format="e", digits=2))), file=con)
  })
  cat("</table></div>\n", file=con)
} else {
  cat("<p class='muted'>No gene-level DE rows found.</p>\n", file=con)
}

# gene plots -> pngs
gene_pca     <- pdf_to_pngs(opt$`gene-pca-pdf`,      "gene_PCA_hero", dpi = 220)
gene_volcano <- pdf_to_pngs(opt$`gene-volcano-pdf`,  "gene_volcano_all", dpi = 220)
gene_heatmap <- pdf_to_pngs(opt$`gene-heatmap-pdf`,  "gene_heatmap_all", dpi = 220)
gene_go      <- pdf_to_pngs(opt$`gene-go-pdf`,       "gene_ORA_GO", dpi = 220)
gene_kegg    <- pdf_to_pngs(opt$`gene-kegg-pdf`,     "gene_ORA_KEGG", dpi = 220)
gene_react   <- pdf_to_pngs(opt$`gene-reactome-pdf`, "gene_ORA_REACT", dpi = 220)
gene_gsea_go    <- pdf_to_pngs(opt$`gene-gsea-go-pdf`,       "gene_GSEA_GO", dpi = 220)
gene_gsea_kegg  <- pdf_to_pngs(opt$`gene-gsea-kegg-pdf`,     "gene_GSEA_KEGG", dpi = 220)
gene_gsea_react <- pdf_to_pngs(opt$`gene-gsea-reactome-pdf`, "gene_GSEA_REACT", dpi = 220)

# hero PCA (smaller)
if (length(gene_pca)) {
  cat("<div class='hero'><img src='gene_PCA_hero/", esc(gene_pca[1]), "'/></div>\n", file=con, sep="")
} else {
  cat("<p class='muted'>PCA (gene) not found.</p>\n", file=con)
}

# Volcano & Heatmap as 3‑per‑row grids (no sliders)
grid3("Volcano (Gene)",  "gene_volcano_all", gene_volcano)
grid3("Heatmap (Gene)",  "gene_heatmap_all", gene_heatmap)

# ORA tiles (3-per-row)
grid3("ORA – GO (Gene)",       "gene_ORA_GO",    gene_go)
grid3("ORA – KEGG (Gene)",     "gene_ORA_KEGG",  gene_kegg)
grid3("ORA – Reactome (Gene)", "gene_ORA_REACT", gene_react)

# GSEA: GO as 3-per-row, then KEGG+Reactome centered pair at SAME size
grid3("GSEA – GO (Gene)", "gene_GSEA_GO", gene_gsea_go)
pair_centered("GSEA – KEGG (Gene) – GSEA – Reactome (Gene)", "gene_GSEA_KEGG", gene_gsea_kegg, "gene_GSEA_REACT", gene_gsea_react)

# ----------------- Transcript section -----------------
cat("<h2>Transcript-level (DESeq2 / ORA / GSEA)</h2>\n", file=con)
tx_meta_guess <- if (has_file(opt$`tx-meta`)) opt$`tx-meta` else {
  if (has_file(opt$`tx-filtered`)) {
    f <- file.path(dirname(dirname(opt$`tx-filtered`)), "transcript_metadata.csv")
    if (has_file(f)) f else NA_character_
  } else NA_character_
}
ttab <- read_top_table(opt$`tx-filtered`, tx_meta_guess, mode="transcript", topn = opt$`topn-tx`)
if (!is.null(ttab) && nrow(ttab)) {
  cat("<div class='tablewrap'><table><tr><th>Transcript/Gene</th><th>id</th><th>log2FC</th><th>padj</th></tr>\n", file=con)
  apply(ttab, 1, function(r) {
    cat(sprintf("<tr><td>%s</td><td><code>%s</code></td><td>%.3f</td><td>%s</td></tr>\n",
                esc(as.character(r[["display_name"]])),
                esc(as.character(r[["id"]])),
                as.numeric(r[["log2FoldChange"]]),
                ifelse(is.na(as.numeric(r[["padj"]])), "NA",
                       formatC(as.numeric(r[["padj"]]), format="e", digits=2))), file=con)
  })
  cat("</table></div>\n", file=con)
} else {
  cat("<p class='muted'>No transcript-level DE rows found.</p>\n", file=con)
}

tx_pca     <- pdf_to_pngs(opt$`tx-pca-pdf`,      "tx_PCA_hero", dpi = 220)
tx_volcano <- pdf_to_pngs(opt$`tx-volcano-pdf`,  "tx_volcano_all", dpi = 220)
tx_heatmap <- pdf_to_pngs(opt$`tx-heatmap-pdf`,  "tx_heatmap_all", dpi = 220)
tx_go      <- pdf_to_pngs(opt$`tx-go-pdf`,       "tx_ORA_GO", dpi = 220)
tx_kegg    <- pdf_to_pngs(opt$`tx-kegg-pdf`,     "tx_ORA_KEGG", dpi = 220)
tx_react   <- pdf_to_pngs(opt$`tx-reactome-pdf`, "tx_ORA_REACT", dpi = 220)
tx_gsea_go    <- pdf_to_pngs(opt$`tx-gsea-go-pdf`,       "tx_GSEA_GO", dpi = 220)
tx_gsea_kegg  <- pdf_to_pngs(opt$`tx-gsea-kegg-pdf`,     "tx_GSEA_KEGG", dpi = 220)
tx_gsea_react <- pdf_to_pngs(opt$`tx-gsea-reactome-pdf`, "tx_GSEA_REACT", dpi = 220)

if (length(tx_pca)) {
  cat("<div class='hero'><img src='tx_PCA_hero/", esc(tx_pca[1]), "'/></div>\n", file=con, sep="")
} else cat("<p class='muted'>PCA (tx) not found.</p>\n", file=con)

grid3("Volcano (Transcript)", "tx_volcano_all", tx_volcano)
grid3("Heatmap (Transcript)", "tx_heatmap_all", tx_heatmap)

grid3("ORA – GO (Transcript)",       "tx_ORA_GO",    tx_go)
grid3("ORA – KEGG (Transcript)",     "tx_ORA_KEGG",  tx_kegg)
grid3("ORA – Reactome (Transcript)", "tx_ORA_REACT", tx_react)

grid3("GSEA – GO (Transcript)", "tx_GSEA_GO", tx_gsea_go)
pair_centered("GSEA – KEGG (Transcript) – GSEA – Reactome (Transcript)", "tx_GSEA_KEGG", tx_gsea_kegg, "tx_GSEA_REACT", tx_gsea_react)

# ----------------- Isoform switch -----------------
cat("<h2>Isoform switch</h2>\n", file=con)
sw <- pdf_to_pngs(opt$`switch-qc-pdf`, "switch_QC", dpi = 220)
if (length(sw)) {
  # First page: smaller hero
  cat("<div class='hero'><img src='switch_QC/", esc(sw[1]), "'/></div>\n", file=con, sep="")
  # Remaining pages as 2-per-row grid, centered by design
  if (length(sw) > 1) {
    cat("<div class='grid grid2'>\n", file=con)
    for (b in sw[-1]) cat("<div class='tile'><img src='switch_QC/", esc(b), "'/></div>\n", file=con, sep="")
    cat("</div>\n", file=con)
  }
} else {
  cat("<p class='muted'>Isoform switch QC PDF not found.</p>\n", file=con)
}

cat("</body></html>\n", file=con)
