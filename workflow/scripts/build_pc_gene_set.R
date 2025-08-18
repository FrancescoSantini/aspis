#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
})

opt_list <- list(
  make_option("--gtf",    type="character", help="Reference GTF (from config['reference_annotation'])"),
  make_option("--outdir", type="character", help="Output dir, e.g. results/deg/<bioproject>"),
  make_option("--pc-labels", type="character", default="protein coding",
              help="Biotype labels to treat as protein-coding (comma-separated, case-insensitive). Default: 'protein coding'")
)
opt <- parse_args(OptionParser(option_list = opt_list))
stopifnot(!is.null(opt$gtf), file.exists(opt$gtf), !is.null(opt$outdir))
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# -------- helpers --------
extract_attr <- function(x, key){
  m <- str_match(x, paste0('(?:^|;)[[:space:]]*', key, '[[:space:]]+"([^"]+)"'))
  m[,2]
}
norm_biotype <- function(x){
  x <- tolower(x)
  x <- gsub("[-_]", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

# -------- read GTF robustly (strip comment lines first) --------
all_lines <- readLines(opt$gtf, warn = FALSE)
data_lines <- all_lines[!grepl("^#", all_lines)]
if (!length(data_lines)) stop("No non-comment lines found in GTF.")
gtf <- fread(text = data_lines, sep = "\t", header = FALSE, data.table = TRUE,
             quote = "", fill = TRUE, showProgress = FALSE)
if (ncol(gtf) < 9) stop("Parsed GTF has fewer than 9 columns after comment removal.")
setnames(gtf, paste0("V", 1:9))
gtf <- gtf[V3 == "gene"]                        # only gene features

attr <- gtf$V9
gid  <- extract_attr(attr, "gene_id")
bt   <- extract_attr(attr, "gene_biotype")
if (all(is.na(bt))) bt <- extract_attr(attr, "gene_type")  # fallback

dt <- data.table(gene_id = gid, biotype = norm_biotype(bt))
dt <- dt[!is.na(gene_id) & gene_id != ""]

# write full table
biotype_tsv <- file.path(opt$outdir, "ref_gene_biotype.tsv")
fwrite(unique(dt), biotype_tsv, sep = "\t")

# decide which labels are "protein-coding"
pc_labels <- norm_biotype(unlist(strsplit(opt$`pc-labels`, ",")))
is_pc <- dt$biotype %in% pc_labels
pc <- unique(dt$gene_id[is_pc])

pc_txt <- file.path(opt$outdir, "pc_gene_ids.txt")
writeLines(pc, pc_txt)

# quick summary
sum_txt <- file.path(opt$outdir, "ref_gene_biotype.summary.txt")
counts <- dt[, .N, by = biotype][order(-N)]
sink(sum_txt)
cat(sprintf("GTF: %s\nGenes: %d\nUnique biotypes: %d\n\n", opt$gtf, nrow(dt), nrow(counts)))
print(counts)
cat("\nProtein-coding labels considered:\n  - ", paste(pc_labels, collapse = "\n  - "), "\n", sep = "")
cat(sprintf("\nProtein-coding gene_id count: %d\n", length(pc)))
sink()
