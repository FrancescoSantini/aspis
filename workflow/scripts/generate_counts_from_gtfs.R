#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(future.apply)
  library(tidyr)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("--gtfs"), type = "character", help = "Comma-separated list of GTF files"),
  make_option(c("--phenodata"), type = "character", help = "Phenodata CSV with sample names"),
  make_option(c("--gene-count-output"), type = "character", help = "Path to write gene count matrix (TSV)"),
  make_option(c("--gene-metadata-output"), type = "character", help = "Path to write gene metadata CSV"),
  make_option(c("--transcript-count-output"), type = "character", help = "Path to write transcript count matrix (TSV)"),
  make_option(c("--transcript-metadata-output"), type = "character", help = "Path to write transcript metadata CSV"),
  make_option(c("--tmap"), type = "character", default = NULL, help = "Optional .tmap file from gffcompare"),
  make_option(c("--known-codes-strict"),  type = "character", default = "=,j",
              help = "Comma-separated class codes counting as Known (strict) [default: %default]"),
  make_option(c("--known-codes-lenient"), type = "character", default = "=,j,c,o",
              help = "Comma-separated class codes counting as Known (lenient) [default: %default]"),
  make_option(c("--gene-type-view"), type = "character", default = "strict",
              help = "Which classification to write into 'gene_type': 'strict' or 'lenient' [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

strict_codes  <- unique(strsplit(opt$`known-codes-strict`,  "\\s*,\\s*")[[1]])
lenient_codes <- unique(strsplit(opt$`known-codes-lenient`, "\\s*,\\s*")[[1]])
view <- tolower(opt$`gene-type-view`)
if (!view %in% c("strict","lenient")) stop("Invalid --gene-type-view: ", view)

cat("[INFO] Known codes (strict):  ", paste(strict_codes,  collapse = ","), "\n", sep = "")
cat("[INFO] Known codes (lenient): ", paste(lenient_codes, collapse = ","), "\n", sep = "")
cat("[INFO] gene_type view: ", view, "\n", sep = "")

# ---------- Helpers ----------
looks_ensg <- function(x) !is.na(x) & grepl("^ENSG\\d+", x)
`%||%` <- function(a, b) if (length(a) && !is.na(a) && a != "") a else b
mode_or_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return(NA_character_)
  tt <- sort(table(x), decreasing = TRUE)
  names(tt)[1]
}

# ---------- Read phenodata (sanity only) ----------
suppressWarnings( read_csv(opt$phenodata, show_col_types = FALSE) )

# ---------- Parse GTFs (parallel) ----------
gtf_paths <- strsplit(opt$gtfs, ",")[[1]]
plan(multisession, workers = min(length(gtf_paths), 4))

parse_gtf <- function(gtf_file) {
  sample_id   <- sub(".*/([^/]+)\\.gtf$", "\\1", gtf_file)
  bioproject  <- basename(dirname(dirname(gtf_file)))
  bam_path    <- file.path("results/sorted", bioproject, paste0(sample_id, "_sorted.bam"))
  sample_lab  <- sub("\\.bam$", "", sub("_sorted$", "", basename(bam_path)))

  message(sprintf("[INFO] Parsing: %s", gtf_file))
  gtf <- fread(cmd = paste("grep '\ttranscript\t'", shQuote(gtf_file)),
               sep = "\t", header = FALSE, data.table = FALSE)
  if (nrow(gtf) == 0) return(NULL)

  colnames(gtf) <- c("chr","source","feature","start","end","score","strand","frame","attribute")
  attr <- gtf$attribute

  tx_id    <- str_match(attr, 'transcript_id "([^"]+)"')[,2]
  g_id     <- str_match(attr, 'gene_id "([^"]+)"')[,2]
  ref_g_id <- str_match(attr, 'ref_gene_id "([^"]+)"')[,2]
  g_name   <- str_match(attr, 'ref_gene_name "([^"]+)"')[,2]
  cov      <- suppressWarnings(as.numeric(str_match(attr, 'cov "([^"]+)"')[,2]))

  g_id   <- ifelse(!is.na(ref_g_id) & ref_g_id != "", ref_g_id, g_id)
  g_name <- ifelse(!is.na(g_name)  & g_name  != "", g_name, NA)

  valid <- !is.na(tx_id) & !is.na(g_id) & !is.na(cov)
  if (!any(valid)) return(NULL)
  gtf   <- gtf[valid, ]
  tx_id <- tx_id[valid]; g_id <- g_id[valid]; g_name <- g_name[valid]; cov <- cov[valid]

  len <- gtf$end - gtf$start + 1
  cov[is.na(cov)] <- 0
  counts <- ceiling(cov * len / 75)

  tibble(
    sample      = sample_lab,
    gene_id     = g_id,
    gene_name   = g_name,
    transcript  = tx_id,
    counts      = counts,
    chr         = gtf$chr,
    start       = gtf$start,
    end         = gtf$end,
    strand      = gtf$strand
  )
}

parsed_dfs <- future_lapply(gtf_paths, parse_gtf)
parsed_dfs <- parsed_dfs[!sapply(parsed_dfs, is.null)]
stopifnot(length(parsed_dfs) > 0)
all_data <- bind_rows(parsed_dfs) %>% filter(!is.na(gene_id), !is.na(transcript))

# ---- NEW: normalize sample names globally ----
all_data <- all_data %>%
  mutate(sample = basename(sample),
         sample = sub("\\.bam$", "", sample),
         sample = sub("_sorted$", "", sample))

# ---------- GENE COUNT MATRIX ----------
coord_df <- all_data %>%
  group_by(Geneid = gene_id) %>%
  summarise(
    Chr = mode_or_na(chr),
    Start = min(start),
    End = max(end),
    Strand = mode_or_na(strand),
    Length = End - Start + 1,
    .groups = "drop"
  )

gene_matrix <- all_data %>%
  mutate(counts = as.numeric(counts)) %>%
  group_by(Geneid = gene_id, sample) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = counts, values_fill = 0L) %>%
  arrange(Geneid)

gene_matrix <- left_join(coord_df, gene_matrix, by = "Geneid") %>%
  select(Geneid, Chr, Start, End, Strand, Length, everything())

# ---------- TRANSCRIPT COUNT MATRIX ----------
tx_matrix <- all_data %>%
  mutate(counts = as.numeric(counts)) %>%
  group_by(transcript_id = transcript, sample) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = counts, values_fill = 0L) %>%
  arrange(transcript_id)

write_tsv(gene_matrix, opt$`gene-count-output`)
write_tsv(tx_matrix,   opt$`transcript-count-output`)
cat(sprintf("[INFO] Written %d genes to %s\n", nrow(gene_matrix), opt$`gene-count-output`))
cat(sprintf("[INFO] Written %d transcripts to %s\n", nrow(tx_matrix), opt$`transcript-count-output`))

# ---------- METADATA ----------
shared_base <- all_data %>%
  transmute(
    gene_id,
    gene_name,
    transcript_id = transcript
  ) %>%
  distinct()

swap_idx <- looks_ensg(shared_base$gene_name) & !looks_ensg(shared_base$gene_id)
if (any(swap_idx)) {
  tmp <- shared_base$gene_id[swap_idx]
  shared_base$gene_id[swap_idx]   <- shared_base$gene_name[swap_idx]
  shared_base$gene_name[swap_idx] <- tmp
}

gene_name_by_id <- shared_base %>%
  group_by(gene_id) %>%
  summarise(gene_name = mode_or_na(gene_name), .groups = "drop")

tx_meta <- shared_base %>%
  transmute(
    feature_id    = transcript_id,
    transcript_id = transcript_id,
    gene_id
  ) %>%
  left_join(gene_name_by_id, by = "gene_id")

gene_meta <- gene_name_by_id %>%
  transmute(feature_id = gene_id, gene_id, gene_name)

# ----- .tmap integration (class_code per transcript) -----
class_code_levels <- c("=", "j", "c", "o", "e", "i", "u", "x", "s", "p", "r", "y", "k")
if (!is.null(opt$tmap) && file.exists(opt$tmap)) {
  cat("[INFO] Reading .tmap: ", opt$tmap, "\n", sep = "")
  tmap <- read_tsv(opt$tmap, comment = "#", col_types = cols())
  if ("qry_id" %in% names(tmap)) names(tmap)[names(tmap) == "qry_id"] <- "q_id"

  if (all(c("q_id","class_code") %in% names(tmap))) {
    tmap_class <- tmap %>%
      select(transcript_id = q_id, class_code) %>%
      filter(!is.na(transcript_id), !is.na(class_code)) %>%
      group_by(transcript_id) %>%
      summarise(
        class_code = {
          cc <- unique(class_code)
          cc <- cc[cc %in% class_code_levels]
          if (length(cc) == 0) NA_character_ else cc[order(match(cc, class_code_levels))][1]
        },
        .groups = "drop"
      )
    tx_meta <- tx_meta %>% left_join(tmap_class, by = "transcript_id")
  } else {
    warning("[WARN] .tmap lacks q_id/class_code; skipping class_code assignment.")
  }
} else {
  cat("[INFO] No .tmap provided; class_code will be NA and fallbacks will apply.\n")
}

# ----- strict/lenient at transcript level -----
is_ensg_tx  <- looks_ensg(tx_meta$gene_id)
has_code_tx <- !is.na(tx_meta$class_code)

is_known_strict_tx  <- has_code_tx & (tx_meta$class_code %in% strict_codes)  & is_ensg_tx
is_known_lenient_tx <- has_code_tx & (tx_meta$class_code %in% lenient_codes) & is_ensg_tx

# fallback when no class_code: treat Ensembl genes as Known so downstream doesn’t break
is_known_strict_tx  <- ifelse(is.na(is_known_strict_tx),  looks_ensg(tx_meta$gene_id), is_known_strict_tx)
is_known_lenient_tx <- ifelse(is.na(is_known_lenient_tx), looks_ensg(tx_meta$gene_id), is_known_lenient_tx)

tx_meta$gene_type_strict  <- ifelse(is_known_strict_tx,  "Known", "Novel")
tx_meta$gene_type_lenient <- ifelse(is_known_lenient_tx, "Known", "Novel")
tx_meta <- tx_meta %>% mutate(gene_type = if (view == "strict") gene_type_strict else gene_type_lenient) %>%
  arrange(transcript_id)

# ----- gene-level aggregation (any Known -> Known) -----
gene_from_tx <- tx_meta %>%
  group_by(gene_id) %>%
  summarise(
    gene_type_strict  = ifelse(any(gene_type_strict  == "Known"), "Known", "Novel"),
    gene_type_lenient = ifelse(any(gene_type_lenient == "Known"), "Known", "Novel"),
    .groups = "drop"
  )

gene_meta <- gene_meta %>%
  left_join(gene_from_tx, by = "gene_id")

# fallback if still NA
fill_idx <- is.na(gene_meta$gene_type_strict)
gene_meta$gene_type_strict[fill_idx]  <- ifelse(looks_ensg(gene_meta$gene_id[fill_idx]), "Known", "Novel")
gene_meta$gene_type_lenient[fill_idx] <- gene_meta$gene_type_strict[fill_idx]
gene_meta$gene_type <- if (view == "strict") gene_meta$gene_type_strict else gene_meta$gene_type_lenient
gene_meta <- gene_meta %>% arrange(gene_id)

# ----- write CSV metadata -----
write_csv(gene_meta, opt$`gene-metadata-output`)
write_csv(tx_meta,   opt$`transcript-metadata-output`)

cat(sprintf("[INFO] Written %d gene annotations to %s\n", nrow(gene_meta), opt$`gene-metadata-output`))
cat(sprintf("[INFO] Written %d transcript annotations to %s\n", nrow(tx_meta),   opt$`transcript-metadata-output`))

# brief summary
cat("[INFO] Summary (gene_meta): strict=", sum(gene_meta$gene_type_strict=="Known", na.rm=TRUE), "/",
    nrow(gene_meta), " Known; lenient=", sum(gene_meta$gene_type_lenient=="Known", na.rm=TRUE), "\n", sep="")
cat("[INFO] Summary (tx_meta):   strict=", sum(tx_meta$gene_type_strict=="Known", na.rm=TRUE), "/",
    nrow(tx_meta), " Known; lenient=", sum(tx_meta$gene_type_lenient=="Known", na.rm=TRUE), "\n", sep="")
