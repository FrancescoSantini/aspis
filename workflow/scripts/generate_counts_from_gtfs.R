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

# ---------- CONSTANTS ----------
# Keep 75 to match your current pipeline assumption and avoid downstream shifts.
READ_LENGTH <- 75L

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
cat("[INFO] Using READ_LENGTH = ", READ_LENGTH, " bp for count derivation\n", sep = "")

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
if (!is.null(opt$phenodata) && file.exists(opt$phenodata)) {
  suppressWarnings(readr::read_csv(opt$phenodata, show_col_types = FALSE))
}

# ---------- Parse GTFs (parallel, robust) ----------
gtf_paths <- strsplit(opt$gtfs, ",", fixed = TRUE)[[1]]
gtf_paths <- trimws(gtf_paths)
gtf_paths <- gtf_paths[gtf_paths != ""]
if (!length(gtf_paths)) stop("No --gtfs provided")

future::plan(future::multisession, workers = min(length(gtf_paths), 4))

parse_gtf <- function(gtf_file) {
  message(sprintf("[INFO] Parsing: %s", gtf_file))
  # Derive a clean sample label from the gtf path
  sample_lab <- sub("\\.gtf$", "", basename(gtf_file))

  # Read the GTF; drop comments and keep only 'transcript' rows
  gtf <- tryCatch(
    data.table::fread(
      gtf_file,
      sep = "\t", header = FALSE, data.table = FALSE,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(gtf) || !nrow(gtf)) return(NULL)
  if (ncol(gtf) < 9) return(NULL)

  # Discard comment lines if fread kept them
  if (is.character(gtf[[1]])) {
    gtf <- gtf[!startsWith(gtf[[1]], "#"), , drop = FALSE]
  }
  if (!nrow(gtf)) return(NULL)

  colnames(gtf) <- c("chr","source","feature","start","end","score","strand","frame","attribute")
  gtf <- gtf[gtf$feature == "transcript", , drop = FALSE]
  if (!nrow(gtf)) return(NULL)

  attr <- gtf$attribute

  tx_id    <- stringr::str_match(attr, 'transcript_id "([^"]+)"')[,2]
  g_id     <- stringr::str_match(attr, 'gene_id "([^"]+)"')[,2]
  ref_g_id <- stringr::str_match(attr, 'ref_gene_id "([^"]+)"')[,2]
  g_name   <- stringr::str_match(attr, 'ref_gene_name "([^"]+)"')[,2]
  cov_raw  <- stringr::str_match(attr, 'cov "([^"]+)"')[,2]

  # Prefer ref_gene_id when present
  g_id <- ifelse(!is.na(ref_g_id) & ref_g_id != "", ref_g_id, g_id)
  # Gene name is optional
  g_name[is.na(g_name)] <- NA_character_

  # Coverage and length → counts
  cov_num <- suppressWarnings(as.numeric(cov_raw))
  cov_num[is.na(cov_num)] <- 0
  len <- as.numeric(gtf$end) - as.numeric(gtf$start) + 1L
  len[is.na(len) | len < 0] <- 0
  counts <- ceiling(cov_num * len / READ_LENGTH)
  counts[!is.finite(counts) | counts < 0] <- 0

  # Keep valid rows
  valid <- !is.na(tx_id) & !is.na(g_id)
  if (!any(valid)) return(NULL)

  tibble::tibble(
    sample      = sample_lab,
    gene_id     = g_id[valid],
    gene_name   = g_name[valid],
    transcript  = tx_id[valid],
    counts      = as.integer(counts[valid]),
    chr         = gtf$chr[valid],
    start       = as.integer(gtf$start[valid]),
    end         = as.integer(gtf$end[valid]),
    strand      = gtf$strand[valid]
  )
}

parsed_dfs <- future.apply::future_lapply(gtf_paths, parse_gtf)
parsed_dfs <- parsed_dfs[!vapply(parsed_dfs, is.null, logical(1))]
if (!length(parsed_dfs)) stop("No usable transcripts found in provided GTFs")

all_data <- dplyr::bind_rows(parsed_dfs) |>
  dplyr::filter(!is.na(gene_id), !is.na(transcript))

# Normalize sample names coherently (no _sorted/.bam)
all_data <- all_data |>
  dplyr::mutate(
    sample = basename(sample),
    sample = sub("\\.bam$", "", sample),
    sample = sub("_sorted$", "", sample)
  )

# ---------- GENE COUNT MATRIX ----------
mode_or_na <- get("mode_or_na", inherits = FALSE)

coord_df <- all_data |>
  dplyr::group_by(Geneid = gene_id) |>
  dplyr::summarise(
    Chr    = mode_or_na(chr),
    Start  = min(start, na.rm = TRUE),
    End    = max(end,   na.rm = TRUE),
    Strand = mode_or_na(strand),
    Length = End - Start + 1L,
    .groups = "drop"
  )

gene_matrix <- all_data |>
  dplyr::mutate(counts = as.numeric(counts)) |>
  dplyr::group_by(Geneid = gene_id, sample) |>
  dplyr::summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = sample, values_from = counts, values_fill = 0) |>
  dplyr::arrange(Geneid)

gene_matrix <- dplyr::left_join(coord_df, gene_matrix, by = "Geneid") |>
  dplyr::select(Geneid, Chr, Start, End, Strand, Length, dplyr::everything())

# ---------- TRANSCRIPT COUNT MATRIX ----------
tx_matrix <- all_data |>
  dplyr::mutate(counts = as.numeric(counts)) |>
  dplyr::group_by(transcript_id = transcript, sample) |>
  dplyr::summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = sample, values_from = counts, values_fill = 0) |>
  dplyr::arrange(transcript_id)

readr::write_tsv(gene_matrix, opt$`gene-count-output`)
readr::write_tsv(tx_matrix,   opt$`transcript-count-output`)
cat(sprintf("[INFO] Written %d genes to %s\n", nrow(gene_matrix), opt$`gene-count-output`))
cat(sprintf("[INFO] Written %d transcripts to %s\n", nrow(tx_matrix), opt$`transcript-count-output`))

# ---------- METADATA ----------
shared_base <- all_data |>
  dplyr::transmute(
    gene_id,
    gene_name,
    transcript_id = transcript
  ) |>
  dplyr::distinct()

# If ref_gene_name accidentally carried an ENSG and gene_id didn't, swap
looks_ensg <- get("looks_ensg", inherits = FALSE)
swap_idx <- looks_ensg(shared_base$gene_name) & !looks_ensg(shared_base$gene_id)
if (any(swap_idx, na.rm = TRUE)) {
  tmp <- shared_base$gene_id[swap_idx]
  shared_base$gene_id[swap_idx]   <- shared_base$gene_name[swap_idx]
  shared_base$gene_name[swap_idx] <- tmp
}

gene_name_by_id <- shared_base |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(gene_name = mode_or_na(gene_name), .groups = "drop")

tx_meta <- shared_base |>
  dplyr::transmute(
    feature_id    = transcript_id,
    transcript_id = transcript_id,
    gene_id
  ) |>
  dplyr::left_join(gene_name_by_id, by = "gene_id")

gene_meta <- gene_name_by_id |>
  dplyr::transmute(feature_id = gene_id, gene_id, gene_name)

# ----- .tmap integration (class_code per transcript) -----
class_code_levels <- c("=", "j", "c", "o", "e", "i", "u", "x", "s", "p", "r", "y", "k")
if (!is.null(opt$tmap) && file.exists(opt$tmap)) {
  cat("[INFO] Reading .tmap: ", opt$tmap, "\n", sep = "")
  tmap <- suppressWarnings(readr::read_tsv(opt$tmap, comment = "#", show_col_types = FALSE))
  if ("qry_id" %in% names(tmap)) names(tmap)[names(tmap) == "qry_id"] <- "q_id"

  if (all(c("q_id","class_code") %in% names(tmap))) {
    tmap_class <- tmap |>
      dplyr::select(transcript_id = q_id, class_code) |>
      dplyr::filter(!is.na(transcript_id), !is.na(class_code)) |>
      dplyr::group_by(transcript_id) |>
      dplyr::summarise(
        class_code = {
          cc <- unique(class_code)
          cc <- cc[cc %in% class_code_levels]
          if (length(cc) == 0) NA_character_ else cc[order(match(cc, class_code_levels))][1]
        },
        .groups = "drop"
      )
    tx_meta <- dplyr::left_join(tx_meta, tmap_class, by = "transcript_id")
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
tx_meta <- tx_meta |>
  dplyr::mutate(gene_type = if (view == "strict") gene_type_strict else gene_type_lenient) |>
  dplyr::arrange(transcript_id)

# ----- gene-level aggregation (any Known -> Known) -----
gene_from_tx <- tx_meta |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    gene_type_strict  = ifelse(any(gene_type_strict  == "Known", na.rm = TRUE), "Known", "Novel"),
    gene_type_lenient = ifelse(any(gene_type_lenient == "Known", na.rm = TRUE), "Known", "Novel"),
    .groups = "drop"
  )

gene_meta <- gene_meta |>
  dplyr::left_join(gene_from_tx, by = "gene_id")

# fallback if still NA
fill_idx <- is.na(gene_meta$gene_type_strict)
gene_meta$gene_type_strict[fill_idx]  <- ifelse(looks_ensg(gene_meta$gene_id[fill_idx]), "Known", "Novel")
gene_meta$gene_type_lenient[fill_idx] <- gene_meta$gene_type_strict[fill_idx]
gene_meta$gene_type <- if (view == "strict") gene_meta$gene_type_strict else gene_meta$gene_type_lenient
gene_meta <- gene_meta |>
  dplyr::arrange(gene_id)

# ----- write CSV metadata -----
readr::write_csv(gene_meta, opt$`gene-metadata-output`)
readr::write_csv(tx_meta,   opt$`transcript-metadata-output`)

cat(sprintf("[INFO] Written %d gene annotations to %s\n", nrow(gene_meta), opt$`gene-metadata-output`))
cat(sprintf("[INFO] Written %d transcript annotations to %s\n", nrow(tx_meta),   opt$`transcript-metadata-output`))

# brief summary
cat("[INFO] Summary (gene_meta): strict=", sum(gene_meta$gene_type_strict=="Known", na.rm=TRUE), "/",
    nrow(gene_meta), " Known; lenient=", sum(gene_meta$gene_type_lenient=="Known", na.rm=TRUE), "\n", sep="")
cat("[INFO] Summary (tx_meta):   strict=", sum(tx_meta$gene_type_strict=="Known", na.rm=TRUE), "/",
    nrow(tx_meta), " Known; lenient=", sum(tx_meta$gene_type_lenient=="Known", na.rm=TRUE), "\n", sep="")
