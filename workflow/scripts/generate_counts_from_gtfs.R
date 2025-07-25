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

# --- CLI parsing ---
option_list <- list(
  make_option(c("--gtfs"), type = "character", help = "Comma-separated list of GTF files"),
  make_option(c("--phenodata"), type = "character", help = "Phenodata CSV with sample names"),
  make_option(c("--gene-output"), type = "character", help = "Path to write gene count matrix"),
  make_option(c("--transcript-output"), type = "character", help = "Path to write transcript count matrix")
)
opt <- parse_args(OptionParser(option_list = option_list))

# --- Read sample names ---
phenodata <- read_csv(opt$phenodata, show_col_types = FALSE)
sample_names <- rownames(read.csv(opt$phenodata, row.names = 1))

# --- Parse GTFs (in parallel) ---
gtf_paths <- strsplit(opt$gtfs, ",")[[1]]
plan(multisession, workers = max(1, availableCores() - 1))  # parallel setup

parse_gtf <- function(gtf_file) {
  sample_id <- sub(".*/([^/]+)\\.gtf$", "\\1", gtf_file)
  bioproject <- basename(dirname(dirname(gtf_file)))
  sample_bam_path <- file.path("results/sorted", bioproject, paste0(sample_id, "_sorted.bam"))

  message(sprintf("[INFO] Parsing: %s", gtf_file))

  gtf <- fread(cmd = paste("grep '\ttranscript\t'", gtf_file),
               sep = "\t", header = FALSE, data.table = FALSE)
  if (nrow(gtf) == 0) return(NULL)

  colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  attr <- gtf$attribute

  tx_id   <- str_match(attr, 'transcript_id "([^"]+)"')[,2]
  g_id    <- str_match(attr, 'gene_id "([^"]+)"')[,2]
  g_name  <- str_match(attr, 'ref_gene_name "([^"]+)"')[,2]
  cov     <- as.numeric(str_match(attr, 'cov "([^"]+)"')[,2])

  # Filter valid rows
  valid <- !is.na(tx_id) & !is.na(g_id) & !is.na(cov)
  gtf <- gtf[valid, ]
  tx_id <- tx_id[valid]; g_id <- g_id[valid]; g_name <- g_name[valid]; cov <- cov[valid]

  len <- gtf$end - gtf$start + 1
  counts <- ceiling(cov * len / 75)
  gene_key <- ifelse(g_name != "", paste0(g_id, "|", g_name), g_id)

  tibble(
    sample = sample_bam_path,
    gene = gene_key,
    transcript = tx_id,
    counts = counts,
    chr = gtf$chr,
    start = gtf$start,
    end = gtf$end,
    strand = gtf$strand
  )
}

parsed_dfs <- future_lapply(gtf_paths, parse_gtf)
parsed_dfs <- parsed_dfs[!sapply(parsed_dfs, is.null)]
all_data <- bind_rows(parsed_dfs)

# --- Remove rows with NA gene or transcript identifiers ---
all_data <- all_data %>%
  filter(!is.na(gene), !is.na(transcript))

# --- Gene matrix ---
gene_subset <- all_data %>%
  select(sample, gene, counts) %>%
  mutate(counts = as.numeric(counts))

gene_matrix <- gene_subset %>%
  mutate(counts = replace_na(counts, 0)) %>%
  group_by(gene, sample) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = counts, values_fill = 0) %>%
  rename(Geneid = gene)

# --- Transcript matrix ---
tx_subset <- all_data %>%
  select(sample, transcript, counts) %>%
  mutate(counts = as.numeric(counts))

tx_matrix <- tx_subset %>%
  mutate(counts = replace_na(counts, 0)) %>%
  group_by(transcript, sample) %>%
  summarise(counts = sum(counts), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = counts, values_fill = 0) %>%
  rename(transcript_id = transcript)

# --- Coordinates (gene-level) ---
coord_df <- all_data %>%
  group_by(Geneid = gene) %>%
  summarise(
    Chr = first(chr),
    Start = min(start),
    End = max(end),
    Strand = first(strand),
    Length = End - Start + 1,
    .groups = "drop"
  )

# --- Merge coordinates before select ---
gene_matrix <- left_join(coord_df, gene_matrix, by = "Geneid")

# --- Final sort and export ---
gene_matrix <- gene_matrix %>%
  select(Geneid, Chr, Start, End, Strand, Length, everything()) %>%
  arrange(Geneid)

tx_matrix <- tx_matrix %>%
  select(transcript_id, everything()) %>%
  arrange(transcript_id)

write_tsv(gene_matrix, opt$`gene-output`)
write_tsv(tx_matrix, opt$`transcript-output`)

cat(sprintf("[INFO] Written %d genes to %s\n", nrow(gene_matrix), opt$`gene-output`))
cat(sprintf("[INFO] Written %d transcripts to %s\n", nrow(tx_matrix), opt$`transcript-output`))
