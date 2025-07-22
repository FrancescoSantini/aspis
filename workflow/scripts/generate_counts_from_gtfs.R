#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})

# ---- CLI parsing ----
option_list <- list(
  make_option(c("--gtfs"), type = "character", help = "Comma-separated list of GTF files"),
  make_option(c("--phenodata"), type = "character", help = "Phenodata CSV with sample names"),
  make_option(c("--gene-output"), type = "character", help = "Path to write gene count matrix"),
  make_option(c("--transcript-output"), type = "character", help = "Path to write transcript count matrix")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Load inputs ----
gtf_paths <- strsplit(opt$gtfs, ",")[[1]]
phenodata <- read_csv(opt$phenodata, show_col_types = FALSE)
sample_names <- rownames(read.csv(opt$phenodata, row.names = 1))

# ---- Helper regex ----
extract <- function(pattern, text) {
  res <- str_match(text, pattern)[,2]
  ifelse(is.na(res), "", res)
}

# Containers
gene_data <- list()
tx_data <- list()
gene_coords <- list()

# ---- Process each GTF ----
for (gtf_file in gtf_paths) {
  sample_id <- sub(".*/([^/]+)\\.gtf$", "\\1", gtf_file)
  lines <- readLines(gtf_file)
  tx_lines <- lines[grepl("\ttranscript\t", lines)]

  gene_counts <- list()
  tx_counts <- list()

  for (line in tx_lines) {
    fields <- strsplit(line, "\t")[[1]]
    attr <- fields[9]
    chr <- fields[1]
    start <- as.integer(fields[4])
    end <- as.integer(fields[5])
    strand <- fields[7]
    len <- end - start + 1

    tx_id <- extract('transcript_id "([^"]+)"', attr)
    g_id <- extract('gene_id "([^"]+)"', attr)
    g_name <- extract('gene_name "([^"]+)"', attr)
    coverage <- as.numeric(extract('cov "([^"]+)"', attr))

    if (is.na(coverage)) next
    count <- ceiling(coverage * len / 75)

    gene_key <- if (g_name != "") paste0(g_id, "|", g_name) else g_id

    # Transcript-level
    tx_counts[[tx_id]] <- count

    # Gene-level
    gene_counts[[gene_key]] <- ifelse(!is.null(gene_counts[[gene_key]]), gene_counts[[gene_key]] + count, count)

    # Coordinates per gene
    if (is.null(gene_coords[[gene_key]])) {
      gene_coords[[gene_key]] <- list(chr = chr, start = start, end = end, strand = strand)
    } else {
      gene_coords[[gene_key]]$start <- min(gene_coords[[gene_key]]$start, start)
      gene_coords[[gene_key]]$end <- max(gene_coords[[gene_key]]$end, end)
    }
  }

  gene_data[[sample_id]] <- gene_counts
  tx_data[[sample_id]] <- tx_counts
}

# ---- Build matrices ----
build_matrix <- function(data_list, id_name) {
  df <- bind_rows(lapply(data_list, ~as_tibble_row(.)))
  df[is.na(df)] <- 0
  rownames(df) <- names(data_list)
  df <- df[, sort(colnames(df)), drop = FALSE]
  df <- as.data.frame(t(df))
  df <- rownames_to_column(df, id_name)
  return(df)
}

gene_matrix <- gene_matrix[order(gene_matrix$Geneid), ]
tx_matrix <- tx_matrix[order(tx_matrix$transcript_id), ]

# ---- Add coordinates to gene matrix ----
coord_df <- tibble(
  Geneid = names(gene_coords),
  Chr = sapply(gene_coords, `[[`, "chr"),
  Start = sapply(gene_coords, `[[`, "start"),
  End = sapply(gene_coords, `[[`, "end"),
  Strand = sapply(gene_coords, `[[`, "strand")
) %>%
  mutate(Length = End - Start + 1)

gene_matrix <- left_join(coord_df, gene_matrix, by = "Geneid")

# ---- Write outputs ----
write_csv(gene_matrix, opt$`gene-output`)
write_csv(tx_matrix, opt$`transcript-output`)

cat(sprintf("[INFO] Written %d genes to %s\n", nrow(gene_matrix), opt$`gene-output`))
cat(sprintf("[INFO] Written %d transcripts to %s\n", nrow(tx_matrix), opt$`transcript-output`))
