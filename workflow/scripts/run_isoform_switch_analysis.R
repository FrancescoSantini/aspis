#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DRIMSeq)
  library(IsoformSwitchAnalyzeR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(readr)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(ggplot2)
})

# --- CLI parsing ---
option_list <- list(
  make_option("--counts", type = "character", help = "Transcript count matrix (CSV)"),
  make_option("--phenodata", type = "character", help = "Phenodata file (CSV)"),
  make_option("--gtf", type = "character", help = "Merged GTF file"),
  make_option("--outdir", type = "character", help = "Output directory"),
  make_option("--gene-expr", type = "double", default = 1, help = "Gene expression cutoff"),
  make_option("--isoform-expr", type = "double", default = 1, help = "Isoform expression cutoff"),
  make_option("--condition", type = "character", help = "Experimental condition to compare vs control"),
  make_option("--covariate", type = "character", help = "Covariate1 value to subset samples"),
  make_option("--padj", type = "double", default = 0.05, help = "Adjusted p-value cutoff for plotting switches"),
  make_option("--dIF", type = "double", default = 0.2, help = "Minimum |dIF| (effect size) for plotting"),
  make_option("--max_genes", type = "integer", default = 20, help = "Maximum number of genes to visualize")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- thresholds used throughout plotting/selection ----
padj_cutoff      <- opt$padj
dIF_cutoff       <- opt$dIF
max_genes_to_plot <- opt$max_genes

cat("[DEBUG] IsoformSwitchAnalyzeR version: ", as.character(packageVersion("IsoformSwitchAnalyzeR")), "\n")

# --- Load counts ---
cat("[INFO] Loading transcript counts...\n")
counts <- read_tsv(opt$counts)
colnames(counts)[1] <- "isoform_id"
colnames(counts) <- gsub(".*/", "", colnames(counts))
colnames(counts) <- gsub("\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("_sorted$", "", colnames(counts))
cat("[DEBUG] colnames(counts):\n")
print(colnames(counts))

# --- Load phenodata ---
cat("[INFO] Loading phenodata...\n")
pheno <- read.csv(opt$phenodata, row.names = 1)
pheno$sampleID <- rownames(pheno)

if (!is.null(opt$covariate)) {
  cat(sprintf("[INFO] Subsetting to covariate1 == %s\n", opt$covariate))
  pheno <- pheno[pheno$covariate1 == opt$covariate, ]
}
if (!is.null(opt$condition)) {
  cat(sprintf("[INFO] Subsetting to condition in (control, %s)\n", opt$condition))
  pheno <- pheno[pheno$condition %in% c("control", opt$condition), ]
  pheno$condition <- factor(pheno$condition, levels = c("control", opt$condition))
}

# Clean design matrix
is_constant <- function(col) length(unique(na.omit(col))) <= 1
to_keep <- !sapply(pheno, is_constant) | names(pheno) %in% c("sampleID", "condition")
pheno <- pheno[, to_keep, drop = FALSE]
pheno <- pheno[, c("sampleID", "condition")]

cat("[INFO] Design matrix after filtering:\n")
print(str(pheno))
cat("[DEBUG] colnames(pheno):\n")
print(colnames(pheno))

# --- Match columns ---
samples <- intersect(colnames(counts)[-1], pheno$sampleID)
counts <- counts[, c("isoform_id", samples)]
pheno <- pheno[pheno$sampleID %in% samples, ]
cat(sprintf("[INFO] Matched %d samples\n", length(samples)))

# --- Construct switchAnalyzeRlist ---
cat("[INFO] Building switchAnalyzeRlist...\n")
switchList <- importRdata(
  isoformCountMatrix = counts,
  designMatrix = pheno,
  isoformExonAnno = opt$gtf,
  comparisonsToMake = NULL,
  ignoreAfterBar = FALSE,
  showProgress = TRUE
)

if (all(is.na(switchList$isoformFeatures$gene_name))) {
  cat("[INFO] gene_name missing in import — recovering from GTF directly\n")

  gtf <- rtracklayer::import(opt$gtf)
  gtf_df <- as.data.frame(gtf)
  gtf_tx <- gtf_df[gtf_df$type == "transcript", ]
  tx2name <- gtf_tx[, c("transcript_id", "gene_name")]
  tx2name <- tx2name[!is.na(tx2name$gene_name), ]

  matched_names <- tx2name$gene_name[match(
    switchList$isoformFeatures$isoform_id,
    tx2name$transcript_id
  )]
  switchList$isoformFeatures$gene_name <- matched_names
}

cat("[DEBUG] Top 5 gene_id values in switchAnalyzeRlist:\n")
print(head(switchList$isoformFeatures$gene_id, 5))

cat("[DEBUG] Top 5 gene_name values:\n")
print(head(switchList$isoformFeatures$gene_name, 5))

# --- Manual expression estimation ---
cat("[INFO] Estimating expression manually...\n")
count_mat <- as.data.frame(counts)
rownames(count_mat) <- count_mat$isoform_id
count_mat <- count_mat[, -1]
lib_sizes <- colSums(count_mat)
cpm <- t(t(count_mat) / lib_sizes * 1e6)
isoform_expr <- rowMeans(cpm)
switchList$isoformFeatures$isoform_expression <- isoform_expr[match(switchList$isoformFeatures$isoform_id, rownames(cpm))]

gene_expr <- switchList$isoformFeatures %>%
  group_by(gene_id) %>%
  summarise(gene_expression = sum(isoform_expression, na.rm = TRUE))
switchList$isoformFeatures <- left_join(switchList$isoformFeatures, gene_expr, by = "gene_id")
cat("[INFO] Expression estimation complete.\n")

cat("[DEBUG] Columns in switchList$isoformFeatures:\n")
print(colnames(switchList$isoformFeatures))

# --- Write expression summary ---
expr_out <- file.path(opt$outdir, "expression_summary.txt")
sink(expr_out)
cat("[DEBUG] Summary of gene_expression:\n")
print(summary(switchList$isoformFeatures$gene_expression))
cat("\n[DEBUG] Summary of isoform_expression:\n")
print(summary(switchList$isoformFeatures$isoform_expression))
sink()

# --- Pre-filter ---
switchListFiltered <- preFilter(switchList, removeSingleIsoformGenes = TRUE)

# --- Manual filtering ---
cat("[INFO] Applying expression filters: gene >= ", opt$`gene-expr`, ", isoform >= ", opt$`isoform-expr`, "\n")
expr_pass <- with(switchListFiltered$isoformFeatures,
  gene_expression >= opt$`gene-expr` & isoform_expression >= opt$`isoform-expr`
)

if (sum(expr_pass) == 0) {
  cat("[WARN] No isoforms pass filters — writing empty result.\n")
  write.csv(data.frame(), file.path(opt$outdir, "isoform_switch_summary.csv"), row.names = FALSE)
  file.create(file.path(opt$outdir, "isoform_switch_qc.pdf"))
  file.create(file.path(opt$outdir, "switch_list.rds"))
  quit(save = "no", status = 0)
}

switchListFiltered <- subsetSwitchAnalyzeRlist(switchListFiltered, expr_pass)

# --- ORF prediction ---
cat("[INFO] Running ORF prediction (longest)...\n")
switchListFiltered <- analyzeORF(switchListFiltered, orfMethod = "longest", genomeObject = BSgenome.Hsapiens.UCSC.hg38)
cat("[INFO] ORF prediction complete.\n")

# --- Isoform switch testing ---
cat("[INFO] Running DEXSeq isoform switch test...\n")
switchListTested <- isoformSwitchTestDEXSeq(switchListFiltered)

# --- DEBUG: Check for dIF column ---
cat("[DEBUG] Checking isoformSwitchAnalysis for dIF column...\n")
print(colnames(switchListTested$isoformSwitchAnalysis))

if (!"dIF" %in% colnames(switchListTested$isoformSwitchAnalysis)) {
  cat("[ERROR] dIF column is missing in isoformSwitchAnalysis!\n")
} else {
  cat("[INFO] dIF column is present. Summary:\n")
  print(summary(switchListTested$isoformSwitchAnalysis$dIF))
}

# --- Extract ORF sequences (AA) ---
cat("[INFO] Extracting amino acid sequences...\n")
switchListTested <- extractSequence(switchListTested, genomeObject = BSgenome.Hsapiens.UCSC.hg38)
cat("[INFO] Amino acid sequences extracted.\n")

# --- Intron retention classification ---
cat("[INFO] Classifying intron retention...\n")
switchListTested <- analyzeIntronRetention(switchListTested)
cat("[INFO] Intron retention classification complete.\n")

# --- Consequence analysis ---
cat("[INFO] Analyzing consequences of isoform switches...\n")
switchListTested <- analyzeSwitchConsequences(
  switchListTested,
  consequencesToAnalyze = c("intron_retention", "NMD_status", "ORF_seq_similarity")
)
cat("[INFO] Consequence analysis complete.\n")

# --- Output comparison summary (unchanged) ---
summary_path <- file.path(opt$outdir, "isoform_switch_summary.csv")
summary_df <- extractSwitchSummary(switchListTested)
write.csv(summary_df, summary_path, row.names = FALSE)

# --- Build 'joined' isoform-level table (analysis + gene names) ---
joined <- left_join(
  switchListTested$isoformSwitchAnalysis,
  switchListTested$isoformFeatures[, c("isoform_id", "gene_id", "gene_name")],
  by = "isoform_id"
)

# --- Long-format consequences (isoform, direction, consequence), then COLLAPSE ---
conseq_df <- switchListTested$switchConsequence
if (!is.data.frame(conseq_df) ||
    !all(c("isoformUpregulated","isoformDownregulated","switchConsequence") %in% names(conseq_df))) {
  cat("[WARN] switchConsequence lacked expected columns; writing NA-only consequences.\n")
  collapsed_conseq <- tibble(isoform_id = character(), direction = character(), consequence = character())
} else {
  up_conseq <- conseq_df %>%
    filter(!is.na(isoformUpregulated)) %>%
    transmute(isoform_id = isoformUpregulated,
              direction = "up",
              consequence = switchConsequence)

  down_conseq <- conseq_df %>%
    filter(!is.na(isoformDownregulated)) %>%
    transmute(isoform_id = isoformDownregulated,
              direction = "down",
              consequence = switchConsequence)

  conseq_long <- bind_rows(up_conseq, down_conseq)

  # collapse: per (isoform_id, direction), keep unique non-NA consequences; drop lone NA duplicates
  collapsed_conseq <- conseq_long %>%
    mutate(consequence = ifelse(is.na(consequence) | consequence == "", NA_character_, consequence)) %>%
    group_by(isoform_id, direction) %>%
    summarise(
      consequence = {
        uniq <- sort(unique(na.omit(consequence)))
        if (length(uniq) == 0) NA_character_ else paste(uniq, collapse = "; ")
      },
      .groups = "drop"
    )
}

# --- Write out detailed and COLLAPSED CSVs for downstream sanity checks ---
detail_path      <- file.path(opt$outdir, "isoform_switch_detailed.csv")
detail_collapsed <- file.path(opt$outdir, "isoform_switch_detailed_collapsed.csv")

# detailed with one row per (isoform_id, direction) after collapsing consequences
joined_collapsed <- joined %>%
  # Attach both 'up' and 'down' (may create up to 2 rows per isoform)
  left_join(collapsed_conseq, by = "isoform_id")

write.csv(joined,            detail_path,      row.names = FALSE)
write.csv(joined_collapsed,  detail_collapsed, row.names = FALSE)
cat("[INFO] Wrote detailed tables:\n  - ", detail_path, "\n  - ", detail_collapsed, "\n", sep = "")

cat("[DEBUG] Columns in extractSwitchSummary:\n")
print(colnames(summary_df))

# --- Output filtered isoform table ---
filtered_table_path <- file.path(opt$outdir, "filtered_isoforms.csv")
write.csv(switchListTested$isoformSwitchAnalysis, filtered_table_path, row.names = FALSE)

# --- Output dIF distribution plot ---
cat("[INFO] Plotting dIF distribution...\n")
difs <- switchListTested$isoformSwitchAnalysis
dif_plot_path <- file.path(opt$outdir, "dif_distribution.pdf")

if (!"dIF" %in% colnames(difs)) {
  cat("[ERROR] dIF column is missing in isoformSwitchAnalysis — skipping plot.\n")
  file.create(dif_plot_path)
} else {
  cat("[DEBUG] Summary of dIF column:\n")
  print(summary(difs$dIF))
  cat("[DEBUG] dIF quantiles:\n")
  print(quantile(difs$dIF, probs = seq(0, 1, 0.1)))
  cat("[DEBUG] First few rows with dIF:\n")
  print(head(difs[, c("isoform_id", "dIF")]))

  pdf(dif_plot_path)
  print(
    ggplot(difs, aes(x = dIF)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black", size = 0.25) +
      theme_minimal() +
      labs(
        title = "Distribution of dIF values",
        x = "dIF (Delta Isoform Fraction)",
        y = "Count"
      )
  )
  dev.off()

  cat("[INFO] dIF distribution plot saved.\n")
}

# --- Plot top 10 isoform switches ---
cat("[TRACE] === Generating isoform_switch_qc.pdf ===\n")
qc_plot_path <- file.path(opt$outdir, "isoform_switch_qc.pdf")

# Join for gene names
joined <- left_join(
  switchListTested$isoformSwitchAnalysis,
  switchListTested$isoformFeatures[, c("isoform_id", "gene_id", "gene_name")],
  by = "isoform_id"
)

# Write full table for logging
write.csv(joined, file.path(opt$outdir, "isoform_switch_detailed.csv"), row.names = FALSE)

# Select top genes with significant isoform switching
padj_cutoff <- opt$padj
dIF_cutoff <- opt$dIF
max_genes_to_plot <- opt$max_genes

# --- Apply filtering ---
filtered <- joined %>%
  filter(!is.na(padj), padj < padj_cutoff, abs(dIF) >= dIF_cutoff) %>%
  arrange(desc(abs(dIF))) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  slice_head(n = max_genes_to_plot)

top_genes <- filtered$gene_id

cat(sprintf("[INFO] Selected %d genes with padj < %.3f and |dIF| >= %.2f\n",
            length(top_genes), padj_cutoff, dIF_cutoff))

if (length(top_genes) == 0) {
  pdf(qc_plot_path)
  plot.new()
  title("No genes to plot")
  dev.off()
  cat("[TRACE] No top genes found — empty plot written.\n")
} else {
  plots_rendered <- 0
  pdf(qc_plot_path, width = 8, height = 10)

  # --- Page 1: Title and dIF distribution ---
  par(mfrow = c(2, 1), mar = c(2, 4, 5, 1))
  plot.new()
  title("Top Isoform Switches by dIF", cex.main = 1.8, line = -1)
  mtext(paste("Comparison:", opt$condition, "vs control"), side = 3, line = 2, cex = 1.2)
  mtext(paste("Samples:", length(unique(pheno$sampleID)), 
              "- Genes analyzed:", length(unique(switchListTested$isoformFeatures$gene_id))), 
        side = 3, line = 3, cex = 1)
  mtext(paste("Cutoffs - Gene expr:", opt$`gene-expr`, ", Isoform expr:", opt$`isoform-expr`), 
        side = 3, line = 4, cex = 1)
  mtext(paste("Switch thresholds - padj <", opt$padj, 
            ", |dIF| >=", opt$dIF, 
            ", max genes:", opt$max_genes), 
      side = 3, line = 5, cex = 1)

  hist(switchListTested$isoformSwitchAnalysis$dIF,
       main = "Distribution of dIF Values",
       xlab = "dIF (Delta Isoform Fraction)",
       col = "steelblue", border = "black", breaks = 30)

  par(mfrow = c(1, 1))  # Reset layout

  # --- Per-gene plots and tables ---
  for (g in top_genes) {
    cat("[TRACE] Rendering gene: ", g, "\n")

    p <- tryCatch({
      suppressMessages(switchPlot(switchListTested, gene = g))
    }, error = function(e) {
      cat("[ERROR] Failed to plot gene: ", g, ": ", e$message, "\n")
      NULL
    })

    if (!is.null(p)) {
      print(p)

      # isoforms for this gene
      iso_ids <- switchListTested$isoformFeatures %>%
        filter(gene_id == g) %>%
        pull(isoform_id)

      # base table: one row per isoform in the gene
      tab <- switchListTested$isoformSwitchAnalysis %>%
        filter(isoform_id %in% iso_ids) %>%
        select(isoform_id, dIF, padj) %>%
        mutate(
          dIF  = round(dIF, 4),
          padj = ifelse(
            is.na(padj), NA_character_,
            ifelse(padj < 1e-4,
                   formatC(padj, format = "e", digits = 2),
                   formatC(padj, format = "f", digits = 4))
          )
        )

      # attach consequences (collapse duplicates)
      conseq_df <- switchListTested$switchConsequence
      cat("[DEBUG] switchConsequence is a", class(conseq_df), "with", nrow(conseq_df), "rows\n")
      cat("[DEBUG] switchConsequence column names:\n")
      print(colnames(conseq_df))

      if (is.data.frame(conseq_df) &&
          all(c("isoformUpregulated", "isoformDownregulated", "switchConsequence") %in% colnames(conseq_df))) {
        cat("[INFO] Adding isoform consequences from both upregulated and downregulated isoforms\n")

        up_conseq <- conseq_df %>%
          filter(!is.na(isoformUpregulated)) %>%
          transmute(isoform_id = isoformUpregulated,
                    direction  = "up",
                    consequence = switchConsequence)

        down_conseq <- conseq_df %>%
          filter(!is.na(isoformDownregulated)) %>%
          transmute(isoform_id = isoformDownregulated,
                    direction  = "down",
                    consequence = switchConsequence)

        # combine and collapse duplicates:
        # - keep unique non-NA consequences
        # - if only NA exists for (isoform_id, direction), keep a single NA
        conseq_collapsed <- bind_rows(up_conseq, down_conseq) %>%
          group_by(isoform_id, direction) %>%
          summarise(
            consequence = {
              vals <- unique(consequence)
              non_na <- vals[!is.na(vals) & vals != ""]
              if (length(non_na) == 0) NA_character_
              else paste(sort(non_na), collapse = "; ")
            },
            .groups = "drop"
          )

        tab <- left_join(tab, conseq_collapsed, by = "isoform_id") %>%
          # nicer ordering: down first, then up; strongest |dIF| first
          mutate(direction = factor(direction, levels = c("down", "up"))) %>%
          arrange(direction, desc(abs(as.numeric(dIF))))
      } else {
        cat("[INFO] Skipping consequence annotation: required isoform-level columns missing.\n")
      }

      # draw table
      table_grob <- gridExtra::tableGrob(tab)
      grid::grid.draw(table_grob)
      plots_rendered <- plots_rendered + 1
    } else {
      cat("[TRACE] Skipping gene ", g, " — plot is NULL.\n")
    }
  }

  dev.off()
  cat("[TRACE] Closed isoform_switch_qc.pdf with", plots_rendered, "plots.\n")
}

# --- Save switch consequences summary ---
conseq_out <- file.path(opt$outdir, "switch_consequences_summary.csv")
write.csv(switchListTested$switchConsequence, conseq_out, row.names = FALSE)
cat("[INFO] Consequence summary saved.\n")

# --- Save FASTA outputs (redirected and moved manually) ---
cat("[INFO] Extracting sequences with output redirection...\n")

extractSequence(
  switchListTested,
  genomeObject = BSgenome.Hsapiens.UCSC.hg38,
  outputPrefix = file.path(opt$outdir, "isoformSwitchAnalyzeR")
)

# Manually move the internal hardcoded FASTAs from root to outdir
internal_nt <- "isoformSwitchAnalyzeR_isoform_nt.fasta"
internal_aa <- "isoformSwitchAnalyzeR_isoform_AA.fasta"
out_nt <- file.path(opt$outdir, internal_nt)
out_aa <- file.path(opt$outdir, internal_aa)

if (file.exists(internal_nt)) {
  file.rename(internal_nt, out_nt)
  cat("[INFO] Moved nt FASTA to: ", out_nt, "\n")
} else {
  cat("[WARN] nt FASTA not found in root.\n")
}

if (file.exists(internal_aa)) {
  file.rename(internal_aa, out_aa)
  cat("[INFO] Moved AA FASTA to: ", out_aa, "\n")
} else {
  cat("[WARN] AA FASTA not found in root.\n")
}

cat("[INFO] All FASTA sequences saved and moved to outdir.\n")

# --- Save full R object ---
saveRDS(switchListTested, file = file.path(opt$outdir, "switch_list.rds"))
cat("[INFO] Analysis completed and results saved.\n")
