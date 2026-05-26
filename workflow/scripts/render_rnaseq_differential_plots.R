#!/usr/bin/env Rscript

# Render lightweight post-DESeq2 plots from an RNA-seq differential report plan.
# This intentionally uses base R only: the heavy statistical work has already
# happened upstream, and this layer should be cheap to run on login/local nodes.

parse_args <- function(argv) {
  out <- list(top_n = "50", padj = "0.1", log2fc = "1.0")
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (i == length(argv)) {
      stop("Missing value for argument: ", key)
    }
    name <- chartr("-", "_", substring(key, 3))
    out[[name]] <- argv[[i + 1]]
    i <- i + 2
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || value == "") {
    stop("Missing required argument: --", name)
  }
  value
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Input does not exist: ", path)
  }
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

ensure_columns <- function(table, columns, label) {
  missing <- setdiff(columns, colnames(table))
  if (length(missing)) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

ensure_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

blank_pdf <- function(path, title, message) {
  ensure_parent(path)
  grDevices::pdf(path)
  plot.new()
  title(main = title)
  text(0.5, 0.5, message)
  grDevices::dev.off()
}

numeric_column <- function(table, column) {
  suppressWarnings(as.numeric(table[[column]]))
}

plot_volcano <- function(results, path, title, padj_cutoff, log2fc_cutoff) {
  ensure_columns(results, c("log2FoldChange", "padj"), "DESeq2 results")
  log2fc <- numeric_column(results, "log2FoldChange")
  padj <- numeric_column(results, "padj")
  keep <- !is.na(log2fc) & !is.na(padj)
  if (!any(keep)) {
    blank_pdf(path, title, "No valid padj/log2FoldChange rows")
    return()
  }
  log2fc <- log2fc[keep]
  y <- -log10(pmax(padj[keep], .Machine$double.xmin))
  significant <- padj[keep] < padj_cutoff & abs(log2fc) >= log2fc_cutoff

  ensure_parent(path)
  grDevices::pdf(path)
  plot(
    log2fc,
    y,
    pch = 16,
    col = ifelse(significant, "#b2182b", "#4d4d4d"),
    xlab = "log2 fold change",
    ylab = "-log10 adjusted p-value",
    main = title
  )
  abline(v = c(-log2fc_cutoff, log2fc_cutoff), col = "#737373", lty = 2)
  abline(h = -log10(padj_cutoff), col = "#737373", lty = 2)
  grDevices::dev.off()
}

read_normalized_counts <- function(path) {
  counts <- read_tsv(path)
  if (ncol(counts) < 2) {
    stop("Normalized counts table needs one feature column and at least one sample column: ", path)
  }
  feature_id_column <- colnames(counts)[[1]]
  feature_ids <- counts[[feature_id_column]]
  matrix <- as.matrix(counts[, -1, drop = FALSE])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- feature_ids
  matrix
}

write_transformed_counts <- function(matrix, path) {
  transformed <- log2(matrix + 1)
  output <- data.frame(feature_id = rownames(transformed), transformed, check.names = FALSE)
  ensure_parent(path)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
  transformed
}

plot_pca <- function(transformed, path, title) {
  if (ncol(transformed) < 2 || nrow(transformed) < 2) {
    blank_pdf(path, title, "PCA requires at least two samples and two features")
    return()
  }
  pca <- stats::prcomp(t(transformed), center = TRUE, scale. = FALSE)
  if (ncol(pca$x) < 2 || sum(pca$sdev^2) == 0) {
    blank_pdf(path, title, "PCA requires at least two variable dimensions")
    return()
  }
  variance <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  ensure_parent(path)
  grDevices::pdf(path)
  plot(
    pca$x[, 1],
    pca$x[, 2],
    pch = 16,
    xlab = paste0("PC1 (", variance[[1]], "%)"),
    ylab = paste0("PC2 (", variance[[2]], "%)"),
    main = title
  )
  text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.7)
  grDevices::dev.off()
}

plot_heatmap <- function(transformed, path, title, top_n) {
  if (ncol(transformed) < 2 || nrow(transformed) < 2) {
    blank_pdf(path, title, "Heatmap requires at least two samples and two features")
    return()
  }
  variances <- apply(transformed, 1, stats::var)
  keep <- head(order(variances, decreasing = TRUE), min(top_n, nrow(transformed)))
  matrix <- transformed[keep, , drop = FALSE]
  ensure_parent(path)
  grDevices::pdf(path)
  stats::heatmap(matrix, scale = "row", margins = c(8, 8), main = title)
  grDevices::dev.off()
}

manifest_columns <- c(
  "project",
  "level",
  "contrast_id",
  "status",
  "reason",
  "volcano_pdf",
  "pca_pdf",
  "heatmap_pdf",
  "vst_tsv",
  "n_features",
  "n_significant"
)

render_row <- function(row, top_n, padj_cutoff, log2fc_cutoff) {
  if (!identical(row[["status"]], "ready")) {
    return(data.frame(
      project = row[["project"]],
      level = row[["level"]],
      contrast_id = row[["contrast_id"]],
      status = "blocked",
      reason = row[["reason"]],
      volcano_pdf = row[["volcano_pdf"]],
      pca_pdf = row[["pca_pdf"]],
      heatmap_pdf = row[["heatmap_pdf"]],
      vst_tsv = row[["vst_tsv"]],
      n_features = 0,
      n_significant = 0,
      check.names = FALSE
    ))
  }

  title <- paste(row[["project"]], row[["level"]], row[["contrast_id"]])
  results <- read_tsv(row[["results"]])
  filtered <- read_tsv(row[["filtered"]])
  normalized <- read_normalized_counts(row[["normalized_counts"]])
  transformed <- write_transformed_counts(normalized, row[["vst_tsv"]])
  plot_volcano(results, row[["volcano_pdf"]], paste(title, "volcano"), padj_cutoff, log2fc_cutoff)
  plot_pca(transformed, row[["pca_pdf"]], paste(title, "PCA"))
  plot_heatmap(transformed, row[["heatmap_pdf"]], paste(title, "heatmap"), top_n)

  data.frame(
    project = row[["project"]],
    level = row[["level"]],
    contrast_id = row[["contrast_id"]],
    status = "ok",
    reason = "",
    volcano_pdf = row[["volcano_pdf"]],
    pca_pdf = row[["pca_pdf"]],
    heatmap_pdf = row[["heatmap_pdf"]],
    vst_tsv = row[["vst_tsv"]],
    n_features = nrow(results),
    n_significant = nrow(filtered),
    check.names = FALSE
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  plan_path <- required_arg(args, "plan")
  manifest_path <- required_arg(args, "manifest")
  done_path <- required_arg(args, "done")
  top_n <- as.integer(args[["top_n"]])
  if (is.na(top_n) || top_n < 1) {
    stop("--top-n/--top_n must be a positive integer")
  }
  padj_cutoff <- as.numeric(args[["padj"]])
  if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
    stop("--padj must be a number in (0, 1]")
  }
  log2fc_cutoff <- as.numeric(args[["log2fc"]])
  if (is.na(log2fc_cutoff) || log2fc_cutoff < 0) {
    stop("--log2fc must be a non-negative number")
  }

  plan <- read_tsv(plan_path)
  ensure_columns(
    plan,
    c(
      "project",
      "level",
      "contrast_id",
      "status",
      "reason",
      "results",
      "filtered",
      "normalized_counts",
      "volcano_pdf",
      "pca_pdf",
      "heatmap_pdf",
      "vst_tsv"
    ),
    "Differential report plan"
  )
  if (nrow(plan) == 0) {
    stop("Differential report plan has no rows")
  }

  rows <- list()
  for (i in seq_len(nrow(plan))) {
    row <- as.list(plan[i, , drop = FALSE])
    rows[[i]] <- tryCatch(
      render_row(row, top_n, padj_cutoff, log2fc_cutoff),
      error = function(err) {
        data.frame(
          project = row[["project"]],
          level = row[["level"]],
          contrast_id = row[["contrast_id"]],
          status = "failed",
          reason = conditionMessage(err),
          volcano_pdf = row[["volcano_pdf"]],
          pca_pdf = row[["pca_pdf"]],
          heatmap_pdf = row[["heatmap_pdf"]],
          vst_tsv = row[["vst_tsv"]],
          n_features = 0,
          n_significant = 0,
          check.names = FALSE
        )
      }
    )
  }
  manifest <- do.call(rbind, rows)
  ensure_parent(manifest_path)
  write.table(manifest[, manifest_columns], manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

  ok <- sum(manifest$status == "ok")
  blocked <- sum(manifest$status == "blocked")
  failed <- sum(manifest$status == "failed")
  status <- if (failed > 0) "failed" else if (ok > 0 && blocked == 0) "ok" else "blocked"
  ensure_parent(done_path)
  done <- data.frame(
    status = status,
    plots_ok = ok,
    plots_blocked = blocked,
    plots_failed = failed,
    plots_total = nrow(manifest),
    check.names = FALSE
  )
  write.table(done, done_path, sep = "\t", quote = FALSE, row.names = FALSE)
  if (failed > 0) {
    failed_ids <- paste(manifest$contrast_id[manifest$status == "failed"], collapse = ", ")
    stop("Differential plot rendering failed for contrast(s): ", failed_ids)
  }
}

main()
