#!/usr/bin/env Rscript

# Render lightweight post-DESeq2 plots from an RNA-seq differential report plan.
# This intentionally uses base R only: the heavy statistical work has already
# happened upstream, and this layer should be cheap to run on login/local nodes.

parse_args <- function(argv) {
  out <- list(
    top_n = "50",
    padj = "0.1",
    log2fc = "1.0",
    transcript_plot_groups = "all,known_compatible,novel_isoform,novel_locus,ambiguous,artifact"
  )
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

parse_list_arg <- function(value) {
  if (is.null(value) || is.na(value) || value == "") {
    return(character())
  }
  unique(trimws(unlist(strsplit(value, "[,[:space:]]+"))))
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
  blank_panel(title, message)
  grDevices::dev.off()
}

blank_panel <- function(title, message) {
  plot.new()
  title(main = title)
  text(0.5, 0.5, message)
}

feature_id_column <- function(table) {
  stat_columns <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  candidates <- setdiff(colnames(table), stat_columns)
  if (length(candidates)) {
    return(candidates[[1]])
  }
  colnames(table)[[1]]
}

numeric_column <- function(table, column) {
  suppressWarnings(as.numeric(table[[column]]))
}

group_label <- function(results, group_name) {
  if (group_name == "all") {
    return("All features")
  }
  if ("transcript_plot_label" %in% colnames(results) && nrow(results) > 0) {
    labels <- unique(results$transcript_plot_label[results$transcript_plot_label != ""])
    if (length(labels)) {
      return(labels[[1]])
    }
  }
  labels <- c(
    known_compatible = "Known/reference-compatible",
    novel_isoform = "Novel isoform",
    novel_locus = "Novel locus",
    ambiguous = "Ambiguous overlap",
    artifact = "Artifact/repeat"
  )
  if (group_name %in% names(labels)) {
    return(labels[[group_name]])
  }
  group_name
}

plot_groups_for_results <- function(results, requested_groups) {
  feature_column <- feature_id_column(results)
  groups <- list(list(
    name = "all",
    label = "All features",
    results = results,
    feature_ids = as.character(results[[feature_column]])
  ))
  if (!("transcript_plot_group" %in% colnames(results))) {
    return(groups)
  }
  requested <- requested_groups[requested_groups != ""]
  if (!length(requested)) {
    requested <- unique(results$transcript_plot_group[results$transcript_plot_group != ""])
  }
  requested <- unique(requested[requested != "all"])
  for (group_name in requested) {
    subset <- results[results$transcript_plot_group == group_name, , drop = FALSE]
    if (nrow(subset) == 0) {
      next
    }
    groups[[length(groups) + 1]] <- list(
      name = group_name,
      label = group_label(subset, group_name),
      results = subset,
      feature_ids = as.character(subset[[feature_column]])
    )
  }
  groups
}

plot_volcano_panel <- function(results, title, padj_cutoff, log2fc_cutoff) {
  ensure_columns(results, c("log2FoldChange", "padj"), "DESeq2 results")
  log2fc <- numeric_column(results, "log2FoldChange")
  padj <- numeric_column(results, "padj")
  keep <- !is.na(log2fc) & !is.na(padj)
  if (!any(keep)) {
    blank_panel(title, "No valid padj/log2FoldChange rows")
    return()
  }
  log2fc <- log2fc[keep]
  y <- -log10(pmax(padj[keep], .Machine$double.xmin))
  significant <- padj[keep] < padj_cutoff & abs(log2fc) >= log2fc_cutoff
  kept_results <- results[keep, , drop = FALSE]

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
  if (any(significant)) {
    label_column <- feature_id_column(kept_results)
    label_candidates <- which(significant)
    label_candidates <- label_candidates[order(padj[keep][label_candidates], decreasing = FALSE)]
    label_candidates <- head(label_candidates, 8)
    text(
      log2fc[label_candidates],
      y[label_candidates],
      labels = kept_results[[label_column]][label_candidates],
      pos = 3,
      cex = 0.55,
      col = "#111111"
    )
  }
}

plot_volcano <- function(groups, path, title, padj_cutoff, log2fc_cutoff) {
  ensure_parent(path)
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (group in groups) {
    plot_volcano_panel(group$results, paste(title, "-", group$label), padj_cutoff, log2fc_cutoff)
  }
}

plot_ma <- function(results, path, title, padj_cutoff, log2fc_cutoff) {
  ensure_columns(results, c("baseMean", "log2FoldChange", "padj"), "DESeq2 results")
  base_mean <- numeric_column(results, "baseMean")
  log2fc <- numeric_column(results, "log2FoldChange")
  padj <- numeric_column(results, "padj")
  keep <- !is.na(base_mean) & !is.na(log2fc) & !is.na(padj)
  if (!any(keep)) {
    blank_pdf(path, title, "No valid baseMean/log2FoldChange/padj rows")
    return()
  }
  x <- log10(base_mean[keep] + 1)
  y <- log2fc[keep]
  significant <- padj[keep] < padj_cutoff & abs(y) >= log2fc_cutoff

  ensure_parent(path)
  grDevices::pdf(path)
  plot(
    x,
    y,
    pch = 16,
    col = ifelse(significant, "#2166ac", "#4d4d4d"),
    xlab = "log10 mean normalized count + 1",
    ylab = "log2 fold change",
    main = title
  )
  abline(h = c(-log2fc_cutoff, 0, log2fc_cutoff), col = c("#737373", "#bdbdbd", "#737373"), lty = c(2, 1, 2))
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

read_coldata <- function(path, sample_ids) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    return(NULL)
  }
  coldata <- read_tsv(path)
  id_column <- intersect(c("library_id", "sample_id", "sample", "id"), colnames(coldata))
  if (!length(id_column)) {
    return(NULL)
  }
  id_column <- id_column[[1]]
  matched <- coldata[match(sample_ids, coldata[[id_column]]), , drop = FALSE]
  if (any(is.na(matched[[id_column]]))) {
    return(NULL)
  }
  rownames(matched) <- matched[[id_column]]
  matched
}

plot_pca <- function(transformed, path, title, coldata = NULL) {
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
  groups <- NULL
  if (!is.null(coldata) && "condition" %in% colnames(coldata)) {
    groups <- as.factor(coldata[rownames(pca$x), "condition"])
  }
  colors <- rep("#1f78b4", nrow(pca$x))
  if (!is.null(groups)) {
    palette <- grDevices::hcl.colors(length(levels(groups)), "Dark 3")
    colors <- palette[as.integer(groups)]
  }
  ensure_parent(path)
  grDevices::pdf(path)
  plot(
    pca$x[, 1],
    pca$x[, 2],
    pch = 16,
    col = colors,
    xlab = paste0("PC1 (", variance[[1]], "%)"),
    ylab = paste0("PC2 (", variance[[2]], "%)"),
    main = title
  )
  text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.7)
  if (!is.null(groups)) {
    legend("topright", legend = levels(groups), col = palette, pch = 16, bty = "n", cex = 0.8)
  }
  grDevices::dev.off()
}

plot_heatmap_panel <- function(transformed, title, top_n, coldata = NULL) {
  if (ncol(transformed) < 2 || nrow(transformed) < 2) {
    blank_panel(title, "Heatmap requires at least two samples and two features")
    return()
  }
  variances <- apply(transformed, 1, stats::var)
  keep <- head(order(variances, decreasing = TRUE), min(top_n, nrow(transformed)))
  matrix <- transformed[keep, , drop = FALSE]
  if (!is.null(coldata) && "condition" %in% colnames(coldata)) {
    order_index <- order(coldata[colnames(matrix), "condition"], colnames(matrix))
    matrix <- matrix[, order_index, drop = FALSE]
  }
  stats::heatmap(matrix, scale = "row", margins = c(8, 8), main = title)
}

plot_heatmap <- function(transformed, groups, path, title, top_n, coldata = NULL) {
  ensure_parent(path)
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (group in groups) {
    group_matrix <- transformed
    if (group$name != "all") {
      keep <- intersect(group$feature_ids, rownames(transformed))
      group_matrix <- transformed[keep, , drop = FALSE]
    }
    plot_heatmap_panel(group_matrix, paste(title, "-", group$label), top_n, coldata)
  }
}

manifest_columns <- c(
  "project",
  "level",
  "contrast_id",
  "status",
  "reason",
  "volcano_pdf",
  "ma_pdf",
  "pca_pdf",
  "heatmap_pdf",
  "vst_tsv",
  "n_features",
  "n_significant"
)

render_row <- function(row, top_n, padj_cutoff, log2fc_cutoff, transcript_plot_groups) {
  if (!identical(row[["status"]], "ready")) {
    return(data.frame(
      project = row[["project"]],
      level = row[["level"]],
      contrast_id = row[["contrast_id"]],
      status = "blocked",
      reason = row[["reason"]],
      volcano_pdf = row[["volcano_pdf"]],
      ma_pdf = row[["ma_pdf"]],
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
  coldata <- read_coldata(row[["coldata"]], colnames(transformed))
  plot_groups <- plot_groups_for_results(results, transcript_plot_groups)
  plot_volcano(plot_groups, row[["volcano_pdf"]], paste(title, "volcano"), padj_cutoff, log2fc_cutoff)
  plot_ma(results, row[["ma_pdf"]], paste(title, "MA"), padj_cutoff, log2fc_cutoff)
  plot_pca(transformed, row[["pca_pdf"]], paste(title, "PCA"), coldata)
  plot_heatmap(transformed, plot_groups, row[["heatmap_pdf"]], paste(title, "heatmap"), top_n, coldata)

  data.frame(
    project = row[["project"]],
    level = row[["level"]],
    contrast_id = row[["contrast_id"]],
    status = "ok",
    reason = "",
    volcano_pdf = row[["volcano_pdf"]],
    ma_pdf = row[["ma_pdf"]],
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
  transcript_plot_groups <- parse_list_arg(args[["transcript_plot_groups"]])

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
      "ma_pdf",
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
      render_row(row, top_n, padj_cutoff, log2fc_cutoff, transcript_plot_groups),
      error = function(err) {
        data.frame(
          project = row[["project"]],
          level = row[["level"]],
          contrast_id = row[["contrast_id"]],
          status = "failed",
          reason = conditionMessage(err),
          volcano_pdf = row[["volcano_pdf"]],
          ma_pdf = row[["ma_pdf"]],
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
