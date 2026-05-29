#!/usr/bin/env Rscript

# Render lightweight post-DESeq2 plots from an RNA-seq differential report plan.
# This intentionally uses base R only: the heavy statistical work has already
# happened upstream, and this layer should be cheap to run on login/local nodes.

parse_args <- function(argv) {
  out <- list(
    top_n = "50",
    padj = "0.1",
    log2fc = "1.0",
    pca_color_columns = "condition,time,time_h,batch,batch_id,biospecimen,biospecimen_id,replicate,replicate_id",
    transcript_plot_groups = "all,known_compatible,novel_isoform,novel_locus,ambiguous,artifact",
    gene_biotype_plot_groups = "protein_coding,lncRNA,pseudogene,snoRNA,snRNA,miRNA",
    transcript_biotype_plot_groups = "protein_coding,lncRNA,pseudogene,snoRNA,snRNA,miRNA",
    mirna_plot_groups = "all,up,down,arm,target_source,target_source_type,target_evidence_type",
    heatmap_modes = "significant,variable",
    heatmap_feature_lists = "",
    heatmap_significant_fallback = "variable"
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

safe_group_name <- function(value) {
  cleaned <- gsub("[^A-Za-z0-9]+", "_", tolower(trimws(value)))
  cleaned <- gsub("^_+|_+$", "", cleaned)
  ifelse(cleaned == "", "unclassified", cleaned)
}

biotype_label <- function(value) {
  if (is.null(value) || is.na(value) || value == "") {
    return("unclassified")
  }
  value
}

first_existing_column <- function(table, columns) {
  matched <- columns[columns %in% colnames(table)]
  if (length(matched)) {
    return(matched[[1]])
  }
  ""
}

biotype_column_for_level <- function(results, level) {
  if (level == "transcript") {
    return(first_existing_column(
      results,
      c("transcript_biotype", "transcript_type", "gene_biotype", "gene_type", "biotype", "feature_type")
    ))
  }
  first_existing_column(results, c("gene_biotype", "gene_type", "biotype", "feature_type"))
}

requested_biotypes_present <- function(values, requested_biotypes) {
  clean_values <- values[!is.na(values) & values != ""]
  if (!length(clean_values)) {
    return(character())
  }
  if (!length(requested_biotypes)) {
    return(unique(clean_values))
  }
  requested_keys <- safe_group_name(requested_biotypes)
  value_keys <- safe_group_name(clean_values)
  selected <- vapply(seq_along(clean_values), function(index) {
    value_key <- value_keys[[index]]
    any(value_key %in% requested_keys | grepl(paste(requested_keys, collapse = "|"), value_key))
  }, logical(1))
  unique(clean_values[selected])
}

append_plot_group <- function(groups, name, label, type, subset, feature_column) {
  if (nrow(subset) == 0 || name %in% vapply(groups, function(group) group$name, character(1))) {
    return(groups)
  }
  groups[[length(groups) + 1]] <- list(
    name = name,
    label = label,
    type = type,
    results = subset,
    feature_ids = as.character(subset[[feature_column]])
  )
  groups
}

mirna_arm <- function(feature_ids) {
  ids <- tolower(as.character(feature_ids))
  ifelse(
    grepl("(^|[-_])5p$", ids),
    "5p",
    ifelse(grepl("(^|[-_])3p$", ids), "3p", "unannotated_arm")
  )
}

split_group_values <- function(values) {
  values <- unique(trimws(unlist(strsplit(as.character(values), "[,;|]"))))
  values[!is.na(values) & values != ""]
}

read_mirna_target_rows <- function(path) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) {
    return(data.frame())
  }
  read_tsv(path)
}

mirna_target_feature_ids <- function(targets) {
  if (nrow(targets) == 0) {
    return(character())
  }
  id_column <- first_existing_column(
    targets,
    c("mirna_id", "mature_mirna_id", "miRNA", "mirna", "mature_id", "Geneid", "feature_id")
  )
  if (id_column == "") {
    return(character())
  }
  as.character(targets[[id_column]])
}

add_mirna_target_groups <- function(groups, results, feature_column, target_rows, group_column, group_type, label_prefix) {
  if (!nrow(target_rows) || !(group_column %in% colnames(target_rows))) {
    return(groups)
  }
  result_ids <- as.character(results[[feature_column]])
  for (value in sort(split_group_values(target_rows[[group_column]]))) {
    members <- unique(mirna_target_feature_ids(target_rows[target_rows[[group_column]] == value, , drop = FALSE]))
    matched_ids <- intersect(members, result_ids)
    if (!length(matched_ids)) {
      next
    }
    subset <- results[result_ids %in% matched_ids, , drop = FALSE]
    groups <- append_plot_group(
      groups,
      paste(group_type, safe_group_name(value), sep = "__"),
      paste(label_prefix, value, sep = ": "),
      group_type,
      subset,
      feature_column
    )
  }
  groups
}

add_mirna_plot_groups <- function(groups, results, requested_groups, padj_cutoff, log2fc_cutoff, plan_row) {
  requested <- requested_groups[requested_groups != ""]
  if (!length(requested)) {
    requested <- c("all", "up", "down", "arm", "target_source", "target_source_type", "target_evidence_type")
  }
  feature_column <- feature_id_column(results)
  log2fc <- if ("log2FoldChange" %in% colnames(results)) numeric_column(results, "log2FoldChange") else rep(NA_real_, nrow(results))
  padj <- if ("padj" %in% colnames(results)) numeric_column(results, "padj") else rep(NA_real_, nrow(results))
  significant <- !is.na(padj) & padj < padj_cutoff & !is.na(log2fc) & abs(log2fc) >= log2fc_cutoff

  if ("up" %in% requested || "direction" %in% requested || "regulated" %in% requested) {
    groups <- append_plot_group(
      groups,
      "mirna_up",
      "Upregulated miRNAs",
      "mirna_direction",
      results[significant & log2fc > 0, , drop = FALSE],
      feature_column
    )
  }
  if ("down" %in% requested || "direction" %in% requested || "regulated" %in% requested) {
    groups <- append_plot_group(
      groups,
      "mirna_down",
      "Downregulated miRNAs",
      "mirna_direction",
      results[significant & log2fc < 0, , drop = FALSE],
      feature_column
    )
  }

  if ("arm" %in% requested || "mature_arm" %in% requested) {
    arm <- if ("arm" %in% colnames(results)) as.character(results[["arm"]]) else mirna_arm(results[[feature_column]])
    for (value in c("5p", "3p")) {
      subset <- results[arm == value, , drop = FALSE]
      groups <- append_plot_group(
        groups,
        paste("mirna_arm", value, sep = "__"),
        paste("Mature arm", value, sep = ": "),
        "mirna_arm",
        subset,
        feature_column
      )
    }
  }

  target_rows <- read_mirna_target_rows(plan_row[["mirna_targets"]])
  if ("target_source" %in% requested) {
    groups <- add_mirna_target_groups(
      groups, results, feature_column, target_rows, "target_source", "mirna_target_source", "Target source"
    )
  }
  if ("target_source_type" %in% requested) {
    groups <- add_mirna_target_groups(
      groups, results, feature_column, target_rows, "target_source_type", "mirna_target_source_type", "Target source type"
    )
  }
  if ("target_evidence_type" %in% requested) {
    groups <- add_mirna_target_groups(
      groups, results, feature_column, target_rows, "target_evidence_type", "mirna_target_evidence_type", "Target evidence"
    )
  }

  groups
}

plot_groups_for_results <- function(
  results,
  level,
  requested_groups,
  gene_biotype_groups,
  transcript_biotype_groups,
  mirna_plot_groups,
  padj_cutoff,
  log2fc_cutoff,
  plan_row = list()
) {
  feature_column <- feature_id_column(results)
  groups <- list(list(
    name = "all",
    label = "All features",
    type = "all",
    results = results,
    feature_ids = as.character(results[[feature_column]])
  ))

  if (level == "mirna") {
    return(add_mirna_plot_groups(groups, results, mirna_plot_groups, padj_cutoff, log2fc_cutoff, plan_row))
  }

  biotype_column <- biotype_column_for_level(results, level)
  if (level != "transcript" && biotype_column == "") {
    return(groups)
  }

  if (level == "transcript" && "transcript_plot_group" %in% colnames(results)) {
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
        type = "transcript_novelty",
        results = subset,
        feature_ids = as.character(subset[[feature_column]])
      )
    }
  }

  if (biotype_column == "") {
    return(groups)
  }
  requested_biotypes <- if (level == "transcript") transcript_biotype_groups else gene_biotype_groups
  present_biotypes <- requested_biotypes_present(as.character(results[[biotype_column]]), requested_biotypes)
  for (biotype in present_biotypes) {
    subset <- results[safe_group_name(results[[biotype_column]]) == safe_group_name(biotype), , drop = FALSE]
    if (nrow(subset) == 0) {
      next
    }
    group_type <- ifelse(level == "transcript", "transcript_biotype", "gene_biotype")
    groups[[length(groups) + 1]] <- list(
      name = paste(group_type, safe_group_name(biotype), sep = "__"),
      label = paste(ifelse(level == "transcript", "Transcript biotype", "Gene biotype"), biotype_label(biotype), sep = ": "),
      type = group_type,
      results = subset,
      feature_ids = as.character(subset[[feature_column]])
    )
  }

  if (level == "transcript" && "transcript_plot_group" %in% colnames(results) && length(present_biotypes)) {
    requested <- requested_groups[requested_groups != "" & requested_groups != "all"]
    if (!length(requested)) {
      requested <- unique(results$transcript_plot_group[results$transcript_plot_group != ""])
    }
    for (group_name in unique(requested)) {
      novelty_subset <- results[results$transcript_plot_group == group_name, , drop = FALSE]
      if (nrow(novelty_subset) == 0) {
        next
      }
      for (biotype in present_biotypes) {
        subset <- novelty_subset[
          safe_group_name(novelty_subset[[biotype_column]]) == safe_group_name(biotype),
          ,
          drop = FALSE
        ]
        if (nrow(subset) == 0) {
          next
        }
        groups[[length(groups) + 1]] <- list(
          name = paste("transcript_novelty_biotype", group_name, safe_group_name(biotype), sep = "__"),
          label = paste(group_label(subset, group_name), biotype_label(biotype), sep = " / "),
          type = "transcript_novelty_biotype",
          results = subset,
          feature_ids = as.character(subset[[feature_column]])
        )
      }
    }
  }

  groups
}

write_plot_group_manifest <- function(path, level, contrast_id, groups) {
  ensure_parent(path)
  rows <- do.call(rbind, lapply(groups, function(group) {
    data.frame(
      level = level,
      contrast_id = contrast_id,
      plot_group = group$name,
      plot_group_type = group$type,
      plot_label = group$label,
      n_features = nrow(group$results),
      check.names = FALSE
    )
  }))
  write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
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

plot_ma_panel <- function(results, title, padj_cutoff, log2fc_cutoff) {
  ensure_columns(results, c("baseMean", "log2FoldChange", "padj"), "DESeq2 results")
  base_mean <- numeric_column(results, "baseMean")
  log2fc <- numeric_column(results, "log2FoldChange")
  padj <- numeric_column(results, "padj")
  keep <- !is.na(base_mean) & !is.na(log2fc) & !is.na(padj)
  if (!any(keep)) {
    blank_panel(title, "No valid baseMean/log2FoldChange/padj rows")
    return()
  }
  x <- log10(base_mean[keep] + 1)
  y <- log2fc[keep]
  significant <- padj[keep] < padj_cutoff & abs(y) >= log2fc_cutoff

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
}

plot_ma <- function(groups, path, title, padj_cutoff, log2fc_cutoff) {
  ensure_parent(path)
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (group in groups) {
    plot_ma_panel(group$results, paste(title, "-", group$label), padj_cutoff, log2fc_cutoff)
  }
}

read_feature_matrix <- function(path, label) {
  counts <- read_tsv(path)
  if (ncol(counts) < 2) {
    stop(label, " table needs one feature column and at least one sample column: ", path)
  }
  feature_id_column <- colnames(counts)[[1]]
  feature_ids <- counts[[feature_id_column]]
  matrix <- as.matrix(counts[, -1, drop = FALSE])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- feature_ids
  matrix
}

write_feature_matrix <- function(matrix, path) {
  output <- data.frame(feature_id = rownames(matrix), matrix, check.names = FALSE)
  ensure_parent(path)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_transformed_counts <- function(matrix, path) {
  transformed <- log2(matrix + 1)
  write_feature_matrix(transformed, path)
  transformed
}

ordered_result_feature_ids <- function(results) {
  if (nrow(results) == 0) {
    return(character())
  }
  feature_column <- feature_id_column(results)
  padj <- if ("padj" %in% colnames(results)) numeric_column(results, "padj") else rep(NA_real_, nrow(results))
  log2fc <- if ("log2FoldChange" %in% colnames(results)) numeric_column(results, "log2FoldChange") else rep(0, nrow(results))
  order_index <- order(is.na(padj), padj, -abs(log2fc), na.last = TRUE)
  unique(as.character(results[[feature_column]][order_index]))
}

top_variable_feature_ids <- function(matrix, top_n) {
  if (nrow(matrix) == 0) {
    return(character())
  }
  variances <- apply(matrix, 1, stats::var)
  names(variances)[head(order(variances, decreasing = TRUE), min(top_n, length(variances)))]
}

subset_matrix_by_ids <- function(matrix, feature_ids) {
  matched <- intersect(unique(feature_ids), rownames(matrix))
  matrix[matched, , drop = FALSE]
}

read_heatmap_feature_list <- function(path) {
  if (!file.exists(path)) {
    stop("Heatmap feature list does not exist: ", path)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) {
    return(data.frame(
      feature_id = character(),
      feature_list = character(),
      level = character(),
      contrast_id = character(),
      plot_group = character(),
      source = character(),
      check.names = FALSE
    ))
  }
  first_fields <- strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
  known_headers <- c(
    "feature_id",
    "gene_id",
    "Geneid",
    "transcript_id",
    "mirna_id",
    "id",
    "feature_list",
    "list_name",
    "level",
    "contrast_id",
    "plot_group"
  )
  has_header <- any(first_fields %in% known_headers)
  table <- read.delim(path, header = has_header, check.names = FALSE, stringsAsFactors = FALSE)
  if (!has_header) {
    colnames(table) <- paste0("V", seq_len(ncol(table)))
  }
  feature_column <- first_existing_column(
    table,
    c("feature_id", "gene_id", "Geneid", "transcript_id", "mirna_id", "id", "V1")
  )
  if (feature_column == "") {
    stop("Heatmap feature list lacks a feature-id column: ", path)
  }
  list_column <- first_existing_column(table, c("feature_list", "list_name", "name", "label"))
  default_list <- tools::file_path_sans_ext(basename(path))
  data.frame(
    feature_id = as.character(table[[feature_column]]),
    feature_list = if (list_column == "") default_list else as.character(table[[list_column]]),
    level = if ("level" %in% colnames(table)) as.character(table[["level"]]) else "",
    contrast_id = if ("contrast_id" %in% colnames(table)) as.character(table[["contrast_id"]]) else "",
    plot_group = if ("plot_group" %in% colnames(table)) as.character(table[["plot_group"]]) else "",
    source = path,
    check.names = FALSE
  )
}

read_heatmap_feature_lists <- function(paths) {
  paths <- parse_list_arg(paths)
  if (!length(paths)) {
    return(data.frame(
      feature_id = character(),
      feature_list = character(),
      level = character(),
      contrast_id = character(),
      plot_group = character(),
      source = character(),
      check.names = FALSE
    ))
  }
  rows <- lapply(paths, read_heatmap_feature_list)
  do.call(rbind, rows)
}

matching_feature_list_rows <- function(feature_lists, level, contrast_id, group) {
  if (!nrow(feature_lists)) {
    return(feature_lists)
  }
  keep <- rep(TRUE, nrow(feature_lists))
  if ("level" %in% colnames(feature_lists)) {
    keep <- keep & (feature_lists$level == "" | feature_lists$level == level)
  }
  if ("contrast_id" %in% colnames(feature_lists)) {
    keep <- keep & (feature_lists$contrast_id == "" | feature_lists$contrast_id == contrast_id)
  }
  if ("plot_group" %in% colnames(feature_lists)) {
    allowed <- c("", "all", group$name, group$type, group$label)
    keep <- keep & feature_lists$plot_group %in% allowed
  }
  feature_lists[keep, , drop = FALSE]
}

heatmap_panel <- function(
  level,
  contrast_id,
  group,
  mode,
  label,
  source,
  status,
  reason,
  requested_features,
  matched_features,
  selected_matrix
) {
  list(
    level = level,
    contrast_id = contrast_id,
    plot_group = group$name,
    plot_group_type = group$type,
    plot_label = group$label,
    heatmap_mode = mode,
    heatmap_label = label,
    source = source,
    status = status,
    reason = reason,
    n_requested_features = length(unique(requested_features)),
    n_matched_features = length(unique(matched_features)),
    n_plotted_features = nrow(selected_matrix),
    features = paste(rownames(selected_matrix), collapse = ","),
    matrix = selected_matrix
  )
}

heatmap_panels_for_group <- function(
  transformed,
  filtered,
  group,
  top_n,
  heatmap_modes,
  heatmap_feature_lists,
  fallback_mode,
  level,
  contrast_id
) {
  group_ids <- intersect(group$feature_ids, rownames(transformed))
  group_matrix <- transformed[group_ids, , drop = FALSE]
  panels <- list()

  if ("significant" %in% heatmap_modes) {
    requested <- intersect(ordered_result_feature_ids(filtered), group$feature_ids)
    matched <- intersect(requested, rownames(transformed))
    selected <- subset_matrix_by_ids(transformed, head(matched, top_n))
    status <- "ok"
    reason <- ""
    label <- "Top significant features"
    if (nrow(selected) < 2 && fallback_mode == "variable") {
      fallback_ids <- top_variable_feature_ids(group_matrix, top_n)
      selected <- subset_matrix_by_ids(transformed, fallback_ids)
      status <- if (nrow(selected) >= 2) "fallback" else "blocked"
      reason <- if (status == "fallback") {
        "Fewer than two significant features matched; plotted top variable features instead"
      } else {
        "Fewer than two significant or variable features matched"
      }
      label <- "Top significant features (fallback: top variable)"
    } else if (nrow(selected) < 2) {
      status <- "blocked"
      reason <- "Fewer than two significant features matched"
    }
    panels[[length(panels) + 1]] <- heatmap_panel(
      level,
      contrast_id,
      group,
      "significant",
      label,
      "filtered_results",
      status,
      reason,
      requested,
      matched,
      selected
    )
  }

  if ("variable" %in% heatmap_modes) {
    requested <- group_ids
    selected <- subset_matrix_by_ids(transformed, top_variable_feature_ids(group_matrix, top_n))
    status <- if (nrow(selected) >= 2) "ok" else "blocked"
    reason <- if (status == "ok") "" else "Fewer than two variable features matched"
    panels[[length(panels) + 1]] <- heatmap_panel(
      level,
      contrast_id,
      group,
      "variable",
      "Top variable features",
      "transformed_counts",
      status,
      reason,
      requested,
      rownames(group_matrix),
      selected
    )
  }

  if ("feature_list" %in% heatmap_modes) {
    rows <- matching_feature_list_rows(heatmap_feature_lists, level, contrast_id, group)
    if (nrow(rows)) {
      for (feature_list in unique(rows$feature_list)) {
        list_rows <- rows[rows$feature_list == feature_list, , drop = FALSE]
        requested <- unique(list_rows$feature_id[list_rows$feature_id != ""])
        matched <- intersect(requested, group_ids)
        selected <- subset_matrix_by_ids(transformed, head(matched, top_n))
        status <- if (nrow(selected) >= 2) "ok" else "blocked"
        reason <- if (status == "ok") "" else "Fewer than two configured-list features matched"
        panels[[length(panels) + 1]] <- heatmap_panel(
          level,
          contrast_id,
          group,
          "feature_list",
          paste("Feature list", feature_list, sep = ": "),
          paste(unique(list_rows$source), collapse = ","),
          status,
          reason,
          requested,
          matched,
          selected
        )
      }
    }
  }

  panels
}

write_heatmap_panel_manifest <- function(path, panels) {
  ensure_parent(path)
  columns <- c(
    "level",
    "contrast_id",
    "plot_group",
    "plot_group_type",
    "plot_label",
    "heatmap_mode",
    "heatmap_label",
    "source",
    "status",
    "reason",
    "n_requested_features",
    "n_matched_features",
    "n_plotted_features",
    "features"
  )
  if (!length(panels)) {
    rows <- as.data.frame(setNames(rep(list(character()), length(columns)), columns), check.names = FALSE)
  } else {
    rows <- do.call(rbind, lapply(panels, function(panel) {
      data.frame(
        level = panel$level,
        contrast_id = panel$contrast_id,
        plot_group = panel$plot_group,
        plot_group_type = panel$plot_group_type,
        plot_label = panel$plot_label,
        heatmap_mode = panel$heatmap_mode,
        heatmap_label = panel$heatmap_label,
        source = panel$source,
        status = panel$status,
        reason = panel$reason,
        n_requested_features = panel$n_requested_features,
        n_matched_features = panel$n_matched_features,
        n_plotted_features = panel$n_plotted_features,
        features = panel$features,
        check.names = FALSE
      )
    }))
  }
  write.table(rows[, columns], path, sep = "\t", quote = FALSE, row.names = FALSE)
}

transformed_counts_for_row <- function(row) {
  transformed_path <- row[["transformed_counts"]]
  if (!is.null(transformed_path) && !is.na(transformed_path) && transformed_path != "" && file.exists(transformed_path)) {
    transformed <- read_feature_matrix(transformed_path, "DESeq2 transformed counts")
    write_feature_matrix(transformed, row[["vst_tsv"]])
    return(transformed)
  }
  normalized <- read_feature_matrix(row[["normalized_counts"]], "Normalized counts")
  write_transformed_counts(normalized, row[["vst_tsv"]])
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

usable_pca_color_columns <- function(coldata, requested_columns) {
  if (is.null(coldata)) {
    return(character())
  }
  requested <- requested_columns[requested_columns != ""]
  columns <- requested[requested %in% colnames(coldata)]
  columns[vapply(
    columns,
    function(column) length(unique(na.omit(as.character(coldata[[column]])))) > 1,
    logical(1)
  )]
}

point_colors <- function(values) {
  if (is.numeric(values)) {
    if (length(unique(na.omit(values))) < 2) {
      return(rep("#1f78b4", length(values)))
    }
    palette <- grDevices::hcl.colors(100, "Viridis")
    breaks <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE), length.out = 101)
    index <- findInterval(values, breaks, all.inside = TRUE)
    colors <- palette[pmax(1, pmin(100, index))]
    colors[is.na(values)] <- "#bdbdbd"
    return(colors)
  }
  groups <- as.factor(values)
  palette <- grDevices::hcl.colors(length(levels(groups)), "Dark 3")
  colors <- palette[as.integer(groups)]
  colors[is.na(groups)] <- "#bdbdbd"
  colors
}

draw_pca_panel <- function(pca, variance, title, color_column = "", coldata = NULL) {
  values <- NULL
  colors <- rep("#1f78b4", nrow(pca$x))
  main <- title
  if (color_column != "" && !is.null(coldata) && color_column %in% colnames(coldata)) {
    values <- coldata[rownames(pca$x), color_column]
    colors <- point_colors(values)
    main <- paste(title, "-", color_column)
  }
  plot(
    pca$x[, 1],
    pca$x[, 2],
    pch = 16,
    col = colors,
    xlab = paste0("PC1 (", variance[[1]], "%)"),
    ylab = paste0("PC2 (", variance[[2]], "%)"),
    main = main
  )
  text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.7)
  if (!is.null(values) && !is.numeric(values)) {
    groups <- as.factor(values)
    palette <- grDevices::hcl.colors(length(levels(groups)), "Dark 3")
    legend("topright", legend = levels(groups), col = palette, pch = 16, bty = "n", cex = 0.8)
  }
  if (!is.null(values) && is.numeric(values)) {
    legend_values <- unique(round(stats::quantile(values, probs = c(0, 0.5, 1), na.rm = TRUE), 3))
    palette <- grDevices::hcl.colors(100, "Viridis")
    breaks <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE), length.out = 101)
    index <- findInterval(legend_values, breaks, all.inside = TRUE)
    legend(
      "topright",
      title = color_column,
      legend = legend_values,
      col = palette[pmax(1, pmin(100, index))],
      pch = 16,
      bty = "n",
      cex = 0.8
    )
  }
}

plot_pca <- function(transformed, path, title, coldata = NULL, color_columns = character()) {
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
  pca_columns <- usable_pca_color_columns(coldata, color_columns)
  ensure_parent(path)
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!length(pca_columns)) {
    draw_pca_panel(pca, variance, title)
    return()
  }
  for (column in pca_columns) {
    draw_pca_panel(pca, variance, title, column, coldata)
  }
}

plot_sample_distance <- function(transformed, path, title, coldata = NULL) {
  if (ncol(transformed) < 2 || nrow(transformed) < 2) {
    blank_pdf(path, title, "Sample-distance heatmap requires at least two samples and two features")
    return()
  }
  distances <- as.matrix(stats::dist(t(transformed)))
  if (!is.null(coldata) && "condition" %in% colnames(coldata)) {
    order_index <- order(coldata[colnames(distances), "condition"], colnames(distances))
    distances <- distances[order_index, order_index, drop = FALSE]
  }
  ensure_parent(path)
  grDevices::pdf(path)
  stats::heatmap(
    distances,
    symm = TRUE,
    margins = c(8, 8),
    main = title,
    xlab = "sample",
    ylab = "sample"
  )
  grDevices::dev.off()
}

plot_heatmap_panel <- function(transformed, title, coldata = NULL) {
  if (ncol(transformed) < 2 || nrow(transformed) < 2) {
    blank_panel(title, "Heatmap requires at least two samples and two features")
    return()
  }
  matrix <- transformed
  if (!is.null(coldata) && "condition" %in% colnames(coldata)) {
    order_index <- order(coldata[colnames(matrix), "condition"], colnames(matrix))
    matrix <- matrix[, order_index, drop = FALSE]
  }
  stats::heatmap(matrix, scale = "row", margins = c(8, 8), main = title)
}

plot_heatmap <- function(
  transformed,
  filtered,
  groups,
  path,
  heatmap_panel_tsv,
  title,
  top_n,
  heatmap_modes,
  heatmap_feature_lists,
  fallback_mode,
  level,
  contrast_id,
  coldata = NULL
) {
  panels <- unlist(
    lapply(groups, function(group) {
      heatmap_panels_for_group(
        transformed,
        filtered,
        group,
        top_n,
        heatmap_modes,
        heatmap_feature_lists,
        fallback_mode,
        level,
        contrast_id
      )
    }),
    recursive = FALSE
  )
  write_heatmap_panel_manifest(heatmap_panel_tsv, panels)
  ensure_parent(path)
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (!length(panels)) {
    blank_panel(title, "No heatmap panels requested")
    return()
  }
  for (panel in panels) {
    panel_title <- paste(title, "-", panel$plot_label, "-", panel$heatmap_label)
    if (panel$status == "blocked") {
      blank_panel(panel_title, panel$reason)
      next
    }
    plot_heatmap_panel(panel$matrix, panel_title, coldata)
  }
}

optional_plot_path <- function(row, column, fallback_name) {
  value <- row[[column]]
  if (!is.null(value) && !is.na(value) && value != "") {
    return(value)
  }
  file.path(dirname(row[["pca_pdf"]]), fallback_name)
}

write_pca_metrics <- function(path, status, reason, variance = c(NA, NA), transformed = NULL, color_columns = character()) {
  ensure_parent(path)
  metrics <- data.frame(
    status = status,
    reason = reason,
    pc1_variance_percent = ifelse(length(variance) >= 1, variance[[1]], NA),
    pc2_variance_percent = ifelse(length(variance) >= 2, variance[[2]], NA),
    n_features = ifelse(is.null(transformed), 0, nrow(transformed)),
    n_samples = ifelse(is.null(transformed), 0, ncol(transformed)),
    pca_color_columns = paste(color_columns, collapse = ","),
    check.names = FALSE
  )
  write.table(metrics, path, sep = "\t", quote = FALSE, row.names = FALSE)
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
  "pca_metrics_tsv",
  "sample_distance_pdf",
  "heatmap_pdf",
  "heatmap_panel_tsv",
  "plot_group_tsv",
  "vst_tsv",
  "n_features",
  "n_significant"
)

render_row <- function(
  row,
  top_n,
  padj_cutoff,
  log2fc_cutoff,
  transcript_plot_groups,
  gene_biotype_plot_groups,
  transcript_biotype_plot_groups,
  mirna_plot_groups,
  heatmap_modes,
  heatmap_feature_lists,
  heatmap_significant_fallback,
  pca_color_columns
) {
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
      pca_metrics_tsv = optional_plot_path(row, "pca_metrics_tsv", "pca_metrics.tsv"),
      sample_distance_pdf = optional_plot_path(row, "sample_distance_pdf", "sample_distance.pdf"),
      heatmap_pdf = row[["heatmap_pdf"]],
      heatmap_panel_tsv = optional_plot_path(row, "heatmap_panel_tsv", "heatmap_panels.tsv"),
      plot_group_tsv = optional_plot_path(row, "plot_group_tsv", "plot_groups.tsv"),
      vst_tsv = row[["vst_tsv"]],
      n_features = 0,
      n_significant = 0,
      check.names = FALSE
    ))
  }

  title <- paste(row[["project"]], row[["level"]], row[["contrast_id"]])
  results <- read_tsv(row[["results"]])
  filtered <- read_tsv(row[["filtered"]])
  transformed <- transformed_counts_for_row(row)
  coldata <- read_coldata(row[["coldata"]], colnames(transformed))
  pca_metrics_tsv <- optional_plot_path(row, "pca_metrics_tsv", "pca_metrics.tsv")
  sample_distance_pdf <- optional_plot_path(row, "sample_distance_pdf", "sample_distance.pdf")
  plot_group_tsv <- optional_plot_path(row, "plot_group_tsv", "plot_groups.tsv")
  heatmap_panel_tsv <- optional_plot_path(row, "heatmap_panel_tsv", "heatmap_panels.tsv")
  plot_groups <- plot_groups_for_results(
    results,
    row[["level"]],
    transcript_plot_groups,
    gene_biotype_plot_groups,
    transcript_biotype_plot_groups,
    mirna_plot_groups,
    padj_cutoff,
    log2fc_cutoff,
    row
  )
  write_plot_group_manifest(plot_group_tsv, row[["level"]], row[["contrast_id"]], plot_groups)
  plot_volcano(plot_groups, row[["volcano_pdf"]], paste(title, "volcano"), padj_cutoff, log2fc_cutoff)
  plot_ma(plot_groups, row[["ma_pdf"]], paste(title, "MA"), padj_cutoff, log2fc_cutoff)
  plot_pca(transformed, row[["pca_pdf"]], paste(title, "PCA"), coldata, pca_color_columns)
  if (ncol(transformed) < 2 || nrow(transformed) < 2) {
    write_pca_metrics(pca_metrics_tsv, "blocked", "PCA requires at least two samples and two features", transformed = transformed)
  } else {
    pca <- stats::prcomp(t(transformed), center = TRUE, scale. = FALSE)
    if (ncol(pca$x) < 2 || sum(pca$sdev^2) == 0) {
      write_pca_metrics(pca_metrics_tsv, "blocked", "PCA requires at least two variable dimensions", transformed = transformed)
    } else {
      variance <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
      write_pca_metrics(pca_metrics_tsv, "ok", "", variance, transformed, usable_pca_color_columns(coldata, pca_color_columns))
    }
  }
  plot_sample_distance(transformed, sample_distance_pdf, paste(title, "sample distance"), coldata)
  plot_heatmap(
    transformed,
    filtered,
    plot_groups,
    row[["heatmap_pdf"]],
    heatmap_panel_tsv,
    paste(title, "heatmap"),
    top_n,
    heatmap_modes,
    heatmap_feature_lists,
    heatmap_significant_fallback,
    row[["level"]],
    row[["contrast_id"]],
    coldata
  )

  data.frame(
    project = row[["project"]],
    level = row[["level"]],
    contrast_id = row[["contrast_id"]],
    status = "ok",
    reason = "",
    volcano_pdf = row[["volcano_pdf"]],
    ma_pdf = row[["ma_pdf"]],
    pca_pdf = row[["pca_pdf"]],
    pca_metrics_tsv = pca_metrics_tsv,
    sample_distance_pdf = sample_distance_pdf,
    heatmap_pdf = row[["heatmap_pdf"]],
    heatmap_panel_tsv = heatmap_panel_tsv,
    plot_group_tsv = plot_group_tsv,
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
  gene_biotype_plot_groups <- parse_list_arg(args[["gene_biotype_plot_groups"]])
  transcript_biotype_plot_groups <- parse_list_arg(args[["transcript_biotype_plot_groups"]])
  mirna_plot_groups <- parse_list_arg(args[["mirna_plot_groups"]])
  heatmap_modes <- parse_list_arg(args[["heatmap_modes"]])
  if (!length(heatmap_modes)) {
    heatmap_modes <- c("significant", "variable")
  }
  invalid_heatmap_modes <- setdiff(heatmap_modes, c("significant", "variable", "feature_list"))
  if (length(invalid_heatmap_modes)) {
    stop("Unsupported --heatmap-modes value(s): ", paste(invalid_heatmap_modes, collapse = ", "))
  }
  heatmap_feature_lists <- read_heatmap_feature_lists(args[["heatmap_feature_lists"]])
  if (nrow(heatmap_feature_lists) > 0 && !"feature_list" %in% heatmap_modes) {
    heatmap_modes <- c(heatmap_modes, "feature_list")
  }
  heatmap_significant_fallback <- args[["heatmap_significant_fallback"]]
  if (!heatmap_significant_fallback %in% c("variable", "none")) {
    stop("--heatmap-significant-fallback must be 'variable' or 'none'")
  }
  pca_color_columns <- parse_list_arg(args[["pca_color_columns"]])

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
      render_row(
        row,
        top_n,
        padj_cutoff,
        log2fc_cutoff,
        transcript_plot_groups,
        gene_biotype_plot_groups,
        transcript_biotype_plot_groups,
        mirna_plot_groups,
        heatmap_modes,
        heatmap_feature_lists,
        heatmap_significant_fallback,
        pca_color_columns
      ),
      error = function(err) {
        pca_metrics_tsv <- optional_plot_path(row, "pca_metrics_tsv", "pca_metrics.tsv")
        try(write_pca_metrics(pca_metrics_tsv, "failed", conditionMessage(err)), silent = TRUE)
        data.frame(
          project = row[["project"]],
          level = row[["level"]],
          contrast_id = row[["contrast_id"]],
          status = "failed",
          reason = conditionMessage(err),
          volcano_pdf = row[["volcano_pdf"]],
          ma_pdf = row[["ma_pdf"]],
          pca_pdf = row[["pca_pdf"]],
          pca_metrics_tsv = pca_metrics_tsv,
          sample_distance_pdf = optional_plot_path(row, "sample_distance_pdf", "sample_distance.pdf"),
          heatmap_pdf = row[["heatmap_pdf"]],
          heatmap_panel_tsv = optional_plot_path(row, "heatmap_panel_tsv", "heatmap_panels.tsv"),
          plot_group_tsv = optional_plot_path(row, "plot_group_tsv", "plot_groups.tsv"),
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
