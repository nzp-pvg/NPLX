source("script/R/00_project_config.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

ensure_project_dirs()
dir.create(project_path("results", "cellchat_npl_module"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figures"), recursive = TRUE, showWarnings = FALSE)

core_color <- "#005493"
adj_color <- "#F5BD4D"
obs_color <- "#C34062"
low_color <- adj_color
high_color <- core_color
edge_gray <- "#5C5C5C"

receiver_levels <- c("Endothelial", "SMC", "Fibroblast/Mesenchymal", "T cell")
receiver_map <- c(
  ENDOTHELIAL = "Endothelial",
  SMC = "SMC",
  FIBROBLAST = "Fibroblast/Mesenchymal",
  T_CELL = "T cell"
)

pair_genes <- c(
  "SPP1", "CD44", "ITGAV", "ITGB1", "ITGB3", "ITGA5", "ITGA8", "ITGAX", "ITGB5",
  "APOE", "LRP1", "LRP8", "LRP10", "LRP12", "C1QTNF1",
  "LGALS9", "HAVCR2", "PTPRC", "PTPN6",
  "CCL2", "CCR2", "CCR5", "CXCL8", "CXCR1", "CXCR2", "CXCL12", "CXCR4",
  "TNF", "TNFRSF1A", "TNFRSF1B", "IL1B", "IL1R1", "IL1RAP",
  "MIF", "CD74", "CXCR4", "ACKR3",
  "NAMPT", "ITGA4", "ITGB1",
  "TGFB1", "TGFBR1", "TGFBR2", "TGFBR3",
  "PDGFA", "PDGFB", "PDGFC", "PDGFD", "PDGFRA", "PDGFRB",
  "VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1", "FLT4",
  "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "JAG1", "JAG2",
  "DLL1", "DLL3", "DLL4", "CX3CL1", "CX3CR1"
)

compact_genes <- c("NPL", "FABP5", "GPNMB", "APOC1", "PLA2G7", "SPP1", "CD36", "CYP27A1", "APOE")

find_triplet_dir <- function(root_dir) {
  candidates <- list.files(root_dir, pattern = "matrix.mtx.gz$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) != 1) {
    stop("Expected exactly one matrix.mtx.gz under ", root_dir, ", found ", length(candidates))
  }
  dirname(candidates[[1]])
}

read_10x_triplet <- function(triplet_dir) {
  matrix_path <- file.path(triplet_dir, "matrix.mtx.gz")
  feature_path <- file.path(triplet_dir, "features.tsv.gz")
  barcode_path <- file.path(triplet_dir, "barcodes.tsv.gz")
  mat <- readMM(matrix_path)
  features <- read.delim(gzfile(feature_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(barcode_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  rownames(mat) <- make.unique(features$V2)
  colnames(mat) <- barcodes$V1
  mat
}

normalize_log1p <- function(mat) {
  lib <- Matrix::colSums(mat)
  lib[lib == 0] <- 1
  norm <- t(t(mat) / lib * 10000)
  log1p(norm)
}

split_symbols <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  parts <- unlist(strsplit(x, "\\s*,\\s*|\\s*;\\s*|\\s*\\|\\s*"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}

parse_complex_name <- function(x) {
  if (is.na(x) || !nzchar(x)) {
    return(character(0))
  }
  x <- gsub("_", ",", x)
  split_symbols(x)
}

gene_group_summary <- function(norm_mat, cell_idx) {
  sub <- norm_mat[, cell_idx, drop = FALSE]
  data.frame(
    gene_symbol = rownames(sub),
    mean_expr = Matrix::rowMeans(sub),
    pct_expr = Matrix::rowSums(sub > 0) / ncol(sub),
    stringsAsFactors = FALSE
  )
}

complex_activity <- function(summary_df, genes) {
  genes <- intersect(genes, summary_df$gene_symbol)
  if (length(genes) == 0) {
    return(NA_real_)
  }
  x <- summary_df[match(genes, summary_df$gene_symbol), , drop = FALSE]
  gene_weight <- pmax(x$mean_expr, 0) * pmax(x$pct_expr, 0)
  gene_weight <- gene_weight[is.finite(gene_weight)]
  gene_weight <- gene_weight[gene_weight > 0]
  if (length(gene_weight) == 0) {
    return(0)
  }
  exp(mean(log(gene_weight + 1e-8)))
}

safe_paired_cohen_d <- function(high, low) {
  d <- high - low
  d <- d[is.finite(d)]
  if (length(d) < 2) {
    return(NA_real_)
  }
  s <- stats::sd(d, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(NA_real_)
  }
  mean(d, na.rm = TRUE) / s
}

interaction_theme <- function(pathway_name, ligand, receptor, interaction_name) {
  txt <- paste(pathway_name, ligand, receptor, interaction_name, sep = " | ")
  rules <- list(
    "SPP1-CD44/integrin" = c("SPP1", "CD44", "ITGAV", "ITGB1", "ITGB3", "ITGA5", "ITGA8", "ITGAX", "ITGB5"),
    "ApoE-LRP1" = c("APOE", "LRP1", "LRP8", "LRP10", "LRP12"),
    "Galectin" = c("LGALS9", "LGALS3", "HAVCR2"),
    "CCL2-CCR2" = c("CCL2", "CCR2", "CCR5"),
    "CXCL-CXCR" = c("CXCL", "CXCR", "ACKR3"),
    "TNF" = c("TNF", "TNFRSF1A", "TNFRSF1B"),
    "IL1" = c("IL1", "IL1R1", "IL1RAP"),
    "MIF" = c("MIF", "CD74", "CXCR4", "ACKR3"),
    "NAMPT/Visfatin" = c("NAMPT"),
    "TGF-beta" = c("TGFB", "TGFB", "TGFBR", "ACVR"),
    "PDGF" = c("PDGF", "PDGFRA", "PDGFRB"),
    "VEGF" = c("VEGF", "KDR", "FLT1", "FLT4"),
    "NOTCH/JAG" = c("NOTCH", "JAG", "DLL")
  )
  for (nm in names(rules)) {
    pat <- paste(rules[[nm]], collapse = "|")
    if (grepl(pat, txt, ignore.case = TRUE)) {
      return(nm)
    }
  }
  "other"
}

interaction_module <- function(pathway_theme, ligand, receptor, interaction_name) {
  txt <- paste(pathway_theme, ligand, receptor, interaction_name, sep = " | ")
  rules <- list(
    "ECM / vascular remodeling" = c("SPP1-CD44/integrin", "TGF-beta", "PDGF", "VEGF", "NOTCH/JAG"),
    "Lipid / foam biology" = c("ApoE-LRP1", "Galectin"),
    "Immune / inflammatory signaling" = c("CCL2-CCR2", "CXCL-CXCR", "IL1", "MIF", "TNF")
  )
  for (nm in names(rules)) {
    pat <- paste(rules[[nm]], collapse = "|")
    if (grepl(pat, txt, ignore.case = TRUE)) {
      return(nm)
    }
  }
  "Other"
}

make_node_layout <- function(nodes) {
  n <- length(nodes)
  theta <- seq(pi / 2, 2 * pi + pi / 2, length.out = n + 1)[1:n]
  data.frame(
    node = nodes,
    x = cos(theta),
    y = sin(theta),
    angle = theta,
    stringsAsFactors = FALSE
  )
}

make_edges_for_plot <- function(edge_df, node_layout) {
  edges <- edge_df %>%
    mutate(
      x = node_layout$x[match(source, node_layout$node)],
      y = node_layout$y[match(source, node_layout$node)],
      xend = node_layout$x[match(target, node_layout$node)],
      yend = node_layout$y[match(target, node_layout$node)]
    )
  edges
}

load_cellchat_db <- function() {
  db_cache <- project_path("data", "raw", "reference", "CellChatDB.human.rda")
  dir.create(dirname(db_cache), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(db_cache)) {
    url <- "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.human.rda"
    download.file(url, db_cache, quiet = TRUE, mode = "wb")
  }
  load(db_cache, envir = environment())
  CellChatDB.human
}

build_sample_metadata <- function(pheno, cell_typing) {
  meta <- pheno[, c("sample_id", "group_label", "title"), drop = FALSE]
  meta$patient_id <- sub(".*Patient ([0-9]+).*", "P\\1", meta$title)
  meta$sample_name <- cell_typing$sample_name[match(meta$sample_id, cell_typing$sample_id)]
  meta
}

summarize_communications <- function(interactions, per_sample_groups, target_pairs, db_genes, sample_ids, sample_meta) {
  sample_rows <- list()

  for (sample_id in sample_ids) {
    sample_groups <- per_sample_groups[[sample_id]]
    if (is.null(sample_groups)) {
      next
    }
    sample_name <- sample_meta$sample_name[match(sample_id, sample_meta$sample_id)]
    patient_id <- sample_meta$patient_id[match(sample_id, sample_meta$sample_id)]

    for (interaction_idx in seq_len(nrow(interactions))) {
      row <- interactions[interaction_idx, , drop = FALSE]
      ligand_genes <- parse_complex_name(row$ligand.symbol)
      receptor_genes <- parse_complex_name(row$receptor.symbol)
      if (length(ligand_genes) == 0 || length(receptor_genes) == 0) {
        next
      }
      if (!all(c(ligand_genes, receptor_genes) %in% db_genes)) {
        next
      }

      for (pair in target_pairs) {
        src_group <- pair$source
        tgt_group <- pair$target
        if (!src_group %in% names(sample_groups) || !tgt_group %in% names(sample_groups)) {
          next
        }

        src_sum <- sample_groups[[src_group]]
        tgt_sum <- sample_groups[[tgt_group]]
        ligand_score <- complex_activity(src_sum, ligand_genes)
        receptor_score <- complex_activity(tgt_sum, receptor_genes)
        comm_score <- sqrt(pmax(ligand_score, 0) * pmax(receptor_score, 0))
        if (!is.finite(comm_score)) {
          comm_score <- 0
        }

        sample_rows[[length(sample_rows) + 1]] <- data.frame(
          sample_id = sample_id,
          sample_name = sample_name,
          patient_id = patient_id,
          direction = pair$direction,
          macrophage_group = pair$macrophage_group,
          macrophage_side = ifelse(pair$direction == "outgoing", "sender", "receiver"),
          other_cell_type = pair$receiver_cell_type,
          source_group = src_group,
          target_group = tgt_group,
          interaction_name = row$interaction_name,
          pathway_name = row$pathway_name,
          ligand = row$ligand.symbol,
          receptor = row$receptor.symbol,
          ligand_genes = paste(ligand_genes, collapse = "|"),
          receptor_genes = paste(receptor_genes, collapse = "|"),
          ligand_score = ligand_score,
          receptor_score = receptor_score,
          communication_score = comm_score,
          pathway_theme = interaction_theme(row$pathway_name, row$ligand.symbol, row$receptor.symbol, row$interaction_name),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, sample_rows)
}

aggregate_summary <- function(sample_scores) {
  wide <- sample_scores %>%
    select(
      sample_id, sample_name, patient_id, direction, macrophage_side, other_cell_type,
      interaction_name, pathway_name, ligand, receptor, ligand_genes, receptor_genes,
      pathway_theme, macrophage_group, communication_score
    ) %>%
    distinct() %>%
    pivot_wider(
      names_from = macrophage_group,
      values_from = communication_score
    )

  wide %>%
    group_by(direction, macrophage_side, other_cell_type, interaction_name, pathway_name, ligand, receptor, ligand_genes, receptor_genes, pathway_theme) %>%
    summarise(
      n_pairs = sum(!is.na(high) & !is.na(low)),
      mean_score_high = mean(high, na.rm = TRUE),
      mean_score_low = mean(low, na.rm = TRUE),
      median_score_high = median(high, na.rm = TRUE),
      median_score_low = median(low, na.rm = TRUE),
      delta_score = mean_score_high - mean_score_low,
      median_delta_score = median_score_high - median_score_low,
      log2fc_score = log2((mean_score_high + 1e-8) / (mean_score_low + 1e-8)),
      paired_wilcoxon_p = tryCatch(wilcox.test(high, low, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_),
      paired_cohen_d = safe_paired_cohen_d(high, low),
      .groups = "drop"
    ) %>%
    mutate(
      mean_score = (mean_score_high + mean_score_low) / 2,
      abs_delta = abs(delta_score)
    ) %>%
    arrange(desc(abs_delta), desc(mean_score))
}

main <- function() {
  pheno <- load_pheno("GSE159677")
  cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
  macro_cont <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_cell_table.tsv"))
  sample_meta <- build_sample_metadata(pheno, cell_typing)
  db <- load_cellchat_db()
  interactions <- db$interaction

  core_sample_ids <- sample_meta$sample_id[sample_meta$group_label == "core"]
  core_sample_ids <- unique(core_sample_ids)

  sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[grepl("^GSM", basename(sample_dirs))]

  per_sample_groups <- list()
  total_score_rows <- list()

  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    sample_id <- sub("_.*$", "", sample_name)
    if (!sample_id %in% core_sample_ids) {
      next
    }

    counts <- read_10x_triplet(find_triplet_dir(sample_dir))
    norm <- normalize_log1p(counts)

    cell_meta <- cell_typing[cell_typing$sample_name == sample_name, , drop = FALSE]
    cell_meta$barcode <- as.character(cell_meta$barcode)
    cell_meta$location <- sample_meta$group_label[match(sample_id, sample_meta$sample_id)]
    cell_meta$patient_id <- sample_meta$patient_id[match(sample_id, sample_meta$sample_id)]

    cell_ids <- intersect(colnames(norm), cell_meta$barcode)
    cell_meta <- cell_meta[match(cell_ids, cell_meta$barcode), , drop = FALSE]
    names(cell_meta$barcode) <- NULL
    cell_meta$cell_id <- paste(sample_name, cell_meta$barcode, sep = "|")

    macro_meta <- macro_cont[macro_cont$sample_id == sample_id & macro_cont$location == "core", , drop = FALSE]
    macro_meta$barcode <- sub("^.*\\|", "", macro_meta$cell_id)
    macro_scores <- macro_meta$npl_module_score
    names(macro_scores) <- macro_meta$barcode
    macro_cell_barcodes <- intersect(cell_meta$barcode[cell_meta$predicted_cell_type == "MACROPHAGE"], names(macro_scores))
    if (length(macro_cell_barcodes) == 0) {
      next
    }
    macro_cell_scores <- macro_scores[macro_cell_barcodes]
    high_cut <- quantile(macro_cell_scores, probs = 0.70, na.rm = TRUE)
    low_cut <- quantile(macro_cell_scores, probs = 0.30, na.rm = TRUE)
    macro_groups <- rep(NA_character_, length(macro_cell_barcodes))
    names(macro_groups) <- macro_cell_barcodes
    macro_groups[macro_cell_scores >= high_cut] <- "high"
    macro_groups[macro_cell_scores <= low_cut] <- "low"

    cell_meta$npl_group <- NA_character_
    cell_meta$npl_group[cell_meta$predicted_cell_type == "MACROPHAGE" & cell_meta$barcode %in% names(macro_groups)] <-
      macro_groups[cell_meta$barcode[cell_meta$predicted_cell_type == "MACROPHAGE" & cell_meta$barcode %in% names(macro_groups)]]

    group_defs <- list(
      `NPL_module_high_macrophage` = cell_meta$barcode[cell_meta$predicted_cell_type == "MACROPHAGE" & cell_meta$npl_group == "high"],
      `NPL_module_low_macrophage` = cell_meta$barcode[cell_meta$predicted_cell_type == "MACROPHAGE" & cell_meta$npl_group == "low"],
      `Endothelial` = cell_meta$barcode[cell_meta$predicted_cell_type == "ENDOTHELIAL"],
      `SMC` = cell_meta$barcode[cell_meta$predicted_cell_type == "SMC"],
      `Fibroblast/Mesenchymal` = cell_meta$barcode[cell_meta$predicted_cell_type == "FIBROBLAST"],
      `T cell` = cell_meta$barcode[cell_meta$predicted_cell_type == "T_CELL"]
    )

    group_defs <- group_defs[lengths(group_defs) > 0]
    group_summaries <- list()
    for (nm in names(group_defs)) {
      idx <- match(group_defs[[nm]], colnames(norm))
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) {
        next
      }
      group_summaries[[nm]] <- gene_group_summary(norm, idx)
    }
    per_sample_groups[[sample_id]] <- group_summaries

    total_score_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      sample_name = sample_name,
      patient_id = sample_meta$patient_id[match(sample_id, sample_meta$sample_id)],
      macrophage_group = rep(c("high", "low"), each = 2 * 4),
      direction = rep(c("outgoing", "incoming"), times = 2 * 4),
      receiver_cell_type = rep(receiver_levels, times = 2),
      total_communication_score = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  target_pairs <- expand.grid(
    macrophage_group = c("high", "low"),
    direction = c("outgoing", "incoming"),
    receiver_cell_type = receiver_levels,
    stringsAsFactors = FALSE
  )
  target_pairs <- target_pairs %>%
    mutate(
      source = ifelse(direction == "outgoing",
                      ifelse(macrophage_group == "high", "NPL_module_high_macrophage", "NPL_module_low_macrophage"),
                      receiver_cell_type),
      target = ifelse(direction == "outgoing",
                      receiver_cell_type,
                      ifelse(macrophage_group == "high", "NPL_module_high_macrophage", "NPL_module_low_macrophage"))
    ) %>%
    select(macrophage_group, direction, receiver_cell_type, source, target)

  db_genes <- unique(c(interactions$ligand.symbol, interactions$receptor.symbol))
  db_genes <- unlist(lapply(db_genes, parse_complex_name))
  db_genes <- unique(db_genes)

  sample_scores <- summarize_communications(
    interactions = interactions,
    per_sample_groups = per_sample_groups,
    target_pairs = split(target_pairs, seq_len(nrow(target_pairs))),
    db_genes = db_genes,
    sample_ids = names(per_sample_groups),
    sample_meta = sample_meta
  )

  if (is.null(sample_scores) || nrow(sample_scores) == 0) {
    stop("No communication scores were generated.")
  }

  all_summary <- aggregate_summary(sample_scores)
  all_summary <- all_summary %>%
    mutate(
      interaction_label = paste0(ligand, " → ", receptor),
      edge_label = ifelse(
        macrophage_side == "sender",
        paste0("Macrophage → ", other_cell_type),
        paste0(other_cell_type, " → Macrophage")
      )
    )

  # Overall total score summaries for panel A
  sample_scores <- sample_scores %>%
    mutate(
      edge_class = case_when(
        grepl("^NPL_module_", source_group) ~ "outgoing",
        grepl("^NPL_module_", target_group) ~ "incoming",
        TRUE ~ "other"
      ),
      edge_group = ifelse(macrophage_group == "high", "high", "low")
    )

  overall_total <- sample_scores %>%
    group_by(sample_id, sample_name, patient_id, macrophage_group, direction) %>%
    summarise(total_communication_score = sum(communication_score, na.rm = TRUE), .groups = "drop")

  panel_a_stats <- overall_total %>%
    group_by(direction) %>%
    summarise(
      n_patients = n_distinct(patient_id),
      delta_mean = mean(total_communication_score[macrophage_group == "high"], na.rm = TRUE) -
        mean(total_communication_score[macrophage_group == "low"], na.rm = TRUE),
      paired_p = tryCatch(wilcox.test(
        total_communication_score[macrophage_group == "high"],
        total_communication_score[macrophage_group == "low"],
        paired = TRUE, exact = FALSE
      )$p.value, error = function(e) NA_real_),
      cohen_d = safe_paired_cohen_d(
        total_communication_score[macrophage_group == "high"],
        total_communication_score[macrophage_group == "low"]
      ),
      y_pos = max(total_communication_score, na.rm = TRUE) * 1.05,
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(
        "n = ", n_patients,
        "\n\u0394mean = ", sprintf("%.2f", delta_mean),
        "\nP = ", ifelse(is.na(paired_p), "NA", formatC(paired_p, format = "e", digits = 2)),
        "\nd = ", ifelse(is.na(cohen_d), "NA", sprintf("%.2f", cohen_d))
      )
    )

  write.csv(all_summary, project_path("results", "cellchat_npl_module", "cellchat_all_interactions.csv"), row.names = FALSE)
  write.csv(
    all_summary %>% arrange(desc(abs_delta), paired_wilcoxon_p),
    project_path("results", "cellchat_npl_module", "cellchat_top_pairs.csv"),
    row.names = FALSE
  )

  pathway_wide <- sample_scores %>%
    group_by(sample_id, sample_name, patient_id, direction, macrophage_side, other_cell_type, pathway_name, pathway_theme, macrophage_group) %>%
    summarise(total_score = sum(communication_score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = macrophage_group, values_from = total_score)

  pathway_summary <- pathway_wide %>%
    group_by(direction, macrophage_side, other_cell_type, pathway_name, pathway_theme) %>%
    summarise(
      n_pairs = sum(!is.na(high) & !is.na(low)),
      score_high = mean(high, na.rm = TRUE),
      score_low = mean(low, na.rm = TRUE),
      median_score_high = median(high, na.rm = TRUE),
      median_score_low = median(low, na.rm = TRUE),
      delta_score = score_high - score_low,
      median_delta_score = median_score_high - median_score_low,
      log2fc_score = log2((score_high + 1e-8) / (score_low + 1e-8)),
      paired_wilcoxon_p = tryCatch(wilcox.test(high, low, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_),
      paired_cohen_d = safe_paired_cohen_d(high, low),
      .groups = "drop"
    ) %>%
    mutate(abs_delta = abs(delta_score)) %>%
    arrange(desc(abs_delta))

  write.csv(pathway_summary, project_path("results", "cellchat_npl_module", "cellchat_pathway_summary.csv"), row.names = FALSE)
  write.csv(sample_scores, project_path("results", "cellchat_npl_module", "cellchat_sample_level_scores.csv"), row.names = FALSE)
  write.csv(overall_total, project_path("results", "cellchat_npl_module", "cellchat_total_strength_by_sample.csv"), row.names = FALSE)

  # Panel A: overall strength
  p_a <- overall_total %>%
    mutate(
      macrophage_group = factor(macrophage_group, levels = c("low", "high")),
      direction = factor(direction, levels = c("outgoing", "incoming")),
      patient_id = factor(patient_id, levels = c("P1", "P2", "P3"))
    ) %>%
    ggplot(aes(x = macrophage_group, y = total_communication_score, group = patient_id, color = patient_id)) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(size = 2.2) +
    geom_boxplot(aes(fill = macrophage_group), alpha = 0.20, width = 0.55, outlier.shape = NA, color = "black") +
    facet_wrap(~direction, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(low = low_color, high = high_color)) +
    scale_color_manual(values = c(P1 = obs_color, P2 = "#2DA44E", P3 = "#6C8FF0")) +
    labs(
      title = "A  Global communication activity in NPL-module-high versus low macrophages",
      subtitle = "Matched plaque-core samples only; sample-aware paired comparison (n = 3 patients)",
      x = NULL,
      y = "Total communication score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "black", color = "black"),
      strip.text = element_text(color = "white", face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    ) +
    geom_label(
      data = panel_a_stats,
      aes(x = 1.5, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      linewidth = 0.25,
      fill = "white",
      alpha = 0.85,
      label.padding = unit(0.12, "lines")
    )

  # Panel B: circle/network plot of aggregated communication
  network_edges <- sample_scores %>%
    group_by(sample_id, patient_id, macrophage_group, direction, source_group, target_group) %>%
    summarise(score = sum(communication_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(macrophage_group, direction, source_group, target_group) %>%
    summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
  network_edges <- network_edges %>%
    mutate(
      source = source_group,
      target = target_group,
      edge_color = ifelse(macrophage_group == "high", high_color, low_color),
      edge_type = direction
    )
  network_nodes <- c("NPL_module_high_macrophage", "NPL_module_low_macrophage", receiver_levels)
  node_layout <- make_node_layout(network_nodes)
  node_layout <- node_layout %>%
    mutate(
      node_size = ifelse(node == "NPL_module_high_macrophage", 11.5, ifelse(node == "NPL_module_low_macrophage", 10, 8.6))
    )
  network_edges_plot <- make_edges_for_plot(network_edges, node_layout)

  p_b <- ggplot() +
    geom_curve(
      data = network_edges_plot,
      aes(
        x = x, y = y, xend = xend, yend = yend,
        linewidth = score,
        color = macrophage_group,
        linetype = edge_type
      ),
      curvature = 0.18,
      arrow = arrow(length = unit(0.11, "inches"), type = "closed"),
      alpha = 0.8
    ) +
    geom_point(
      data = node_layout,
      aes(x = x, y = y, fill = node, size = node_size),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    geom_text(
      data = node_layout,
      aes(
        x = x * 1.14,
        y = y * 1.14,
        label = dplyr::case_when(
          node == "NPL_module_high_macrophage" ~ "NPL-high\nmacrophage",
          node == "NPL_module_low_macrophage" ~ "NPL-low\nmacrophage",
          node == "Fibroblast/Mesenchymal" ~ "Fibroblast/\nmesenchymal",
          TRUE ~ gsub("_", " ", node)
        )
      ),
      size = 3.4,
      fontface = "bold"
    ) +
    scale_color_manual(values = c(high = high_color, low = low_color)) +
    scale_fill_manual(values = c(
      NPL_module_high_macrophage = high_color,
      NPL_module_low_macrophage = low_color,
      Endothelial = "#A6CEE3",
      SMC = "#B2DF8A",
      `Fibroblast/Mesenchymal` = "#FDBF6F",
      `T cell` = "#CAB2D6"
    )) +
    scale_size_identity() +
    scale_linewidth(range = c(0.4, 2.4)) +
    labs(
      title = "B  Communication network centered on NPL-high macrophages",
      subtitle = "Sample-averaged network; edge width reflects communication score and line type reflects direction",
      color = "Macrophage group",
      linewidth = "Mean score",
      linetype = "Direction"
    ) +
    guides(fill = "none", size = "none") +
    coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), expand = FALSE) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9, hjust = 0),
      legend.box = "vertical"
    )

  # Panel C: bubble plot of top interactions
  top_pairs <- all_summary %>%
    filter(abs_delta > 0) %>%
    mutate(
      module_group = mapply(interaction_module, pathway_theme, ligand, receptor, interaction_name, USE.NAMES = FALSE)
    ) %>%
    filter(module_group != "Other") %>%
    group_by(direction, module_group) %>%
    slice_max(order_by = abs_delta, n = 4, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      interaction_label = interaction_label,
      pair_label = other_cell_type,
      module_group = factor(module_group, levels = c("ECM / vascular remodeling", "Lipid / foam biology", "Immune / inflammatory signaling")),
      direction = factor(direction, levels = c("outgoing", "incoming"))
    )

  top_pairs$interaction_label <- factor(
    top_pairs$interaction_label,
    levels = rev(unique(top_pairs$interaction_label[order(top_pairs$module_group, top_pairs$abs_delta)]))
  )
  top_pairs$pair_label <- factor(top_pairs$pair_label, levels = unique(top_pairs$pair_label))

  p_c <- ggplot(top_pairs, aes(x = pair_label, y = interaction_label)) +
    geom_point(aes(size = mean_score, color = delta_score), alpha = 0.9) +
    facet_grid(module_group ~ direction, scales = "free_y", space = "free_y") +
    scale_color_gradient2(low = low_color, mid = "white", high = high_color, midpoint = 0) +
    scale_size(range = c(1.5, 7)) +
    labs(
      title = "C  Module-organized ligand-receptor pairs linked to NPL-high macrophages",
      subtitle = "Bubble color indicates high minus low communication score; size indicates mean score",
      x = NULL,
      y = NULL,
      color = "Delta score",
      size = "Mean score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "black", color = "black"),
      strip.text = element_text(color = "white", face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  # Panel D: pathway heatmap
  pathway_top <- pathway_summary %>%
    group_by(pathway_name) %>%
    summarise(abs_delta_total = sum(abs(delta_score), na.rm = TRUE), .groups = "drop") %>%
    slice_max(order_by = abs_delta_total, n = 8, with_ties = FALSE)

  pathway_plot_df <- pathway_summary %>%
    semi_join(pathway_top, by = "pathway_name") %>%
    mutate(
      pathway_name = factor(pathway_name, levels = rev(pathway_top$pathway_name)),
      other_cell_type = factor(other_cell_type, levels = receiver_levels),
      direction = factor(direction, levels = c("outgoing", "incoming"))
    )

  p_d <- ggplot(pathway_plot_df, aes(x = other_cell_type, y = pathway_name, fill = delta_score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    facet_wrap(~direction, ncol = 1) +
    scale_fill_gradient2(low = low_color, mid = "white", high = high_color, midpoint = 0) +
    labs(
      title = "D  Strongest pathway-level communication shifts",
      subtitle = "Top 8 pathways ranked by total absolute shift; values are high minus low",
      x = NULL,
      y = NULL,
      fill = "Delta score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "black", color = "black"),
      strip.text = element_text(color = "white", face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  pathway_summary_plot <- pathway_summary %>%
    group_by(direction, pathway_name) %>%
    summarise(delta_total = mean(delta_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(direction) %>%
    slice_max(order_by = abs(delta_total), n = 6, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      direction = factor(direction, levels = c("outgoing", "incoming")),
      pathway_name = factor(pathway_name, levels = rev(unique(pathway_name[order(direction, abs(delta_total))])))
    )

  p_e <- ggplot(pathway_summary_plot, aes(x = delta_total, y = pathway_name, fill = delta_total > 0)) +
    geom_col(width = 0.72) +
    facet_wrap(~direction, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(`TRUE` = high_color, `FALSE` = low_color), guide = "none") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    labs(
      title = "E  Pathway summary of gains and losses with NPL-high macrophages",
      subtitle = "Top direction-specific pathways by mean delta score",
      x = "Mean delta score (high - low)",
      y = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "black", color = "black"),
      strip.text = element_text(color = "white", face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )

  combined <- (p_a | p_b) / p_c / (p_d | p_e) + plot_layout(heights = c(1.02, 1.2, 1.05))
  combined <- combined + plot_annotation(
    title = "Communication remodeling associated with NPL-module-high macrophages in GSE159677",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 10))
    )
  )

  ggsave(
    filename = project_path("figures", "SuppFig_cellchat_npl_module.png"),
    plot = combined,
    width = 22,
    height = 26,
    dpi = 300
  )

  cat("CellChat-like communication analysis complete.\n")
  cat("Saved results to:", project_path("results", "cellchat_npl_module"), "\n")
  cat("Saved figure to:", project_path("figures", "SuppFig_cellchat_npl_module.png"), "\n")
}

main()
