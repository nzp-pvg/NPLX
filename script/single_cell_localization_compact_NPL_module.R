resolve_this_file <- function() {
  frame_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(frame_file) && nzchar(frame_file)) {
    return(normalizePath(frame_file, winslash = "/", mustWork = TRUE))
  }
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(cmd_file) > 0) {
    return(normalizePath(sub("^--file=", "", cmd_file[[1]]), winslash = "/", mustWork = TRUE))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

this_file <- resolve_this_file()
root_dir <- normalizePath(file.path(dirname(this_file), "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(root_dir, "script", "R", "00_project_config.R"))

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

ensure_project_dirs()
dir.create(project_path("results", "figure3"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figures", "main"), recursive = TRUE, showWarnings = FALSE)

compact_genes <- c("NPL", "FABP5", "GPNMB", "APOC1", "PLA2G7", "SPP1", "CD36", "CYP27A1", "APOE")
compact_gene_order_plot <- c("NPL", "FABP5", "GPNMB", "APOC1", "APOE", "CD36", "CYP27A1", "PLA2G7", "SPP1")

celltype_levels <- c("B_CELL", "ENDOTHELIAL", "FIBROBLAST", "MACROPHAGE", "MAST", "SMC", "T_CELL", "UNRESOLVED")
celltype_labels <- c(
  B_CELL = "B cell",
  ENDOTHELIAL = "Endothelial",
  FIBROBLAST = "Fibroblast/\nmesenchymal",
  MACROPHAGE = "Macrophage",
  MAST = "Mast",
  SMC = "SMC",
  T_CELL = "T cell",
  UNRESOLVED = "Unresolved"
)

state_levels <- c("MACROPHAGE_C1Q", "MACROPHAGE_INFLAMMATORY", "MACROPHAGE_IFN", "MACROPHAGE_FOAM_TREM2")
state_labels <- c(
  MACROPHAGE_C1Q = "C1Q-like",
  MACROPHAGE_INFLAMMATORY = "Inflammatory",
  MACROPHAGE_IFN = "IFN-like",
  MACROPHAGE_FOAM_TREM2 = "FOAM/TREM2-like"
)

core_color <- "#005493"
adj_color <- "#E6AB2A"
accent_color <- "#C34062"
mid_blue <- "#4582B0"
mid_gray <- "#8F979A"
light_gray <- "#D2D2D2"

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

score_gene_set <- function(norm_mat, genes) {
  genes <- intersect(genes, rownames(norm_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(norm_mat)))
  }
  Matrix::colMeans(norm_mat[genes, , drop = FALSE])
}

paired_stats <- function(df, value_col, pair_col = "patient_id", group_col = "location", group_a = "adjacent", group_b = "core") {
  required <- c(value_col, pair_col, group_col)
  if (!all(required %in% colnames(df))) {
    return(data.frame(
      n_pairs = 0,
      mean_adjacent = NA_real_,
      mean_core = NA_real_,
      delta_core_minus_adjacent = NA_real_,
      wilcoxon_p = NA_real_,
      paired_cohen_d = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  wide <- reshape(df[, required, drop = FALSE], idvar = pair_col, timevar = group_col, direction = "wide")
  a_col <- paste0(value_col, ".", group_a)
  b_col <- paste0(value_col, ".", group_b)
  if (!(a_col %in% colnames(wide)) || !(b_col %in% colnames(wide))) {
    return(data.frame(
      n_pairs = 0,
      mean_adjacent = mean(df[[value_col]][df[[group_col]] == group_a], na.rm = TRUE),
      mean_core = mean(df[[value_col]][df[[group_col]] == group_b], na.rm = TRUE),
      delta_core_minus_adjacent = NA_real_,
      wilcoxon_p = NA_real_,
      paired_cohen_d = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  a <- wide[[a_col]]
  b <- wide[[b_col]]
  keep <- complete.cases(a, b)
  a <- a[keep]
  b <- b[keep]
  d <- b - a
  p <- if (length(d) >= 2) tryCatch(wilcox.test(b, a, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
  cd <- if (length(d) >= 2 && stats::sd(d, na.rm = TRUE) > 0) mean(d, na.rm = TRUE) / stats::sd(d, na.rm = TRUE) else NA_real_
  data.frame(
    n_pairs = length(d),
    mean_adjacent = mean(a, na.rm = TRUE),
    mean_core = mean(b, na.rm = TRUE),
    delta_core_minus_adjacent = mean(d, na.rm = TRUE),
    wilcoxon_p = p,
    paired_cohen_d = cd,
    stringsAsFactors = FALSE
  )
}

unpaired_stats <- function(df, value_col, group_col = "location", group_a = "adjacent", group_b = "core") {
  a <- df[[value_col]][df[[group_col]] == group_a]
  b <- df[[value_col]][df[[group_col]] == group_b]
  p <- if (length(a) >= 1 && length(b) >= 1) tryCatch(wilcox.test(b, a, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
  data.frame(
    n_adjacent = length(a),
    n_core = length(b),
    mean_adjacent = mean(a, na.rm = TRUE),
    mean_core = mean(b, na.rm = TRUE),
    delta_core_minus_adjacent = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE),
    wilcoxon_p = p,
    stringsAsFactors = FALSE
  )
}

theme_main <- theme_bw(base_family = "Arial", base_size = 11) +
  theme(
    panel.grid.major = element_line(color = light_gray, linewidth = 0.30),
    panel.grid.minor = element_line(color = light_gray, linewidth = 0.18),
    axis.title = element_text(color = "black", size = 13),
    axis.text = element_text(color = "black", size = 11),
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 11),
    plot.title = element_text(color = "black", face = "bold", size = 12),
    plot.subtitle = element_text(color = "black", size = 9),
    strip.text = element_text(color = "black", face = "bold", size = 11),
    strip.background = element_rect(fill = "#F4F4F4", color = mid_gray)
  )

pheno <- load_pheno("GSE159677")
sample_celltype <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_sample_celltype_summary.tsv"))
score_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_celltype_restricted_program_scores.tsv"))
macro_cont <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_cell_table.tsv"))

sample_meta <- pheno[, c("sample_id", "group_label", "title"), drop = FALSE]
sample_meta$patient_id <- sub(".*Patient ([0-9]+).*", "P\\1", sample_meta$title)
sample_meta$location <- sample_meta$group_label

## Panel A data
panel_a_sample <- merge(sample_celltype, sample_meta[, c("sample_id", "patient_id"), drop = FALSE], by = "sample_id", all.x = TRUE)
panel_a_sample$cell_type <- factor(panel_a_sample$cell_type, levels = celltype_levels)
panel_a_sample$cell_type_label <- factor(celltype_labels[as.character(panel_a_sample$cell_type)], levels = unname(celltype_labels[celltype_levels]))
panel_a_sample$location <- factor(panel_a_sample$location, levels = c("adjacent", "core"))

panel_a_mean <- aggregate(fraction ~ location + cell_type_label, data = panel_a_sample, FUN = mean)
panel_a_mean$n_samples <- aggregate(sample_id ~ location + cell_type_label, data = panel_a_sample, FUN = function(x) length(unique(x)))$sample_id

panel_a_stats <- do.call(rbind, lapply(split(panel_a_sample, panel_a_sample$cell_type_label), function(df) {
  out <- paired_stats(df, value_col = "fraction")
  out$cell_type_label <- as.character(df$cell_type_label[1])
  out
}))
panel_a_export <- merge(panel_a_sample, panel_a_stats, by = "cell_type_label", all.x = TRUE)
write.csv(panel_a_export, project_path("results", "figure3", "figure3_panelA_cell_composition.csv"), row.names = FALSE)

## Panel B data
panel_b_sample <- score_df[score_df$gene_set == "NPL_FOAM_MACROPHAGE_COMPACT", , drop = FALSE]
panel_b_sample <- merge(panel_b_sample, sample_meta[, c("sample_id", "patient_id"), drop = FALSE], by = "sample_id", all.x = TRUE)
panel_b_sample$cell_type <- factor(panel_b_sample$cell_type, levels = celltype_levels)
panel_b_sample$cell_type_label <- factor(celltype_labels[as.character(panel_b_sample$cell_type)], levels = unname(celltype_labels[celltype_levels]))
panel_b_sample$location <- factor(panel_b_sample$location, levels = c("adjacent", "core"))

panel_b_stats <- do.call(rbind, lapply(split(panel_b_sample, panel_b_sample$cell_type_label), function(df) {
  out <- paired_stats(df, value_col = "mean_score")
  out$cell_type_label <- as.character(df$cell_type_label[1])
  out
}))
panel_b_counts <- aggregate(sample_id ~ cell_type_label, data = panel_b_sample, FUN = function(x) length(unique(x)))
names(panel_b_counts)[2] <- "n_samples"
panel_b_counts$facet_label <- paste0(panel_b_counts$cell_type_label, "\n(n = ", panel_b_counts$n_samples, " samples)")
panel_b_label_map <- setNames(panel_b_counts$facet_label, panel_b_counts$cell_type_label)
panel_b_sample$facet_label <- factor(panel_b_label_map[as.character(panel_b_sample$cell_type_label)], levels = panel_b_counts$facet_label)
panel_b_export <- merge(panel_b_sample, panel_b_stats, by = "cell_type_label", all.x = TRUE)
write.csv(panel_b_export, project_path("results", "figure3", "figure3_panelB_celltype_module_scores.csv"), row.names = FALSE)

## Panel C data
panel_c_sample <- aggregate(
  npl_module_score ~ sample_id + sample_name + patient_id + location + dominant_state,
  data = macro_cont[macro_cont$dominant_state %in% state_levels & !is.na(macro_cont$npl_module_score), , drop = FALSE],
  FUN = mean
)
panel_c_counts <- aggregate(
  npl_module_score ~ sample_id + sample_name + patient_id + location + dominant_state,
  data = macro_cont[macro_cont$dominant_state %in% state_levels & !is.na(macro_cont$npl_module_score), , drop = FALSE],
  FUN = length
)
names(panel_c_counts)[names(panel_c_counts) == "npl_module_score"] <- "n_cells"
panel_c_sample <- merge(panel_c_sample, panel_c_counts, by = c("sample_id", "sample_name", "patient_id", "location", "dominant_state"))
panel_c_sample$state_label <- factor(state_labels[panel_c_sample$dominant_state], levels = unname(state_labels[state_levels]))
panel_c_sample$location <- factor(panel_c_sample$location, levels = c("adjacent", "core"))

panel_c_stats <- do.call(rbind, lapply(split(panel_c_sample, panel_c_sample$state_label), function(df) {
  out <- unpaired_stats(df, value_col = "npl_module_score")
  out$state_label <- as.character(df$state_label[1])
  out
}))
panel_c_n <- aggregate(sample_id ~ state_label, data = panel_c_sample, FUN = function(x) length(unique(x)))
names(panel_c_n)[2] <- "n_samples"
panel_c_n$facet_label <- paste0(panel_c_n$state_label, "\n(n = ", panel_c_n$n_samples, " samples)")
panel_c_label_map <- setNames(panel_c_n$facet_label, panel_c_n$state_label)
panel_c_sample$facet_label <- factor(panel_c_label_map[as.character(panel_c_sample$state_label)], levels = panel_c_n$facet_label)
panel_c_export <- merge(panel_c_sample, panel_c_stats, by = "state_label", all.x = TRUE)
write.csv(panel_c_export, project_path("results", "figure3", "figure3_panelC_macrophage_state_module_scores.csv"), row.names = FALSE)

## Panel D data
sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("^GSM", basename(sample_dirs))]
gene_rows <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  sample_id <- sub("_.*$", "", sample_name)
  counts <- read_10x_triplet(find_triplet_dir(sample_dir))
  norm <- normalize_log1p(counts)
  sample_macro <- macro_cont[macro_cont$sample_id == sample_id & macro_cont$dominant_state %in% state_levels, , drop = FALSE]
  if (nrow(sample_macro) == 0) next
  sample_macro$barcode <- sub("^.*\\|", "", sample_macro$cell_id)
  keep <- intersect(sample_macro$barcode, colnames(norm))
  genes_present <- intersect(compact_genes, rownames(norm))
  if (length(keep) == 0 || length(genes_present) == 0) next
  expr_sub <- norm[genes_present, keep, drop = FALSE]
  sample_macro <- sample_macro[match(keep, sample_macro$barcode), , drop = FALSE]
  for (gene in genes_present) {
    gene_rows[[length(gene_rows) + 1]] <- data.frame(
      sample_id = sample_id,
      sample_name = sample_name,
      patient_id = unique(sample_macro$patient_id)[1],
      location = unique(sample_macro$location)[1],
      dominant_state = sample_macro$dominant_state,
      state_label = state_labels[sample_macro$dominant_state],
      gene_symbol = gene,
      expr = as.numeric(expr_sub[gene, ]),
      expressed = as.integer(expr_sub[gene, ] > 0),
      stringsAsFactors = FALSE
    )
  }
}

gene_long <- do.call(rbind, gene_rows)
panel_d_sample <- aggregate(
  cbind(expr, expressed) ~ sample_id + patient_id + location + state_label + gene_symbol,
  data = gene_long,
  FUN = mean
)
names(panel_d_sample)[names(panel_d_sample) == "expressed"] <- "pct_expressing"
panel_d_summary <- aggregate(
  cbind(expr, pct_expressing) ~ state_label + gene_symbol,
  data = panel_d_sample,
  FUN = mean
)
panel_d_summary$n_samples <- aggregate(sample_id ~ state_label + gene_symbol, data = panel_d_sample, FUN = function(x) length(unique(x)))$sample_id
panel_d_summary$z_expr <- ave(panel_d_summary$expr, panel_d_summary$gene_symbol, FUN = function(x) as.numeric(scale(x)))
panel_d_summary$gene_symbol <- factor(panel_d_summary$gene_symbol, levels = compact_gene_order_plot)
panel_d_summary$state_label <- factor(panel_d_summary$state_label, levels = unname(state_labels[state_levels]))
write.csv(panel_d_summary, project_path("results", "figure3", "figure3_panelD_gene_dotplot_summary.csv"), row.names = FALSE)

## Key numbers
mac_fraction <- panel_a_stats[panel_a_stats$cell_type_label == "Macrophage", , drop = FALSE]
mac_module <- panel_b_stats[panel_b_stats$cell_type_label == "Macrophage", , drop = FALSE]
foam_vs_c1q <- panel_c_sample[panel_c_sample$state_label %in% c("FOAM/TREM2-like", "C1Q-like"), , drop = FALSE]
foam_mean <- mean(foam_vs_c1q$npl_module_score[foam_vs_c1q$state_label == "FOAM/TREM2-like"], na.rm = TRUE)
c1q_mean <- mean(foam_vs_c1q$npl_module_score[foam_vs_c1q$state_label == "C1Q-like"], na.rm = TRUE)
foam_p <- tryCatch(
  wilcox.test(
    foam_vs_c1q$npl_module_score[foam_vs_c1q$state_label == "FOAM/TREM2-like"],
    foam_vs_c1q$npl_module_score[foam_vs_c1q$state_label == "C1Q-like"],
    exact = FALSE
  )$p.value,
  error = function(e) NA_real_
)
foam_top <- panel_d_summary[panel_d_summary$state_label == "FOAM/TREM2-like", , drop = FALSE]
foam_top <- foam_top[order(-foam_top$z_expr, -foam_top$pct_expressing), ]
top_gene_txt <- paste(as.character(foam_top$gene_symbol[1:5]), collapse = ", ")

key_lines <- c(
  sprintf("Macrophage fraction: core mean %.3f versus adjacent mean %.3f (delta %.3f; paired Wilcoxon P = %s).",
          mac_fraction$mean_core, mac_fraction$mean_adjacent, mac_fraction$delta_core_minus_adjacent,
          formatC(mac_fraction$wilcoxon_p, format = "e", digits = 2)),
  sprintf("Compact NPL module score in macrophages: core mean %.3f versus adjacent mean %.3f (delta %.3f; paired Wilcoxon P = %s).",
          mac_module$mean_core, mac_module$mean_adjacent, mac_module$delta_core_minus_adjacent,
          formatC(mac_module$wilcoxon_p, format = "e", digits = 2)),
  sprintf("Compact NPL module score across macrophage states: FOAM/TREM2-like mean %.3f versus C1Q-like mean %.3f (descriptive Wilcoxon P = %s).",
          foam_mean, c1q_mean, formatC(foam_p, format = "e", digits = 2)),
  sprintf("Top compact genes enriched in FOAM/TREM2-like macrophages by gene-wise standardized expression: %s.", top_gene_txt)
)
writeLines(key_lines, con = project_path("results", "figure3", "Figure3_key_numbers.txt"))

## Plotting
panel_a_plot <- ggplot(panel_a_mean, aes(x = cell_type_label, y = fraction, fill = location)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "grey20", linewidth = 0.22) +
  geom_point(
    data = panel_a_sample,
    aes(x = cell_type_label, y = fraction, group = location),
    inherit.aes = FALSE,
    position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.72),
    size = 1.4,
    color = mid_gray,
    alpha = 0.85
  ) +
  scale_fill_manual(values = c(adjacent = adj_color, core = core_color), name = NULL) +
  labs(
    x = NULL,
    y = "Mean cell fraction"
  ) +
  theme_main +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

panel_b_plot <- ggplot(panel_b_sample, aes(x = location, y = mean_score, fill = location)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85, color = "grey20", linewidth = 0.22) +
  geom_point(
    aes(group = patient_id),
    position = position_jitter(width = 0.08, height = 0),
    size = 1.5,
    color = mid_gray,
    alpha = 0.9
  ) +
  scale_fill_manual(values = c(adjacent = adj_color, core = core_color), guide = "none") +
  facet_wrap(~facet_label, ncol = 4, scales = "free_y") +
  labs(
    x = NULL,
    y = "Per-sample mean module score"
  ) +
  theme_main +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

panel_c_plot <- ggplot(panel_c_sample, aes(x = location, y = npl_module_score, fill = location)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.88, color = "grey20", linewidth = 0.22) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 1.6,
    color = mid_gray,
    alpha = 0.9
  ) +
  scale_fill_manual(values = c(adjacent = adj_color, core = core_color), guide = "none") +
  facet_wrap(~facet_label, ncol = 2, scales = "free_y") +
  labs(
    x = NULL,
    y = "Per-sample mean module score"
  ) +
  theme_main +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

panel_d_plot <- ggplot(panel_d_summary, aes(x = state_label, y = gene_symbol, size = pct_expressing, color = z_expr)) +
  geom_point(alpha = 0.95) +
  scale_color_gradient2(low = light_gray, mid = mid_blue, high = core_color, midpoint = 0, name = "Mean expression\n(gene-wise z)") +
  scale_size_area(max_size = 8.5, name = "Fraction expressing") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_main +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  )

final_fig <- panel_a_plot / panel_b_plot / (panel_c_plot | panel_d_plot) +
  plot_layout(heights = c(0.90, 1.25, 1.15), widths = c(1, 1.08))

ggsave(
  filename = project_path("figures", "main", "Figure3_single_cell_localization_compact_NPL_module.png"),
  plot = final_fig,
  width = 16,
  height = 18,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = project_path("figures", "main", "Figure3_single_cell_localization_compact_NPL_module.pdf"),
  plot = final_fig,
  width = 16,
  height = 18,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

cat("Saved figure and panel tables to results/figure3 and figures/main\n")
