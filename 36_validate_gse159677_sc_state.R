source("script/R/00_project_config.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

ensure_project_dirs()
dir.create(project_path("results", "single_cell_state_validation"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figures"), recursive = TRUE, showWarnings = FALSE)

compact_genes <- c("NPL", "FABP5", "GPNMB", "APOC1", "PLA2G7", "SPP1", "CD36", "CYP27A1", "APOE")
compact_gene_order_plot <- c("NPL", "FABP5", "GPNMB", "APOC1", "APOE", "CD36", "CYP27A1", "PLA2G7", "SPP1")

celltype_levels <- c("B_CELL", "ENDOTHELIAL", "FIBROBLAST", "MACROPHAGE", "MAST", "SMC", "T_CELL", "UNRESOLVED")
celltype_labels <- c(
  B_CELL = "B cell",
  ENDOTHELIAL = "Endothelial",
  FIBROBLAST = "Fibroblast",
  MACROPHAGE = "Macrophage",
  MAST = "Mast",
  SMC = "SMC",
  T_CELL = "T cell",
  UNRESOLVED = "Unresolved"
)

state_levels <- c("MACROPHAGE_C1Q", "MACROPHAGE_FOAM_TREM2", "MACROPHAGE_INFLAMMATORY", "MACROPHAGE_IFN")
state_labels <- c(
  MACROPHAGE_C1Q = "C1Q-like",
  MACROPHAGE_FOAM_TREM2 = "FOAM/TREM2-like",
  MACROPHAGE_INFLAMMATORY = "Inflammatory",
  MACROPHAGE_IFN = "IFN-like"
)

core_color <- "#005493"
adj_color <- "#F5BD4D"
obs_color <- "#C34062"

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

pair_summary <- function(df, score_col = "mean_score", group_col = "location", pair_col = "patient_id", group_a = "adjacent", group_b = "core") {
  required <- c(pair_col, group_col, score_col)
  if (!all(required %in% colnames(df))) {
    return(data.frame(
      n_pairs = 0,
      n_samples_a = if (group_col %in% colnames(df)) sum(df[[group_col]] == group_a, na.rm = TRUE) else 0,
      n_samples_b = if (group_col %in% colnames(df)) sum(df[[group_col]] == group_b, na.rm = TRUE) else 0,
      mean_group_a = if (all(c(score_col, group_col) %in% colnames(df))) mean(df[[score_col]][df[[group_col]] == group_a], na.rm = TRUE) else NA_real_,
      mean_group_b = if (all(c(score_col, group_col) %in% colnames(df))) mean(df[[score_col]][df[[group_col]] == group_b], na.rm = TRUE) else NA_real_,
      median_group_a = if (all(c(score_col, group_col) %in% colnames(df))) median(df[[score_col]][df[[group_col]] == group_a], na.rm = TRUE) else NA_real_,
      median_group_b = if (all(c(score_col, group_col) %in% colnames(df))) median(df[[score_col]][df[[group_col]] == group_b], na.rm = TRUE) else NA_real_,
      delta_b_minus_a = NA_real_,
      median_delta_b_minus_a = NA_real_,
      paired_cohen_d = NA_real_,
      wilcoxon_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  wide <- reshape(df[, required, drop = FALSE], idvar = pair_col, timevar = group_col, direction = "wide")
  a_name <- paste0(score_col, ".", group_a)
  b_name <- paste0(score_col, ".", group_b)
  if (!(a_name %in% colnames(wide)) || !(b_name %in% colnames(wide))) {
    return(data.frame(
      n_pairs = 0,
      n_samples_a = sum(group_a %in% df[[group_col]]),
      n_samples_b = sum(group_b %in% df[[group_col]]),
      mean_group_a = mean(df[[score_col]][df[[group_col]] == group_a], na.rm = TRUE),
      mean_group_b = mean(df[[score_col]][df[[group_col]] == group_b], na.rm = TRUE),
      median_group_a = median(df[[score_col]][df[[group_col]] == group_a], na.rm = TRUE),
      median_group_b = median(df[[score_col]][df[[group_col]] == group_b], na.rm = TRUE),
      delta_b_minus_a = NA_real_,
      median_delta_b_minus_a = NA_real_,
      paired_cohen_d = NA_real_,
      wilcoxon_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  a <- wide[[a_name]]
  b <- wide[[b_name]]
  keep <- complete.cases(a, b)
  a <- a[keep]
  b <- b[keep]
  diff <- b - a
  pval <- if (length(diff) >= 2) {
    tryCatch(wilcox.test(b, a, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
  data.frame(
    n_pairs = length(diff),
    n_samples_a = sum(!is.na(a)),
    n_samples_b = sum(!is.na(b)),
    mean_group_a = mean(a, na.rm = TRUE),
    mean_group_b = mean(b, na.rm = TRUE),
    median_group_a = median(a, na.rm = TRUE),
    median_group_b = median(b, na.rm = TRUE),
    delta_b_minus_a = mean(diff, na.rm = TRUE),
    median_delta_b_minus_a = median(diff, na.rm = TRUE),
    paired_cohen_d = ifelse(stats::sd(diff, na.rm = TRUE) == 0, NA_real_, mean(diff, na.rm = TRUE) / stats::sd(diff, na.rm = TRUE)),
    wilcoxon_p = pval,
    stringsAsFactors = FALSE
  )
}

group_summary_unpaired <- function(df, score_col = "module_score", group_col = "location", group_a = "adjacent", group_b = "core") {
  if (!all(c(score_col, group_col) %in% colnames(df))) {
    return(data.frame(
      n_samples_a = 0,
      n_samples_b = 0,
      mean_group_a = NA_real_,
      mean_group_b = NA_real_,
      median_group_a = NA_real_,
      median_group_b = NA_real_,
      delta_b_minus_a = NA_real_,
      median_delta_b_minus_a = NA_real_,
      wilcoxon_p = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  a <- df[[score_col]][df[[group_col]] == group_a]
  b <- df[[score_col]][df[[group_col]] == group_b]
  pval <- if (length(a) >= 1 && length(b) >= 1) {
    tryCatch(wilcox.test(b, a, exact = FALSE)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
  data.frame(
    n_samples_a = length(a),
    n_samples_b = length(b),
    mean_group_a = mean(a, na.rm = TRUE),
    mean_group_b = mean(b, na.rm = TRUE),
    median_group_a = median(a, na.rm = TRUE),
    median_group_b = median(b, na.rm = TRUE),
    delta_b_minus_a = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE),
    median_delta_b_minus_a = median(b, na.rm = TRUE) - median(a, na.rm = TRUE),
    wilcoxon_p = pval,
    stringsAsFactors = FALSE
  )
}

pheno <- load_pheno("GSE159677")
cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
macro_cont <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_cell_table.tsv"))
macro_pseudo <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_module_scores.tsv"))

sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("^GSM", basename(sample_dirs))]

sample_meta <- pheno[, c("sample_id", "group_label", "title"), drop = FALSE]
sample_meta$patient_id <- sub(".*Patient ([0-9]+).*", "P\\1", sample_meta$title)
sample_meta$sample_name <- cell_typing$sample_name[match(sample_meta$sample_id, cell_typing$sample_id)]

cell_rows <- list()
macro_rows <- list()
gene_rows <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  sample_id <- sub("_.*$", "", sample_name)
  counts <- read_10x_triplet(find_triplet_dir(sample_dir))
  norm <- normalize_log1p(counts)
  common_genes <- intersect(compact_genes, rownames(norm))
  module_score <- score_gene_set(norm, common_genes)
  cell_barcodes <- intersect(colnames(norm), cell_typing$barcode[cell_typing$sample_name == sample_name])

  if (length(cell_barcodes) == 0) {
    next
  }
  module_score <- module_score[match(cell_barcodes, colnames(norm))]
  cell_meta <- cell_typing[cell_typing$sample_name == sample_name & cell_typing$barcode %in% cell_barcodes, , drop = FALSE]
  cell_meta <- cell_meta[match(cell_barcodes, cell_meta$barcode), , drop = FALSE]
  cell_meta$score <- as.numeric(module_score)
  cell_meta$patient_id <- sample_meta$patient_id[match(sample_id, sample_meta$sample_id)]
  cell_meta$location <- sample_meta$group_label[match(sample_id, sample_meta$sample_id)]
  cell_meta$sample_name <- sample_name
  cell_rows[[sample_id]] <- cell_meta[, c("sample_id", "sample_name", "patient_id", "location", "barcode", "predicted_cell_type", "score"), drop = FALSE]

  macro_meta <- macro_cont[macro_cont$sample_id == sample_id, , drop = FALSE]
  if (nrow(macro_meta) > 0) {
    macro_barcodes <- sub("^.*\\|", "", macro_meta$cell_id)
    keep <- intersect(macro_barcodes, colnames(norm))
    if (length(keep) > 0) {
      pos <- match(keep, colnames(norm))
      genes_present <- intersect(compact_genes, rownames(norm))
      expr_sub <- norm[genes_present, pos, drop = FALSE]
      macro_join <- macro_meta[match(keep, macro_barcodes), , drop = FALSE]
      macro_join$patient_id <- sample_meta$patient_id[match(sample_id, sample_meta$sample_id)]
      macro_join$sample_name <- sample_name
      macro_join$location <- sample_meta$group_label[match(sample_id, sample_meta$sample_id)]
      macro_join$barcode <- keep
      macro_join$module_score <- as.numeric(module_score[pos])
      for (gene in genes_present) {
        gene_rows[[length(gene_rows) + 1]] <- data.frame(
          sample_id = sample_id,
          sample_name = sample_name,
          patient_id = sample_meta$patient_id[match(sample_id, sample_meta$sample_id)],
          location = sample_meta$group_label[match(sample_id, sample_meta$sample_id)],
          barcode = keep,
          dominant_state = macro_join$dominant_state,
          gene_symbol = gene,
          expr = as.numeric(expr_sub[gene, ]),
          expressed = as.integer(expr_sub[gene, ] > 0),
          stringsAsFactors = FALSE
        )
      }
      macro_rows[[sample_id]] <- macro_join[, c("sample_id", "sample_name", "patient_id", "location", "barcode", "dominant_state", "module_score"), drop = FALSE]
    }
  }
}

cell_table <- do.call(rbind, cell_rows)
macro_table <- do.call(rbind, macro_rows)
gene_long <- do.call(rbind, gene_rows)

cell_table$cell_type <- celltype_labels[cell_table$predicted_cell_type]
cell_table$cell_type <- factor(cell_table$cell_type, levels = unname(celltype_labels[celltype_levels]))

## ----------------------------------------------------------------------
## Cell-type level sample-aware summaries
## ----------------------------------------------------------------------
cell_sample_means <- aggregate(
  score ~ sample_id + sample_name + patient_id + location + cell_type,
  data = cell_table[cell_table$predicted_cell_type %in% celltype_levels & !is.na(cell_table$score), , drop = FALSE],
  FUN = mean
)
cell_sample_counts <- aggregate(
  score ~ sample_id + sample_name + patient_id + location + cell_type,
  data = cell_table[cell_table$predicted_cell_type %in% celltype_levels & !is.na(cell_table$score), , drop = FALSE],
  FUN = length
)
names(cell_sample_counts)[names(cell_sample_counts) == "score"] <- "n_cells"
cell_sample <- merge(cell_sample_means, cell_sample_counts, by = c("sample_id", "sample_name", "patient_id", "location", "cell_type"))

cell_stats_rows <- list()
for (ct in levels(factor(cell_sample$cell_type))) {
  sub <- cell_sample[cell_sample$cell_type == ct, , drop = FALSE]
  if (nrow(sub) == 0) next
  ps <- pair_summary(sub, score_col = "score", group_col = "location", pair_col = "patient_id", group_a = "adjacent", group_b = "core")
  cell_stats_rows[[ct]] <- data.frame(
    cell_type = ct,
    group_a = "adjacent",
    group_b = "core",
    ps,
    stringsAsFactors = FALSE
  )
}
celltype_stats <- do.call(rbind, cell_stats_rows)
celltype_stats$cell_type <- factor(celltype_stats$cell_type, levels = unname(celltype_labels[celltype_levels]))
celltype_n <- aggregate(sample_id ~ cell_type, data = cell_sample, FUN = function(x) length(unique(x)))
names(celltype_n)[2] <- "n_samples"
celltype_n$cell_type_label <- sprintf(
  "%s\nn = %s samples",
  unname(celltype_labels[celltype_n$cell_type]),
  celltype_n$n_samples
)
celltype_label_map <- setNames(celltype_n$cell_type_label, celltype_n$cell_type)

## ----------------------------------------------------------------------
## Macrophage state level sample-aware summaries
## ----------------------------------------------------------------------
macro_state <- merge(
  macro_table,
  data.frame(cell_id = macro_cont$cell_id, stringsAsFactors = FALSE),
  by.x = "barcode",
  by.y = "cell_id",
  all.x = TRUE
)
macro_state <- macro_table
macro_state$state_label <- state_labels[macro_state$dominant_state]
macro_state$state_label <- factor(macro_state$state_label, levels = unname(state_labels[state_levels]))

macro_state_sample_means <- aggregate(
  module_score ~ sample_id + sample_name + patient_id + location + state_label,
  data = macro_state[!is.na(macro_state$module_score) & !is.na(macro_state$state_label), , drop = FALSE],
  FUN = mean
)
macro_state_sample_counts <- aggregate(
  module_score ~ sample_id + sample_name + patient_id + location + state_label,
  data = macro_state[!is.na(macro_state$module_score) & !is.na(macro_state$state_label), , drop = FALSE],
  FUN = length
)
names(macro_state_sample_counts)[names(macro_state_sample_counts) == "module_score"] <- "n_cells"
macro_state_sample <- merge(macro_state_sample_means, macro_state_sample_counts, by = c("sample_id", "sample_name", "patient_id", "location", "state_label"))

state_stats_rows <- list()
for (st in levels(factor(macro_state_sample$state_label))) {
  sub <- macro_state_sample[macro_state_sample$state_label == st, , drop = FALSE]
  if (nrow(sub) == 0) next
  ps <- group_summary_unpaired(sub, score_col = "module_score", group_col = "location", group_a = "adjacent", group_b = "core")
  state_stats_rows[[st]] <- data.frame(
    state_label = st,
    dominant_state = names(state_labels)[match(st, state_labels)],
    group_a = "adjacent",
    group_b = "core",
    ps,
    stringsAsFactors = FALSE
  )
}
state_stats <- do.call(rbind, state_stats_rows)
state_stats$state_label <- factor(state_stats$state_label, levels = unname(state_labels[state_levels]))
state_n <- aggregate(sample_id ~ state_label, data = macro_state_sample, FUN = function(x) length(unique(x)))
names(state_n)[2] <- "n_samples"
state_n$state_label_label <- sprintf(
  "%s\nn = %s samples",
  unname(state_labels[state_n$state_label]),
  state_n$n_samples
)
state_label_map <- setNames(state_n$state_label_label, state_n$state_label)

## ----------------------------------------------------------------------
## Macrophage pseudobulk summary
## ----------------------------------------------------------------------
compact_pseudo <- macro_pseudo[macro_pseudo$module_name == "NPL_FOAM_MACROPHAGE_COMPACT", , drop = FALSE]
compact_pseudo$sample_role <- compact_pseudo$location
pair_levels <- unique(compact_pseudo$patient_id)
pseudo_rows <- list()
for (pid in pair_levels) {
  sub <- compact_pseudo[compact_pseudo$patient_id == pid, , drop = FALSE]
  if (!all(c("core", "adjacent") %in% sub$location)) {
    next
  }
  core_score <- sub$module_score[sub$location == "core"][1]
  adj_score <- sub$module_score[sub$location == "adjacent"][1]
  pseudo_rows[[pid]] <- data.frame(
    patient_id = pid,
    core_score = core_score,
    adjacent_score = adj_score,
    delta_core_minus_adjacent = core_score - adj_score,
    stringsAsFactors = FALSE
  )
}
pseudo_stats <- do.call(rbind, pseudo_rows)
pseudo_stats$paired_wilcoxon_p <- if (nrow(pseudo_stats) >= 2) {
  tryCatch(wilcox.test(pseudo_stats$core_score, pseudo_stats$adjacent_score, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
} else {
  NA_real_
}
pseudo_stats$n_pairs <- nrow(pseudo_stats)
pseudo_stats$mean_delta_core_minus_adjacent <- mean(pseudo_stats$delta_core_minus_adjacent, na.rm = TRUE)
pseudo_stats$median_delta_core_minus_adjacent <- median(pseudo_stats$delta_core_minus_adjacent, na.rm = TRUE)
pseudo_stats$paired_cohen_d <- ifelse(
  stats::sd(pseudo_stats$delta_core_minus_adjacent, na.rm = TRUE) == 0,
  NA_real_,
  mean(pseudo_stats$delta_core_minus_adjacent, na.rm = TRUE) / stats::sd(pseudo_stats$delta_core_minus_adjacent, na.rm = TRUE)
)

## ----------------------------------------------------------------------
## Macrophage state gene expression dot plot data
## ----------------------------------------------------------------------
gene_summary <- aggregate(
  cbind(expr, expressed) ~ sample_id + sample_name + patient_id + location + dominant_state + gene_symbol,
  data = gene_long[gene_long$dominant_state %in% state_levels & gene_long$gene_symbol %in% compact_genes, , drop = FALSE],
  FUN = mean
)
names(gene_summary)[names(gene_summary) == "expressed"] <- "pct_expressing"
gene_summary$state_label <- state_labels[gene_summary$dominant_state]
gene_summary$state_label <- factor(gene_summary$state_label, levels = unname(state_labels[state_levels]))

gene_state_summary <- aggregate(
  cbind(expr, pct_expressing) ~ state_label + gene_symbol,
  data = gene_summary,
  FUN = mean
)
gene_state_summary$n_samples <- aggregate(
  sample_id ~ state_label + gene_symbol,
  data = gene_summary,
  FUN = function(x) length(unique(x))
)$sample_id
gene_state_summary$z_expr <- ave(gene_state_summary$expr, gene_state_summary$gene_symbol, FUN = function(x) as.numeric(scale(x)))
gene_state_summary$gene_symbol <- factor(gene_state_summary$gene_symbol, levels = compact_gene_order_plot)
gene_state_summary$state_label <- factor(gene_state_summary$state_label, levels = unname(state_labels[state_levels]))
gene_n <- aggregate(sample_id ~ state_label + gene_symbol, data = gene_summary, FUN = function(x) length(unique(x)))
names(gene_n)[3] <- "n_samples"
gene_state_summary <- merge(gene_state_summary, gene_n, by = c("state_label", "gene_symbol"), all.x = TRUE)

## ----------------------------------------------------------------------
## Save tables
## ----------------------------------------------------------------------
write.csv(celltype_stats, project_path("results", "single_cell_state_validation", "celltype_module_score_stats.csv"), row.names = FALSE, quote = TRUE)
write.csv(state_stats, project_path("results", "single_cell_state_validation", "macrophage_state_module_score_stats.csv"), row.names = FALSE, quote = TRUE)
write.csv(pseudo_stats, project_path("results", "single_cell_state_validation", "pseudobulk_core_adjacent_stats.csv"), row.names = FALSE, quote = TRUE)
write.csv(gene_state_summary, project_path("results", "single_cell_state_validation", "macrophage_state_gene_expression_summary.csv"), row.names = FALSE, quote = TRUE)

## ----------------------------------------------------------------------
## Plot A: cell-type-level module score violin
## ----------------------------------------------------------------------
plot_a_df <- cell_sample
plot_a_df$cell_type <- factor(plot_a_df$cell_type, levels = unname(celltype_labels[celltype_levels]))
plot_a_df$cell_type_plot <- celltype_label_map[as.character(plot_a_df$cell_type)]
plot_a_df$location <- factor(plot_a_df$location, levels = c("adjacent", "core"))

plot_a <- ggplot(plot_a_df, aes(x = location, y = score, fill = location)) +
  geom_violin(scale = "width", trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.08, size = 1.1, alpha = 0.8, color = "#2c2c2c") +
  facet_wrap(~cell_type_plot, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(adjacent = adj_color, core = core_color), name = "Location") +
  labs(
    title = "A  Compact NPL module score across major vascular cell types",
    x = NULL,
    y = "Per-sample mean module score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.margin = margin(8, 10, 4, 8)
  )

## ----------------------------------------------------------------------
## Plot B: macrophage-state-level module score violin
## ----------------------------------------------------------------------
plot_b_df <- macro_state_sample
plot_b_df$state_label <- factor(plot_b_df$state_label, levels = unname(state_labels[state_levels]))
plot_b_df$state_label_plot <- state_label_map[as.character(plot_b_df$state_label)]
plot_b_df$location <- factor(plot_b_df$location, levels = c("adjacent", "core"))

plot_b <- ggplot(plot_b_df, aes(x = location, y = module_score, fill = location)) +
  geom_violin(scale = "width", trim = FALSE, color = NA, alpha = 0.8) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(width = 0.08, size = 1.1, alpha = 0.8, color = "#2c2c2c") +
  facet_wrap(~state_label_plot, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(adjacent = adj_color, core = core_color), name = "Location") +
  labs(
    title = "B  Compact NPL module score across macrophage states",
    x = NULL,
    y = "Per-sample mean module score"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.margin = margin(8, 10, 4, 8)
  )

## ----------------------------------------------------------------------
## Plot C: macrophage pseudobulk core vs adjacent
## ----------------------------------------------------------------------
plot_c <- ggplot(compact_pseudo, aes(x = location, y = module_score, group = patient_id, color = patient_id)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.3) +
  scale_x_discrete(limits = c("adjacent", "core")) +
  labs(
    title = "C  Paired macrophage pseudobulk compact module score",
    x = NULL,
    y = "Pseudobulk module score",
    color = "Patient"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(8, 10, 4, 8)
  )

## ----------------------------------------------------------------------
## Plot D: compact module gene expression across macrophage states
## ----------------------------------------------------------------------
plot_d <- ggplot(gene_state_summary, aes(x = state_label, y = gene_symbol, color = z_expr, size = pct_expressing)) +
  geom_point(alpha = 0.95) +
  scale_color_gradient2(
    low = adj_color,
    mid = "white",
    high = core_color,
    midpoint = 0,
    name = "Mean expression\n(z-score)"
  ) +
  scale_size_area(max_size = 9, name = "Mean expressing\ncells (%)") +
  labs(
    title = "D  Compact-module gene expression across macrophage states",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right",
    plot.margin = margin(8, 10, 4, 8)
  )

supp_fig <- (plot_a | plot_b) / plot_c / plot_d + plot_layout(heights = c(1.15, 0.75, 1.0))

ggsave(
  filename = project_path("figures", "SuppFig_sc_state_validation.png"),
  plot = supp_fig,
  width = 16,
  height = 18,
  dpi = 300
)

cat("Single-cell state validation complete.\n")
print(celltype_stats)
print(state_stats)
print(pseudo_stats)
