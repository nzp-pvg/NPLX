source("script/R/00_project_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

ensure_project_dirs()
dir.create(project_path("results", "cross_cohort_robustness"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("figures"), recursive = TRUE, showWarnings = FALSE)

set.seed(20260705)

compact_genes <- c(
  "NPL", "FABP5", "GPNMB", "APOC1", "PLA2G7",
  "SPP1", "CD36", "CYP27A1", "APOE"
)

dataset_map <- list(
  GSE43292 = list(group_a = "intact", group_b = "plaque", display = "GSE43292"),
  GSE100927 = list(group_a = "control", group_b = "atherosclerotic", display = "GSE100927"),
  GSE28829 = list(group_a = "early", group_b = "advanced", display = "GSE28829")
)

read_bulk_dataset <- function(dataset_id) {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", paste0(dataset_id, "_gene_expr.tsv.gz")))
  pheno <- load_pheno(dataset_id)
  rownames(expr) <- expr$gene_symbol
  expr_mat <- as.matrix(expr[, !(colnames(expr) %in% c("gene_symbol", "feature_id")), drop = FALSE])
  storage.mode(expr_mat) <- "double"
  list(expr = expr_mat, pheno = pheno)
}

rowwise_zscore <- function(expr_mat) {
  z <- t(scale(t(expr_mat)))
  z[is.na(z)] <- 0
  z
}

score_gene_set <- function(z_mat, genes) {
  genes <- intersect(genes, rownames(z_mat))
  if (length(genes) == 0) {
    stop("No overlapping genes for module scoring.")
  }
  as.numeric(colMeans(z_mat[genes, , drop = FALSE], na.rm = TRUE))
}

write_csv_like <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE, na = "", quote = TRUE)
}

cohen_d <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  pooled_sd <- sqrt((stats::var(x) + stats::var(y)) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  (mean(y) - mean(x)) / pooled_sd
}

cliffs_delta <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) {
    return(NA_real_)
  }
  gt <- outer(y, x, ">")
  lt <- outer(y, x, "<")
  (sum(gt) - sum(lt)) / (length(x) * length(y))
}

wilcox_p <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 1 || length(y) < 1) {
    return(NA_real_)
  }
  tryCatch(wilcox.test(y, x, exact = FALSE)$p.value, error = function(e) NA_real_)
}

extract_group_scores <- function(scores, pheno, group_a, group_b) {
  df <- data.frame(
    sample_id = names(scores),
    score = as.numeric(scores),
    stringsAsFactors = FALSE
  )
  keep_cols <- intersect(c("sample_id", "group_label", "group_std"), colnames(pheno))
  df <- merge(df, pheno[, keep_cols, drop = FALSE], by = "sample_id", all.x = TRUE)
  df <- df[df$group_label %in% c(group_a, group_b), , drop = FALSE]

  a <- df$score[df$group_label == group_a]
  b <- df$score[df$group_label == group_b]
  data.frame(
    group_a = group_a,
    group_b = group_b,
    n_group_a = length(a),
    n_group_b = length(b),
    mean_group_a = mean(a, na.rm = TRUE),
    mean_group_b = mean(b, na.rm = TRUE),
    median_group_a = stats::median(a, na.rm = TRUE),
    median_group_b = stats::median(b, na.rm = TRUE),
    delta_b_minus_a = mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE),
    cohen_d = cohen_d(a, b),
    cliffs_delta = cliffs_delta(a, b),
    wilcoxon_p = wilcox_p(a, b),
    stringsAsFactors = FALSE
  )
}

prepare_dataset <- function(dataset_id, group_a, group_b) {
  obj <- read_bulk_dataset(dataset_id)
  common <- intersect(compact_genes, rownames(obj$expr))
  if (length(common) < length(compact_genes)) {
    warning(dataset_id, ": missing compact genes -> ", paste(setdiff(compact_genes, common), collapse = ", "))
  }
  z_mat <- rowwise_zscore(obj$expr)
  score <- score_gene_set(z_mat, common)
  names(score) <- colnames(obj$expr)

  a_scores <- score[obj$pheno$sample_id[obj$pheno$group_label == group_a]]
  b_scores <- score[obj$pheno$sample_id[obj$pheno$group_label == group_b]]
  obs_d <- cohen_d(a_scores, b_scores)

  list(
    expr = obj$expr,
    pheno = obj$pheno,
    z_mat = z_mat,
    score = score,
    common_genes = common,
    observed_d = obs_d,
    group_a = group_a,
    group_b = group_b,
    cohort = dataset_id
  )
}

datasets <- lapply(names(dataset_map), function(id) {
  prepare_dataset(id, dataset_map[[id]]$group_a, dataset_map[[id]]$group_b)
})
names(datasets) <- names(dataset_map)

## ----------------------------------------------------------------------
## Observed module score statistics
## ----------------------------------------------------------------------
score_stats <- do.call(
  rbind,
  lapply(names(datasets), function(dataset_id) {
    obj <- datasets[[dataset_id]]
    stats_row <- extract_group_scores(
      scores = obj$score,
      pheno = obj$pheno,
      group_a = obj$group_a,
      group_b = obj$group_b
    )
    data.frame(
      cohort_id = dataset_id,
      score_name = "NPL_FOAM_MACROPHAGE_qPCR_SCORE",
      score_method = "within_cohort_gene_zscore_mean",
      n_genes = length(obj$common_genes),
      gene_set = paste(obj$common_genes, collapse = ";"),
      stats_row,
      stringsAsFactors = FALSE
    )
  })
)

## ----------------------------------------------------------------------
## Leave-one-gene-out stability
## ----------------------------------------------------------------------
loo_rows <- list()
for (dataset_id in names(datasets)) {
  obj <- datasets[[dataset_id]]
  full_score <- obj$score
  for (drop_gene in obj$common_genes) {
    loo_genes <- setdiff(obj$common_genes, drop_gene)
    loo_score <- score_gene_set(obj$z_mat, loo_genes)
    names(loo_score) <- names(full_score)
    cor_val <- suppressWarnings(cor(full_score, loo_score, method = "spearman", use = "pairwise.complete.obs"))
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      cohort_id = dataset_id,
      dropped_gene = drop_gene,
      loo_genes = paste(loo_genes, collapse = ";"),
      score_cor_spearman = cor_val,
      score_cor_loss = 1 - cor_val,
      stringsAsFactors = FALSE
    )
  }
}
loo_table <- do.call(rbind, loo_rows)

## ----------------------------------------------------------------------
## Random module benchmark
## ----------------------------------------------------------------------
valid_pools <- Reduce(
  intersect,
  lapply(datasets, function(obj) {
    g <- rownames(obj$z_mat)
    sd_vec <- apply(obj$z_mat, 1, stats::sd)
    g[is.finite(sd_vec) & sd_vec > 0]
  })
)
random_pool <- setdiff(valid_pools, compact_genes)
module_size <- length(compact_genes)
n_iter <- 5000

random_rows <- list()
benchmark_summary_rows <- list()
for (dataset_id in names(datasets)) {
  obj <- datasets[[dataset_id]]
  group_a_idx <- obj$pheno$sample_id[obj$pheno$group_label == obj$group_a]
  group_b_idx <- obj$pheno$sample_id[obj$pheno$group_label == obj$group_b]
  group_a_idx <- intersect(group_a_idx, names(obj$score))
  group_b_idx <- intersect(group_b_idx, names(obj$score))
  obs_d <- obj$observed_d
  null_effects <- numeric(n_iter)
  null_gene_sets <- character(n_iter)

  for (iter in seq_len(n_iter)) {
    genes <- sample(random_pool, size = module_size, replace = FALSE)
    null_gene_sets[iter] <- paste(genes, collapse = ";")
    null_score <- score_gene_set(obj$z_mat, genes)
    names(null_score) <- names(obj$score)
    null_effects[iter] <- cohen_d(null_score[group_a_idx], null_score[group_b_idx])
  }

  direction <- if (is.finite(obs_d) && obs_d >= 0) "upper" else "lower"
  if (direction == "upper") {
    empirical_p <- (sum(null_effects >= obs_d, na.rm = TRUE) + 1) / (sum(is.finite(null_effects)) + 1)
    empirical_percentile <- (sum(null_effects <= obs_d, na.rm = TRUE) + 1) / (sum(is.finite(null_effects)) + 1)
  } else {
    empirical_p <- (sum(null_effects <= obs_d, na.rm = TRUE) + 1) / (sum(is.finite(null_effects)) + 1)
    empirical_percentile <- (sum(null_effects >= obs_d, na.rm = TRUE) + 1) / (sum(is.finite(null_effects)) + 1)
  }

  random_rows[[dataset_id]] <- data.frame(
    cohort_id = dataset_id,
    iter = seq_len(n_iter),
    random_genes = null_gene_sets,
    null_effect_d = null_effects,
    observed_effect_d = obs_d,
    empirical_percentile = empirical_percentile,
    empirical_p = empirical_p,
    n_iter = n_iter,
    module_size = module_size,
    stringsAsFactors = FALSE
  )

  benchmark_summary_rows[[dataset_id]] <- data.frame(
    cohort_id = dataset_id,
    observed_effect_d = obs_d,
    null_mean_effect_d = mean(null_effects, na.rm = TRUE),
    null_sd_effect_d = stats::sd(null_effects, na.rm = TRUE),
    empirical_percentile = empirical_percentile,
    empirical_p = empirical_p,
    stringsAsFactors = FALSE
  )
}

random_table <- do.call(rbind, random_rows)
benchmark_summary <- do.call(rbind, benchmark_summary_rows)

## ----------------------------------------------------------------------
## Save tables
## ----------------------------------------------------------------------
write_csv_like(score_stats, project_path("results", "cross_cohort_robustness", "module_score_stats.csv"))
write_csv_like(loo_table, project_path("results", "cross_cohort_robustness", "leave_one_gene_out_stability.csv"))
write_csv_like(random_table, project_path("results", "cross_cohort_robustness", "random_module_benchmark.csv"))

write_csv_like(benchmark_summary, project_path("results", "cross_cohort_robustness", "random_module_benchmark_summary.csv"))

## ----------------------------------------------------------------------
## Supplementary figure
## ----------------------------------------------------------------------
plot_a_df <- do.call(
  rbind,
  lapply(names(datasets), function(dataset_id) {
    obj <- datasets[[dataset_id]]
    df <- data.frame(
      cohort_id = dataset_id,
      sample_id = names(obj$score),
      score = as.numeric(obj$score),
      group_label = obj$pheno$group_label[match(names(obj$score), obj$pheno$sample_id)],
      stringsAsFactors = FALSE
    )
    df$group_role <- ifelse(df$group_label == obj$group_a, "reference", "disease")
    df$group_label <- factor(df$group_label, levels = c(obj$group_a, obj$group_b))
    df
  })
)

plot_a_df$cohort_id <- factor(plot_a_df$cohort_id, levels = names(dataset_map))

plot_a <- ggplot(plot_a_df, aes(x = group_label, y = score, fill = group_role)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.88, color = "#333333") +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.8, color = "#2f2f2f") +
  facet_wrap(~cohort_id, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(reference = "#F5BD4D", disease = "#005493")) +
  guides(fill = guide_legend(title = NULL, override.aes = list(alpha = 0.88))) +
  labs(
    title = "A  Compact NPL module score across discovery and validation cohorts",
    x = NULL,
    y = "Compact module score (mean within-cohort z-score)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(8, 12, 4, 8)
  )

heat_df <- loo_table
heat_df$cohort_id <- factor(heat_df$cohort_id, levels = rev(names(dataset_map)))
heat_df$dropped_gene <- factor(heat_df$dropped_gene, levels = compact_genes)
heat_df$label <- sprintf("%.4f", heat_df$score_cor_spearman)

plot_b <- ggplot(heat_df, aes(x = dropped_gene, y = cohort_id, fill = score_cor_spearman)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.2, color = "black") +
  scale_fill_gradient(
    low = "#F5BD4D",
    high = "#005493",
    limits = range(heat_df$score_cor_spearman, na.rm = TRUE),
    name = "Spearman\ncorrelation"
  ) +
  labs(
    title = "B  Leave-one-gene-out stability",
    x = "Dropped compact-module gene",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(8, 12, 4, 8)
  )

random_plot_df <- random_table
random_plot_df$cohort_id <- factor(random_plot_df$cohort_id, levels = names(dataset_map))
summary_text <- benchmark_summary
summary_text$cohort_id <- factor(summary_text$cohort_id, levels = names(dataset_map))
summary_text$label <- sprintf(
  "obs d = %.2f\nemp. pct = %.3f\nemp. p = %.3f",
  summary_text$observed_effect_d,
  summary_text$empirical_percentile,
  summary_text$empirical_p
)

plot_c <- ggplot(random_plot_df, aes(x = null_effect_d)) +
  geom_histogram(aes(y = after_stat(density), fill = "Random null"), bins = 45, color = "white", alpha = 0.75) +
  geom_vline(
    data = benchmark_summary,
    aes(xintercept = observed_effect_d, color = "Observed effect"),
    linewidth = 0.9
  ) +
  geom_text(
    data = summary_text,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.05, vjust = 1.1,
    size = 3.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~factor(cohort_id, levels = names(dataset_map)), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("Random null" = "#7f8c8d"), name = NULL) +
  scale_color_manual(values = c("Observed effect" = "#C34062"), name = NULL) +
  labs(
    title = "C  Random 9-gene module benchmark",
    x = "Random module effect size (Cohen's d)",
    y = "Density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(8, 12, 4, 8)
  ) +
  coord_cartesian(clip = "off")

supp_fig <- plot_a / plot_b / plot_c + plot_layout(heights = c(1.25, 0.95, 1.0))

ggsave(
  filename = project_path("figures", "SuppFig_cross_cohort_robustness.png"),
  plot = supp_fig,
  width = 16,
  height = 14,
  dpi = 300
)

cat("Cross-cohort robustness analysis complete.\n")
print(score_stats)
print(benchmark_summary)
