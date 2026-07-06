source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

ensure_project_dirs()

gold_col <- "#F5BD4D"
blue_col <- "#005493"
accent_col <- "#C34062"

contrast_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_module_bulk_contrast_summary.tsv"))
native_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_compact_native_eigengene_summary.tsv"))
projected_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_compact_projected_eigengene_summary.tsv"))
bulk_scores_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_module_bulk_scores.tsv"))
random_df <- read.csv(project_path("results", "cross_cohort_robustness", "random_module_benchmark.csv"), stringsAsFactors = FALSE)
random_summary_df <- read.csv(project_path("results", "cross_cohort_robustness", "random_module_benchmark_summary.csv"), stringsAsFactors = FALSE)

dataset_order <- c("GSE100927", "GSE28829", "GSE43292")
dataset_labels <- c(
  GSE100927 = "GSE100927",
  GSE28829 = "GSE28829",
  GSE43292 = "GSE43292"
)
dataset_fill <- c(
  GSE100927 = accent_col,
  GSE28829 = gold_col,
  GSE43292 = blue_col
)

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

contrast_df <- contrast_df[contrast_df$module_name == "NPL_FOAM_MACROPHAGE_COMPACT", , drop = FALSE]
contrast_df$dataset_id <- factor(contrast_df$dataset_id, levels = dataset_order)
contrast_df$label_y <- contrast_df$delta_b_minus_a + 0.06
contrast_df$p_label <- sprintf("P=%.1e", contrast_df$wilcox_p)

native_df$dataset_id <- factor(native_df$dataset_id, levels = dataset_order)
projected_df$dataset_id <- factor(projected_df$dataset_id, levels = dataset_order)
bulk_scores_df$dataset_id <- factor(bulk_scores_df$dataset_id, levels = dataset_order)
bulk_scores_df <- bulk_scores_df[bulk_scores_df$module_name == "NPL_FOAM_MACROPHAGE_COMPACT", , drop = FALSE]

line_df <- rbind(
  data.frame(dataset_id = native_df$dataset_id, type = "Within-cohort eigengene", cohen_d = native_df$cohen_d),
  data.frame(dataset_id = projected_df$dataset_id, type = "Transferred eigengene", cohen_d = projected_df$cohen_d)
)
line_df$type <- factor(line_df$type, levels = c("Within-cohort eigengene", "Transferred eigengene"))

cohort_effect_df <- do.call(
  rbind,
  lapply(split(bulk_scores_df, bulk_scores_df$dataset_id), function(df) {
    row <- contrast_df[contrast_df$dataset_id == as.character(df$dataset_id[1]), , drop = FALSE][1, ]
    a <- df$module_score[df$group_label == row$comparison_a]
    b <- df$module_score[df$group_label == row$comparison_b]
    data.frame(
      dataset_id = as.character(df$dataset_id[1]),
      effect_size = cohen_d(a, b),
      p_value = row$wilcox_p,
      stringsAsFactors = FALSE
    )
  })
)
cohort_effect_df$dataset_id <- factor(cohort_effect_df$dataset_id, levels = dataset_order)
cohort_effect_df$p_label <- ifelse(
  cohort_effect_df$p_value < 1e-4,
  "P < 1×10^-4",
  sprintf("P = %.1e", cohort_effect_df$p_value)
)

random_df$cohort_id <- factor(random_df$cohort_id, levels = dataset_order)
random_summary_df$cohort_id <- factor(random_summary_df$cohort_id, levels = dataset_order)
random_summary_df$label <- sprintf(
  "Observed Cohen's d = %.2f\nEmpirical percentile = %.3f\nEmpirical P = %.3f",
  random_summary_df$observed_effect_d,
  random_summary_df$empirical_percentile,
  random_summary_df$empirical_p
)

common_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 15, face = "bold", color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

p1 <- ggplot(cohort_effect_df, aes(x = dataset_id, y = effect_size, fill = dataset_id)) +
  geom_col(width = 0.62, color = "grey20", linewidth = 0.25) +
  geom_text(aes(y = effect_size + 0.06, label = p_label), size = 4.2, color = "black") +
  scale_fill_manual(values = dataset_fill) +
  scale_x_discrete(labels = dataset_labels) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11)
  ) +
  labs(
    x = NULL,
    y = "Compact module effect size (Cohen's d)",
    title = "Consistent compact NPL module elevation\nacross independent vascular cohorts"
  )

p2 <- ggplot(line_df, aes(x = dataset_id, y = cohen_d, color = type, group = type)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +
  geom_text(aes(label = sprintf("%.2f", cohen_d)), nudge_y = 0.03, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Within-cohort eigengene" = gold_col, "Transferred eigengene" = blue_col)) +
  scale_x_discrete(labels = dataset_labels) +
  common_theme +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 8)
  ) +
  labs(
    x = NULL,
    y = "Cohen's d",
    title = "Stable eigengene transfer",
    color = NULL
  )

p3 <- ggplot(random_df, aes(x = null_effect_d)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45, fill = "#7f8c8d", color = "white", alpha = 0.78) +
  geom_vline(
    data = random_summary_df,
    aes(xintercept = observed_effect_d),
    color = accent_col,
    linewidth = 0.9
  ) +
  geom_text(
    data = random_summary_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.05, vjust = 1.1,
    size = 4.0,
    inherit.aes = FALSE
  ) +
  facet_wrap(~cohort_id, nrow = 1, scales = "free_y") +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold")
  ) +
  labs(
    x = "Random module effect size (Cohen's d)",
    y = "Density",
    title = "Random 9-gene module benchmark"
  ) +
  coord_cartesian(clip = "off")

fig_dir <- project_path("manuscript", "figures_final")
png(file.path(fig_dir, "Figure1_bulk_module_transfer.png"), width = 3600, height = 1500, res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1.22, 0.88), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(editGrob(ggplotGrob(p1), vp = viewport(width = 0.98, height = 0.98)))
grid.text("A", x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(editGrob(ggplotGrob(p2), vp = viewport(width = 0.98, height = 0.98)))
grid.text("B", x = unit(0.06, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

dev.off()

png(file.path(fig_dir, "Figure1_bulk_module_transfer_check.png"), width = 3600, height = 2600, res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2, heights = unit(c(1, 1.05), "null"), widths = unit(c(1.22, 0.88), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(editGrob(ggplotGrob(p1), vp = viewport(width = 0.98, height = 0.98)))
grid.text("A", x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(editGrob(ggplotGrob(p2), vp = viewport(width = 0.98, height = 0.98)))
grid.text("B", x = unit(0.06, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
grid.draw(editGrob(ggplotGrob(p3), vp = viewport(width = 0.99, height = 0.97)))
grid.text("C", x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

dev.off()

cat("Wrote", file.path(fig_dir, "Figure1_bulk_module_transfer.png"), "\n")
cat("Wrote", file.path(fig_dir, "Figure1_bulk_module_transfer_check.png"), "\n")
