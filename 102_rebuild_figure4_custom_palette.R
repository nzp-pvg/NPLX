source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(dplyr)
})

ensure_project_dirs()

gold_col <- "#F5BD4D"
blue_col <- "#005493"
accent_col <- "#C34062"
grey_col <- "#8F8F8F"

patient_palette <- c(P1 = accent_col, P2 = gold_col, P3 = blue_col)

module_scores <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_module_scores.tsv"))
pair_deltas <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_pair_deltas.tsv"))
knk_diff <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_sctenifoldknk_core_macrophage_diff.tsv"))
knk_enrich <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_sctenifoldknk_core_macrophage_enrichment.tsv"))
state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))

module_scores$location <- factor(module_scores$location, levels = c("adjacent", "core"))
module_scores$patient_id <- factor(module_scores$patient_id, levels = c("P1", "P2", "P3"))
module_scores$module_name <- factor(
  module_scores$module_name,
  levels = c("NPL_FOAM_MACROPHAGE_COMPACT", "NPL_NEIGHBORHOOD_RAW")
)
module_labels <- c(
  NPL_FOAM_MACROPHAGE_COMPACT = "Compact module",
  NPL_NEIGHBORHOOD_RAW = "NPL neighborhood"
)

representative_genes <- c("NPL", "FABP5", "GPNMB", "APOC1", "SPP1", "APOE")
gene_df <- pair_deltas[pair_deltas$gene_symbol %in% representative_genes, , drop = FALSE]
gene_long <- rbind(
  data.frame(patient_id = gene_df$patient_id, gene_symbol = gene_df$gene_symbol, location = "adjacent", logcpm = gene_df$logcpm.adjacent),
  data.frame(patient_id = gene_df$patient_id, gene_symbol = gene_df$gene_symbol, location = "core", logcpm = gene_df$logcpm.core)
)
gene_long$location <- factor(gene_long$location, levels = c("adjacent", "core"))
gene_long$patient_id <- factor(gene_long$patient_id, levels = c("P1", "P2", "P3"))
gene_long$gene_symbol <- factor(gene_long$gene_symbol, levels = representative_genes)

knk_enrich_plot <- knk_enrich %>%
  dplyr::filter(
    gene_set %in% c("NPL_FOAM_MACROPHAGE_COMPACT", "NPL_NEIGHBORHOOD_RAW"),
    !grepl("_summary$", gene_set)
  ) %>%
  dplyr::mutate(
    gene_set = factor(
      gene_set,
      levels = c("NPL_FOAM_MACROPHAGE_COMPACT", "NPL_NEIGHBORHOOD_RAW"),
      labels = c("Compact module", "NPL neighborhood")
    ),
    top_n_num = as.numeric(top_n),
    top_n = factor(top_n, levels = c(25, 50, 100, 200)),
    enrichment_fraction = in_top / set_size
  )

total_ranked_genes <- nrow(knk_diff)
knk_enrich_plot$random_expectation <- knk_enrich_plot$top_n_num / total_ranked_genes
random_expectation_df <- unique(knk_enrich_plot[, c("top_n_num", "random_expectation")])

ranked_diff <- knk_diff
ranked_diff$perturb_rank <- seq_len(nrow(ranked_diff))
program_order <- c("MACROPHAGE_FOAM_TREM2", "MACROPHAGE_C1Q", "MACROPHAGE_INFLAMMATORY")
program_labels <- c(
  MACROPHAGE_FOAM_TREM2 = "FOAM/TREM2-like",
  MACROPHAGE_C1Q = "C1Q-like",
  MACROPHAGE_INFLAMMATORY = "Inflammatory"
)
program_palette <- c(
  "FOAM/TREM2-like" = blue_col,
  "C1Q-like" = gold_col,
  "Inflammatory" = accent_col
)

bootstrap_mean_ci <- function(x, n_boot = 2000, conf = 0.95) {
  if (length(x) == 0 || all(is.na(x))) {
    return(c(mean = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  x <- x[is.finite(x)]
  if (length(x) == 1) {
    return(c(mean = x, lower = x, upper = x))
  }
  boot_means <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  alpha <- (1 - conf) / 2
  c(
    mean = mean(x),
    lower = unname(quantile(boot_means, probs = alpha)),
    upper = unname(quantile(boot_means, probs = 1 - alpha))
  )
}

program_rank_df <- do.call(rbind, lapply(program_order, function(gs) {
  genes <- unique(state_sets$gene_symbol[state_sets$state_set == gs])
  sub <- ranked_diff[ranked_diff$gene %in% genes, , drop = FALSE]
  stats <- bootstrap_mean_ci(sub$perturb_rank)
  data.frame(
    program = program_labels[[gs]],
    n_genes = nrow(sub),
    mean_rank = stats[["mean"]],
    lower_ci = stats[["lower"]],
    upper_ci = stats[["upper"]],
    stringsAsFactors = FALSE
  )
}))
program_rank_df$program <- factor(program_rank_df$program, levels = rev(program_labels))

compact_top200 <- subset(knk_enrich_plot, gene_set == "Compact module" & top_n_num == 200)

common_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 15),
    plot.margin = margin(8, 8, 8, 8)
  )

p_module <- ggplot(module_scores, aes(location, module_score, group = patient_id, color = patient_id)) +
  geom_point(size = 2.9) +
  geom_line(linewidth = 1.0, alpha = 0.85) +
  scale_color_manual(values = patient_palette) +
  facet_wrap(~ module_name, scales = "free_y", ncol = 2, labeller = labeller(module_name = module_labels)) +
  common_theme +
  theme(legend.position = "right") +
  labs(
    x = NULL,
    y = "Mean logCPM",
    color = "Patient"
  )

p_gene <- ggplot(gene_long, aes(location, logcpm, group = patient_id, color = patient_id)) +
  geom_point(size = 2.0) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = patient_palette) +
  facet_wrap(~ gene_symbol, scales = "free_y", ncol = 3) +
  common_theme +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "logCPM"
  )

p_enrich <- ggplot(knk_enrich_plot, aes(top_n_num, enrichment_fraction, group = gene_set, color = gene_set)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.1) +
  geom_line(
    data = random_expectation_df,
    aes(x = top_n_num, y = random_expectation, group = 1),
    inherit.aes = FALSE,
    color = grey_col,
    linetype = "dashed",
    linewidth = 1.0
  ) +
  geom_point(
    data = random_expectation_df,
    aes(x = top_n_num, y = random_expectation),
    inherit.aes = FALSE,
    color = grey_col,
    size = 2.4
  ) +
  scale_x_continuous(breaks = c(25, 50, 100, 200)) +
  scale_color_manual(values = c("Compact module" = accent_col, "NPL neighborhood" = blue_col)) +
  common_theme +
  theme(legend.position = "right") +
  labs(
    x = "Top perturbed genes",
    y = "Fraction of gene set recovered",
    color = NULL
  ) +
  coord_cartesian(ylim = c(0, max(knk_enrich_plot$enrichment_fraction, na.rm = TRUE) + 0.1))

p_program <- ggplot(program_rank_df, aes(y = program, x = mean_rank, fill = program)) +
  geom_col(width = 0.62, color = NA, alpha = 0.95) +
  geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci), width = 0.18, linewidth = 0.9, color = "black") +
  scale_fill_manual(values = program_palette) +
  common_theme +
  theme(legend.position = "none") +
  labs(
    x = "Mean perturbation rank (lower = stronger effect)",
    y = NULL
  )

fig_dir <- project_path("manuscript", "figures_final")
out_old <- file.path(fig_dir, "Figure4_pseudobulk_virtualKO.png")
out_new <- file.path(fig_dir, "Figure5_NPL_functional_embedding.png")
out_panel_c_square <- file.path(fig_dir, "Figure5_panel_C_NPL_perturbation_square.png")

draw_figure <- function(out_file) {
  png(out_file, width = 3400, height = 4400, res = 300, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 4, ncol = 1,
    heights = unit(c(0.8, 1.15, 1.5, 1.0), "null")
  )))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(editGrob(ggplotGrob(p_module), vp = viewport(width = 0.98, height = 0.98)))
  grid.text("A", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(editGrob(ggplotGrob(p_gene), vp = viewport(width = 0.98, height = 0.98)))
  grid.text("B", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.draw(editGrob(ggplotGrob(p_enrich), vp = viewport(width = 0.98, height = 0.98)))
  grid.text("C", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
  popViewport()

  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
  grid.draw(editGrob(ggplotGrob(p_program), vp = viewport(width = 0.98, height = 0.98)))
  grid.text("D", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
  popViewport()

  dev.off()
}

draw_figure(out_old)
draw_figure(out_new)

png(out_panel_c_square, width = 2200, height = 2200, res = 300, bg = "white")
print(
  p_enrich +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16)
    )
)
dev.off()

cat("Wrote", out_old, "\n")
cat("Wrote", out_new, "\n")
cat("Wrote", out_panel_c_square, "\n")
