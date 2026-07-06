source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

ensure_project_dirs()

adjacent_col <- "#F5BD4D"
core_col <- "#005493"
accent_col <- "#C34062"

celltype_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_group_celltype_summary.tsv"))
score_df <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_celltype_restricted_program_scores.tsv"))

celltype_df$location <- factor(celltype_df$location, levels = c("adjacent", "core"))
celltype_df$cell_type <- factor(
  celltype_df$cell_type,
  levels = c("B_CELL", "ENDOTHELIAL", "FIBROBLAST", "MACROPHAGE", "MAST", "SMC", "T_CELL", "UNRESOLVED")
)

cell_type_labels <- c(
  B_CELL = "B cell",
  ENDOTHELIAL = "Endothelial",
  FIBROBLAST = "Fibroblast",
  MACROPHAGE = "Macrophage",
  MAST = "Mast",
  SMC = "SMC",
  T_CELL = "T cell",
  UNRESOLVED = "Unresolved"
)

focus_sets <- c(
  "ENDOTHELIAL_ACTIVATION_STRESS",
  "KEAP1_NRF2_RESPONSE",
  "NPL_FOAM_MACROPHAGE_COMPACT",
  "SMC_PHENOTYPE_SWITCH_STRESS"
)

score_df <- score_df[score_df$gene_set %in% focus_sets, , drop = FALSE]
score_df$location <- factor(score_df$location, levels = c("adjacent", "core"))
score_df$cell_type <- factor(
  score_df$cell_type,
  levels = c("B_CELL", "ENDOTHELIAL", "FIBROBLAST", "MACROPHAGE", "MAST", "SMC", "T_CELL", "UNRESOLVED")
)
score_df$gene_set <- factor(score_df$gene_set, levels = focus_sets)

gene_set_labels <- c(
  ENDOTHELIAL_ACTIVATION_STRESS = "Endothelial\nstress",
  KEAP1_NRF2_RESPONSE = "KEAP1/NRF2",
  NPL_FOAM_MACROPHAGE_COMPACT = "NPL foam\nmacrophage",
  SMC_PHENOTYPE_SWITCH_STRESS = "SMC phenotype\nswitch"
)

common_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 11)
  )

p1 <- ggplot(celltype_df, aes(x = cell_type, y = fraction, fill = location)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, color = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = c(adjacent = adjacent_col, core = core_col)) +
  common_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    legend.background = element_blank()
  ) +
  scale_x_discrete(labels = cell_type_labels) +
  labs(
    x = NULL,
    y = "Mean fraction",
    title = "GSE159677 cell-type composition",
    fill = "location"
  )

p2 <- ggplot(score_df, aes(x = location, y = mean_score, fill = location)) +
  geom_boxplot(
    width = 0.72,
    outlier.size = 0.5,
    linewidth = 0.35,
    color = "grey20"
  ) +
  scale_fill_manual(values = c(adjacent = adjacent_col, core = core_col)) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text.y = element_text(angle = 270, color = "white", face = "bold", size = 9),
    strip.text.x = element_text(color = "white", face = "bold", size = 9),
    panel.spacing = unit(0.08, "lines")
  ) +
  scale_x_discrete(labels = c(adjacent = "adjacent", core = "core")) +
  facet_grid(gene_set ~ cell_type, scales = "free_y", labeller = labeller(
    cell_type = cell_type_labels,
    gene_set = gene_set_labels
  )) +
  labs(
    title = "Cell type-restricted program scores in GSE159677",
    x = NULL,
    y = "Mean per-cell score"
  )

fig_dir <- project_path("manuscript", "figures_final")
png(file.path(fig_dir, "Figure2_singlecell_localization.png"), width = 2400, height = 3200, res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(0.95, 1.45), "null"))))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(editGrob(ggplotGrob(p1), vp = viewport(width = 0.96, height = 0.94)))
grid.text("A", x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(editGrob(ggplotGrob(p2), vp = viewport(width = 0.98, height = 0.96)))
grid.text("B", x = unit(0.01, "npc"), y = unit(0.99, "npc"), just = c("left", "top"),
          gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

dev.off()

cat("Wrote", file.path(fig_dir, "Figure2_singlecell_localization.png"), "\n")
