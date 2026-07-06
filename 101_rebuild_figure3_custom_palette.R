source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

ensure_project_dirs()

gold_col <- "#F5BD4D"
blue_col <- "#005493"
accent_col <- "#C34062"
grey_col <- "#7A7A7A"

meta <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_cell_table.tsv"))
meta$location <- factor(meta$location, levels = c("adjacent", "core"))
meta$dominant_state <- factor(
  meta$dominant_state,
  levels = c("MACROPHAGE_C1Q", "MACROPHAGE_FOAM_TREM2", "MACROPHAGE_IFN", "MACROPHAGE_INFLAMMATORY")
)

state_labels <- c(
  MACROPHAGE_C1Q = "C1Q",
  MACROPHAGE_FOAM_TREM2 = "FOAM/TREM2",
  MACROPHAGE_IFN = "IFN",
  MACROPHAGE_INFLAMMATORY = "Inflammatory"
)

state_palette <- c(
  MACROPHAGE_C1Q = gold_col,
  MACROPHAGE_FOAM_TREM2 = blue_col,
  MACROPHAGE_IFN = grey_col,
  MACROPHAGE_INFLAMMATORY = accent_col
)

common_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 17, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )

p_state <- ggplot(meta, aes(pc_1, pc_2, color = dominant_state)) +
  geom_point(size = 0.45, alpha = 0.75) +
  scale_color_manual(values = state_palette, labels = state_labels) +
  common_theme +
  theme(legend.position = "right") +
  labs(
    title = "Macrophage continuum state map",
    subtitle = "GSE159677 plaque-core macrophage states",
    x = "PC1",
    y = "PC2",
    color = "State"
  )

p_pt <- ggplot(meta, aes(pc_1, pc_2, color = pseudotime)) +
  geom_point(size = 0.45, alpha = 0.8) +
  scale_color_gradientn(colours = c(gold_col, accent_col, blue_col)) +
  common_theme +
  theme(legend.position = "right") +
  labs(
    title = "Rooted pseudotime map",
    subtitle = "Higher values align with plaque-core foam remodeling",
    x = "PC1",
    y = "PC2",
    color = "Pseudotime"
  )

p_box <- ggplot(meta, aes(location, pseudotime, fill = location)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.55) +
  geom_boxplot(width = 0.16, outlier.size = 0.2, linewidth = 0.35, color = "grey20") +
  scale_fill_manual(values = c(adjacent = gold_col, core = blue_col)) +
  facet_wrap(~ dominant_state, scales = "free_y", ncol = 2, labeller = labeller(dominant_state = state_labels)) +
  common_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11)
  ) +
  labs(
    title = "Pseudotime by location and dominant state",
    x = NULL,
    y = "Pseudotime"
  )

fig_dir <- project_path("manuscript", "figures_final")
out <- file.path(fig_dir, "Figure3_macrophage_continuum.png")
png(out, width = 3200, height = 3200, res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 2, ncol = 2,
  heights = unit(c(1.0, 1.15), "null"),
  widths = unit(c(1, 1), "null")
)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(editGrob(ggplotGrob(p_state), vp = viewport(width = 0.92, height = 0.92)))
grid.text("A", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(editGrob(ggplotGrob(p_pt), vp = viewport(width = 0.92, height = 0.92)))
grid.text("B", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2))
grid.draw(editGrob(ggplotGrob(p_box), vp = viewport(width = 0.98, height = 0.96)))
grid.text("C", x = unit(0.01, "npc"), y = unit(0.99, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

dev.off()
cat("Wrote", out, "\n")

out_a <- file.path(fig_dir, "Figure4_panel_A_macrophage_continuum_state_map.png")
out_b <- file.path(fig_dir, "Figure4_panel_B_rooted_pseudotime_map.png")
out_a_square <- file.path(fig_dir, "Figure4_panel_A_macrophage_continuum_state_map_square.png")
out_b_square <- file.path(fig_dir, "Figure4_panel_B_rooted_pseudotime_map_square.png")

png(out_a, width = 2400, height = 2400, res = 300, bg = "white")
print(
  p_state +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13)
    )
)
dev.off()

png(out_b, width = 2400, height = 2400, res = 300, bg = "white")
print(
  p_pt +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13)
    )
)
dev.off()

cat("Wrote", out_a, "\n")
cat("Wrote", out_b, "\n")

png(out_a_square, width = 2200, height = 2200, res = 300, bg = "white")
print(
  p_state +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )
)
dev.off()

png(out_b_square, width = 2200, height = 2200, res = 300, bg = "white")
print(
  p_pt +
    theme(
      legend.position = "none",
      aspect.ratio = 1,
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14)
    )
)
dev.off()

cat("Wrote", out_a_square, "\n")
cat("Wrote", out_b_square, "\n")
