source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(jpeg)
})

ensure_project_dirs()

gold_col <- "#F5BD4D"
blue_col <- "#005493"
accent_col <- "#C34062"

score_df <- read.csv("/Users/chnqwe/Science/CVD/CVD_MS_2/wet_lab/qPCR/20260605/20260605_qpcr_scores.csv", check.names = FALSE)
norm_df <- read.csv("/Users/chnqwe/Science/CVD/CVD_MS_2/wet_lab/qPCR/20260605/20260605_qpcr_normalized.csv", check.names = FALSE)
if_path <- "/Users/chnqwe/Science/CVD/CVD_MS_2/wet_lab/IF/data/1/thumbnail.jpg"

group_levels <- c("Healthy", "AMI", "CAD")
group_palette <- c(Healthy = gold_col, AMI = accent_col, CAD = blue_col)

score_df$group <- factor(score_df$group, levels = group_levels)
score_df$score_label <- factor(
  score_df$score,
  levels = c("NPL_FOAM_MACROPHAGE_qPCR_SCORE_8G", "ENDOTHELIAL_ACTIVATION_PARTIAL_2G"),
  labels = c("NPL foam score", "Endothelial score")
)

heat_genes <- c("NPL", "FABP5", "GPNMB", "APOC1", "APOE", "PLA2G7", "SPP1", "CD36", "CCL2", "VWF")
sample_order <- c("107", "133", "137", "139", "145", "150", "241", "118", "119+")
norm_df$sample <- as.character(norm_df$sample)
norm_df$group <- factor(norm_df$group, levels = group_levels)
heat_df <- norm_df[norm_df$gene %in% heat_genes, c("sample", "group", "gene", "z_neg_deltaCt")]
heat_df$sample <- factor(heat_df$sample, levels = sample_order)
heat_df$gene <- factor(heat_df$gene, levels = rev(heat_genes))

common_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 15, face = "bold", color = "black"),
    strip.background = element_rect(fill = "black", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )

p_scores <- ggplot(score_df, aes(group, score_mean_z, fill = group, color = group)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.45, color = "grey40") +
  geom_boxplot(alpha = 0.18, outlier.shape = NA, linewidth = 0.8) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.6) +
  scale_fill_manual(values = group_palette) +
  scale_color_manual(values = group_palette) +
  facet_wrap(~ score_label, ncol = 2, scales = "free_y") +
  common_theme +
  theme(legend.position = "none") +
  labs(
    title = "Exploratory peripheral blood qPCR scores",
    x = NULL,
    y = "Mean z score"
  )

p_heat <- ggplot(heat_df, aes(sample, gene, fill = z_neg_deltaCt)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = blue_col, mid = "white", high = accent_col, midpoint = 0) +
  common_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Gene-level qPCR pattern",
    x = "Sample",
    y = NULL,
    fill = "z(-ΔCt)"
  )

if_img <- jpeg::readJPEG(if_path)
roi_box <- c(x1 = 110, y1 = 215, x2 = 235, y2 = 375)
enhance_if <- function(x, mult = 2.2, gamma = 0.9) {
  x2 <- pmin(x * mult, 1)
  x2 ^ gamma
}
if_overview <- enhance_if(if_img)
if_overview[roi_box["y1"]:roi_box["y2"], roi_box["x1"], ] <- 1
if_overview[roi_box["y1"]:roi_box["y2"], roi_box["x2"], ] <- 1
if_overview[roi_box["y1"], roi_box["x1"]:roi_box["x2"], ] <- 1
if_overview[roi_box["y2"], roi_box["x1"]:roi_box["x2"], ] <- 1
if_roi <- enhance_if(if_img[roi_box["y1"]:roi_box["y2"], roi_box["x1"]:roi_box["x2"], , drop = FALSE], mult = 3.0, gamma = 0.8)

fig_dir <- project_path("manuscript", "figures_final")
out <- file.path(fig_dir, "Figure5_qPCR_IF.png")
png(out, width = 3200, height = 3200, res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 2, ncol = 2,
  heights = unit(c(1.0, 1.0), "null"),
  widths = unit(c(1.0, 1.0), "null")
)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(editGrob(ggplotGrob(p_scores), vp = viewport(width = 0.98, height = 0.98)))
grid.text("A", x = unit(-0.03, "npc"), y = unit(1.02, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(editGrob(ggplotGrob(p_heat), vp = viewport(width = 0.98, height = 0.98)))
grid.text("B", x = unit(-0.03, "npc"), y = unit(1.02, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.text("C", x = unit(-0.03, "npc"), y = unit(1.02, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
grid.text("Whole-slide IF context", x = unit(0.06, "npc"), y = unit(0.94, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold", col = "black"))
pushViewport(viewport(x = 0.5, y = 0.44, width = 0.86, height = 0.74))
grid.raster(if_overview, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
popViewport()
grid.text("White box marks the plaque-rich ROI shown in panel D", x = unit(0.06, "npc"), y = unit(0.07, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 14, col = "black"))
popViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.text("D", x = unit(-0.03, "npc"), y = unit(1.02, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 28, fontface = "bold", col = "black"))
grid.text("Plaque-rich ROI", x = unit(0.06, "npc"), y = unit(0.94, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold", col = "black"))
pushViewport(viewport(x = 0.5, y = 0.44, width = 0.8, height = 0.72))
grid.raster(if_roi, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE)
popViewport()
grid.text("Merged NPL/CD68/GPNMB signal in a plaque-rich region", x = unit(0.06, "npc"), y = unit(0.07, "npc"),
          just = c("left", "top"), gp = gpar(fontsize = 14, col = "black"))
popViewport()

dev.off()
cat("Wrote", out, "\n")
