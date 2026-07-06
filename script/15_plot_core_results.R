source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(ggplot2))

ensure_project_dirs()

bulk_scores <- read_tsv_auto(project_path("res", "tables", "mechanism", "bulk_mechanism_scores.tsv.gz"))
sc_group <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_group_celltype_summary.tsv"))
sc_mech <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_group_mechanism_scores.tsv"))

p1 <- ggplot(bulk_scores, aes(x = group_std, y = score, fill = group_std)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
  facet_grid(dataset_id ~ gene_set, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Mean z-score", title = "Bulk mechanism scores across cohorts")

ggsave(
  filename = project_path("figure", "export", "mechanism", "bulk_mechanism_scores_boxplot.pdf"),
  plot = p1,
  width = 12,
  height = 7
)

p2 <- ggplot(sc_group, aes(x = cell_type, y = fraction, fill = location)) +
  geom_col(position = "dodge") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Mean fraction", title = "GSE159677 cell-type composition")

ggsave(
  filename = project_path("figure", "export", "mechanism", "gse159677_group_celltype_barplot.pdf"),
  plot = p2,
  width = 8,
  height = 5
)

p3 <- ggplot(sc_mech, aes(x = gene_set, y = mean_score, fill = location)) +
  geom_col(position = "dodge") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Mean score", title = "GSE159677 mechanism scores")

ggsave(
  filename = project_path("figure", "export", "mechanism", "gse159677_group_mechanism_barplot.pdf"),
  plot = p3,
  width = 8,
  height = 5
)

cat("Core result figures written to figure/export/mechanism\n")
