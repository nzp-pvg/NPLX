source("script/R/00_project_config.R")

ensure_project_dirs()

marker_sets <- list(
  ENDOTHELIAL = c("PECAM1", "VWF", "KDR", "CLDN5", "ESAM", "EMCN"),
  SMC = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "CALD1"),
  MACROPHAGE = c("LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CTSB"),
  FIBROBLAST = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"),
  T_CELL = c("CD3D", "CD3E", "IL7R", "LTB", "TRBC1", "TRBC2"),
  B_CELL = c("MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA"),
  MAST = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2")
)

rows <- do.call(
  rbind,
  lapply(names(marker_sets), function(celltype) {
    data.frame(cell_type = celltype, gene_symbol = marker_sets[[celltype]], stringsAsFactors = FALSE)
  })
)

write_tsv(rows, project_path("res", "tables", "mechanism", "gse159677_marker_sets.tsv"))
cat("Single-cell marker sets written to res/tables/mechanism/gse159677_marker_sets.tsv\n")
