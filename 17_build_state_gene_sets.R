source("script/R/00_project_config.R")

ensure_project_dirs()

state_sets <- list(
  MACROPHAGE_INFLAMMATORY = c("IL1B", "CXCL2", "CXCL3", "CCL3", "CCL4", "NFKBIA", "TNF", "PTGS2", "FOS", "JUNB"),
  MACROPHAGE_FOAM_TREM2 = c("TREM2", "APOE", "LPL", "SPP1", "GPNMB", "LGALS3", "FABP5", "CTSB", "CTSD", "LIPA"),
  MACROPHAGE_C1Q = c("C1QA", "C1QB", "C1QC", "FCER1G", "TYROBP", "AIF1", "LST1", "SAT1", "FCGR3A"),
  MACROPHAGE_IFN = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFITM3", "MX1", "CXCL10", "STAT1"),
  SMC_CONTRACTILE = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "CALD1", "DES"),
  SMC_FIBROMYOCYTE = c("COL1A1", "COL1A2", "COL3A1", "FN1", "VCAN", "LUM", "DCN", "C7", "FBLN1"),
  SMC_OSTEO_STRESS = c("SPP1", "IBSP", "BMP2", "RUNX2", "SOX9", "ALPL", "COL2A1", "CTSK"),
  SMC_INFLAMMATORY = c("VCAM1", "CCL2", "CXCL8", "IL6", "HMOX1", "KLF4", "JUN", "FOS")
)

rows <- do.call(
  rbind,
  lapply(names(state_sets), function(state_name) {
    data.frame(state_set = state_name, gene_symbol = state_sets[[state_name]], stringsAsFactors = FALSE)
  })
)

write_tsv(rows, project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))
cat("State gene sets written to res/tables/mechanism/vascular_state_gene_sets.tsv\n")
