source("script/R/00_project_config.R")

ensure_project_dirs()

gene_sets <- list(
  OXIDATIVE_STRESS_CORE = c("NFE2L2", "KEAP1", "HMOX1", "NQO1", "GCLC", "GCLM", "TXN", "TXNRD1", "PRDX1", "PRDX2", "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "SRXN1"),
  KEAP1_NRF2_RESPONSE = c("NFE2L2", "KEAP1", "HMOX1", "NQO1", "GCLC", "GCLM", "TXNRD1", "SRXN1", "FTL", "FTH1", "GSR", "ME1", "IDH1", "AKR1C1", "AKR1B10"),
  MITOCHONDRIAL_REDOX_HOMEOSTASIS = c("SOD2", "PRDX3", "PRDX5", "TXN2", "TXNRD2", "GPX4", "CAT", "ETFB", "NDUFA4", "NDUFS1", "COX5A", "ATP5F1A", "MFN2", "PINK1", "PARK7"),
  ENDOTHELIAL_ACTIVATION_STRESS = c("SELE", "SELP", "VCAM1", "ICAM1", "KLF2", "KLF4", "NOS3", "EDN1", "THBD", "PECAM1", "VWF", "IL6", "CXCL8", "CCL2", "HMOX1"),
  SMC_PHENOTYPE_SWITCH_STRESS = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYLK", "COL1A1", "COL3A1", "SPP1", "FN1", "VIM", "LGALS3", "MMP2", "MMP9", "KLF4", "PDGFRB")
)

rows <- do.call(
  rbind,
  lapply(names(gene_sets), function(set_name) {
    data.frame(
      gene_set = set_name,
      gene_symbol = gene_sets[[set_name]],
      stringsAsFactors = FALSE
    )
  })
)

write_tsv(rows, project_path("res", "tables", "mechanism", "athero_mechanism_gene_sets.tsv"))
cat("Mechanism gene sets written to res/tables/mechanism/athero_mechanism_gene_sets.tsv\n")
