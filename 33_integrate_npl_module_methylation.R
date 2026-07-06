source("script/R/00_project_config.R")

ensure_project_dirs()

bulk_cons <- read_tsv_auto(project_path("res", "tables", "bulk", "athero_bulk_discovery_consensus.tsv.gz"))
meth <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE46394_methylation_limma.tsv.gz"))
modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
mechanism_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "athero_mechanism_gene_sets.tsv"))

meth <- meth[!(is.na(meth$gene_symbol) | meth$gene_symbol == ""), ]
meth <- meth[order(meth$adj.P.Val), ]
meth_gene <- meth[!duplicated(meth$gene_symbol), c("gene_symbol", "delta_beta", "adj.P.Val")]
names(meth_gene)[2:3] <- c("delta_beta_meth", "adj_p_meth")

bulk_sub <- bulk_cons[, c(
  "gene_symbol", "integrated_priority_score", "logFC_GSE43292", "logFC_GSE100927",
  "logFC_GSE28829", "discovery_fisher_fdr"
)]

integrated <- merge(bulk_sub, meth_gene, by = "gene_symbol", all.x = TRUE)
integrated$transcriptome_methylation_direction_match <- with(
  integrated,
  ifelse(is.na(delta_beta_meth), NA, sign(logFC_GSE43292 + logFC_GSE100927 + logFC_GSE28829) != sign(delta_beta_meth))
)

module_list <- split(modules$gene_symbol, modules$module_name)
module_list <- c(module_list, split(mechanism_sets$gene_symbol, mechanism_sets$gene_set))

module_rows <- list()
for (module_name in names(module_list)) {
  genes <- unique(module_list[[module_name]])
  sub <- integrated[integrated$gene_symbol %in% genes, , drop = FALSE]
  if (nrow(sub) == 0) {
    next
  }
  module_rows[[module_name]] <- data.frame(
    module_name = module_name,
    n_genes = length(genes),
    n_with_methylation = sum(!is.na(sub$delta_beta_meth)),
    mean_abs_bulk_logfc = mean(abs(c(sub$logFC_GSE43292, sub$logFC_GSE100927, sub$logFC_GSE28829)), na.rm = TRUE),
    mean_delta_beta = mean(sub$delta_beta_meth, na.rm = TRUE),
    median_delta_beta = median(sub$delta_beta_meth, na.rm = TRUE),
    n_direction_match = sum(sub$transcriptome_methylation_direction_match %in% TRUE, na.rm = TRUE),
    n_methylation_fdr = sum(sub$adj_p_meth < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
module_summary <- do.call(rbind, module_rows)

npl_detail <- integrated[integrated$gene_symbol %in% unique(module_list$NPL_FOAM_MACROPHAGE_COMPACT), , drop = FALSE]
npl_detail <- npl_detail[order(-rowMeans(abs(as.matrix(npl_detail[, c("logFC_GSE43292", "logFC_GSE100927", "logFC_GSE28829")])), na.rm = TRUE)), ]

write_tsv(integrated, project_path("res", "tables", "mechanism", "athero_transcriptome_methylation_integrated.tsv"))
write_tsv(module_summary, project_path("res", "tables", "mechanism", "athero_module_methylation_summary.tsv"))
write_tsv(npl_detail, project_path("res", "qc", "mechanism", "npl_module_methylation_gene_detail.tsv"))

cat("Transcriptome-methylation integration tables written to res/tables/mechanism\n")
print(module_summary[module_summary$module_name %in% c("NPL_FOAM_MACROPHAGE_COMPACT", "KEAP1_NRF2_RESPONSE", "ENDOTHELIAL_ACTIVATION_STRESS", "SMC_PHENOTYPE_SWITCH_STRESS"), ])
