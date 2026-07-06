source("script/R/00_project_config.R")

ensure_project_dirs()

cons <- read_tsv_auto(project_path("res", "tables", "bulk", "athero_bulk_discovery_consensus.tsv.gz"))
meth <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE46394_methylation_limma.tsv.gz"))
gene_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "athero_mechanism_gene_sets.tsv"))

meth_gene <- meth[!(is.na(meth$gene_symbol) | meth$gene_symbol == ""), , drop = FALSE]
meth_gene <- meth_gene[order(meth_gene$adj.P.Val, -abs(meth_gene$delta_beta)), ]
meth_gene <- meth_gene[!duplicated(meth_gene$gene_symbol), c("gene_symbol", "feature_id", "delta_beta", "adj.P.Val", "sig_fdr", "sig_fdr_beta")]
names(meth_gene)[-1] <- paste0(names(meth_gene)[-1], "_meth")

candidate <- merge(cons, meth_gene, by = "gene_symbol", all.x = TRUE)
candidate$mechanism_member <- candidate$gene_symbol %in% unique(gene_sets$gene_symbol)
candidate$methylation_support <- ifelse(!is.na(candidate$adj.P.Val_meth) & candidate$adj.P.Val_meth < 0.05, 1, 0)
candidate$methylation_support_strict <- ifelse(!is.na(candidate$adj.P.Val_meth) & candidate$adj.P.Val_meth < 0.05 & abs(candidate$delta_beta_meth) >= 0.10, 1, 0)
candidate$bulk_support <- ifelse(candidate$same_direction_discovery, 1, 0) +
  ifelse(candidate$discovery_fisher_fdr < 0.05, 1, 0) +
  ifelse(candidate$progression_same_direction, 1, 0) +
  ifelse(candidate$progression_support_nominal, 1, 0)
candidate$support_points <- candidate$bulk_support +
  candidate$methylation_support +
  candidate$methylation_support_strict +
  ifelse(candidate$mechanism_member, 1, 0)
candidate$priority_rank <- rank(-candidate$integrated_priority_score, ties.method = "min")

candidate <- candidate[order(-candidate$support_points, candidate$priority_rank, candidate$discovery_fisher_fdr), ]
write_tsv_gz(candidate, project_path("res", "tables", "mechanism", "athero_multiomics_candidate_shortlist.tsv.gz"))

summary_table <- head(
  candidate[, c(
    "gene_symbol", "support_points", "integrated_priority_score", "mechanism_member",
    "adj.P.Val_meth", "delta_beta_meth", "logFC_GSE43292", "logFC_GSE100927", "logFC_GSE28829"
  )],
  20
)
write_tsv(summary_table, project_path("res", "qc", "mechanism", "athero_multiomics_candidate_top20.tsv"))

cat("Atherosclerosis multi-omics shortlist written to res/tables/mechanism/athero_multiomics_candidate_shortlist.tsv.gz\n")
print(summary_table)
