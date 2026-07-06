source("script/R/00_project_config.R")

ensure_project_dirs()

gse43292 <- read_tsv_auto(project_path("res", "tables", "bulk", "GSE43292_limma.tsv.gz"))
gse100927 <- read_tsv_auto(project_path("res", "tables", "bulk", "GSE100927_limma.tsv.gz"))
gse28829 <- read_tsv_auto(project_path("res", "tables", "bulk", "GSE28829_limma.tsv.gz"))

g1 <- gse43292[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "sign_direction")]
names(g1)[-1] <- paste0(names(g1)[-1], "_GSE43292")
g2 <- gse100927[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "sign_direction")]
names(g2)[-1] <- paste0(names(g2)[-1], "_GSE100927")
g3 <- gse28829[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "sign_direction")]
names(g3)[-1] <- paste0(names(g3)[-1], "_GSE28829")

consensus <- merge(g1, g2, by = "gene_symbol", all = FALSE)
consensus$same_direction_discovery <- sign(consensus$logFC_GSE43292) == sign(consensus$logFC_GSE100927)
consensus$discovery_fisher_stat <- -2 * (log(consensus$P.Value_GSE43292) + log(consensus$P.Value_GSE100927))
consensus$discovery_fisher_p <- pchisq(consensus$discovery_fisher_stat, df = 4, lower.tail = FALSE)
consensus$discovery_fisher_fdr <- p.adjust(consensus$discovery_fisher_p, method = "BH")
consensus$mean_abs_logFC <- rowMeans(cbind(abs(consensus$logFC_GSE43292), abs(consensus$logFC_GSE100927)), na.rm = TRUE)
consensus <- merge(consensus, g3, by = "gene_symbol", all.x = TRUE)
consensus$progression_same_direction <- sign(consensus$logFC_GSE43292) == sign(consensus$logFC_GSE28829)
consensus$progression_support_nominal <- !is.na(consensus$P.Value_GSE28829) & consensus$P.Value_GSE28829 < 0.05
consensus$integrated_priority_score <- with(
  consensus,
  mean_abs_logFC *
    (same_direction_discovery * 2 + (discovery_fisher_fdr < 0.05) * 2 + progression_same_direction + progression_support_nominal)
)

consensus <- consensus[order(-consensus$integrated_priority_score, consensus$discovery_fisher_fdr), ]
write_tsv_gz(consensus, project_path("res", "tables", "bulk", "athero_bulk_discovery_consensus.tsv.gz"))

summary_table <- data.frame(
  n_overlap = nrow(consensus),
  n_same_direction_discovery = sum(consensus$same_direction_discovery, na.rm = TRUE),
  n_same_direction_discovery_fdr = sum(consensus$same_direction_discovery & consensus$discovery_fisher_fdr < 0.05, na.rm = TRUE),
  n_progression_same_direction = sum(consensus$progression_same_direction, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, project_path("res", "qc", "bulk", "athero_bulk_consensus_summary.tsv"))

cat("Bulk consensus written to res/tables/bulk/athero_bulk_discovery_consensus.tsv.gz\n")
print(summary_table)
