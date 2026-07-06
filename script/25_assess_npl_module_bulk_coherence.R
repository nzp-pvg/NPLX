source("script/R/00_project_config.R")

ensure_project_dirs()

modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
compact_genes <- modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"]

score_dataset_coherence <- function(dataset_id, group_a, group_b) {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", paste0(dataset_id, "_gene_expr.tsv.gz")))
  pheno <- load_pheno(dataset_id)
  rownames(expr) <- expr$gene_symbol
  expr_mat <- as.matrix(expr[, !(colnames(expr) %in% c("gene_symbol", "feature_id")), drop = FALSE])
  storage.mode(expr_mat) <- "double"

  genes <- intersect(compact_genes, rownames(expr_mat))
  sub_expr <- expr_mat[genes, , drop = FALSE]
  rows <- list()
  edge_rows <- list()

  for (group_label in c(group_a, group_b)) {
    keep <- pheno$sample_id[pheno$group_label == group_label]
    keep <- intersect(colnames(sub_expr), keep)
    group_expr <- sub_expr[, keep, drop = FALSE]
    cor_mat <- cor(t(group_expr), method = "pearson")
    cor_vals <- cor_mat[upper.tri(cor_mat)]
    rows[[length(rows) + 1]] <- data.frame(
      dataset_id = dataset_id,
      group_label = group_label,
      n_samples = length(keep),
      n_genes = nrow(group_expr),
      mean_pairwise_cor = mean(cor_vals, na.rm = TRUE),
      median_pairwise_cor = median(cor_vals, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    edge_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
    edge_rows[[length(edge_rows) + 1]] <- data.frame(
      dataset_id = dataset_id,
      group_label = group_label,
      gene_a = rownames(cor_mat)[edge_idx[, 1]],
      gene_b = colnames(cor_mat)[edge_idx[, 2]],
      pearson_r = cor_mat[edge_idx],
      stringsAsFactors = FALSE
    )
  }

  list(summary = do.call(rbind, rows), edges = do.call(rbind, edge_rows))
}

res_43292 <- score_dataset_coherence("GSE43292", "intact", "plaque")
res_100927 <- score_dataset_coherence("GSE100927", "control", "atherosclerotic")
res_28829 <- score_dataset_coherence("GSE28829", "early", "advanced")

summary_table <- rbind(res_43292$summary, res_100927$summary, res_28829$summary)
edge_table <- rbind(res_43292$edges, res_100927$edges, res_28829$edges)

npl_edges <- edge_table[edge_table$gene_a == "NPL" | edge_table$gene_b == "NPL", ]
npl_edges <- npl_edges[order(npl_edges$dataset_id, npl_edges$group_label, -npl_edges$pearson_r), ]

write_tsv(summary_table, project_path("res", "tables", "mechanism", "npl_module_bulk_coherence_summary.tsv"))
write_tsv(edge_table, project_path("res", "tables", "mechanism", "npl_module_bulk_edge_correlations.tsv"))
write_tsv(npl_edges, project_path("res", "qc", "mechanism", "npl_module_bulk_npl_edges.tsv"))

cat("NPL module bulk coherence tables written to res/tables/mechanism and res/qc/mechanism\n")
print(summary_table)
cat("\nTop NPL edges:\n")
print(head(npl_edges, 30))
