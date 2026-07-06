source("script/R/00_project_config.R")

ensure_project_dirs()

collapse_to_gene_level <- function(expr, mapping) {
  sample_ids <- setdiff(colnames(expr), "feature_id")
  merged <- merge(expr, mapping, by = "feature_id", all.x = TRUE, sort = FALSE)
  merged <- merged[!(is.na(merged$gene_symbol) | merged$gene_symbol == ""), , drop = FALSE]

  if (nrow(merged) == 0) {
    stop("No mapped features retained after joining expression matrix to platform mapping.")
  }

  expr_only <- merged[, sample_ids, drop = FALSE]
  row_variance <- apply(expr_only, 1, var, na.rm = TRUE)
  merged$row_variance <- row_variance

  split_index <- split(seq_len(nrow(merged)), merged$gene_symbol)
  keep_index <- vapply(split_index, function(idx) idx[which.max(merged$row_variance[idx])], integer(1))

  collapsed <- merged[keep_index, c("gene_symbol", "feature_id", sample_ids), drop = FALSE]
  collapsed <- collapsed[order(collapsed$gene_symbol), , drop = FALSE]
  rownames(collapsed) <- NULL
  collapsed
}

build_summary_row <- function(dataset_id, original_expr, gene_expr) {
  data.frame(
    dataset_id = dataset_id,
    n_original_features = nrow(original_expr),
    n_gene_level_rows = nrow(gene_expr),
    n_samples = ncol(gene_expr) - 2,
    stringsAsFactors = FALSE
  )
}

datasets <- list(
  list(dataset_id = "GSE43292", expr_file = "GSE43292_expr.tsv.gz", mapping_file = "GPL6244_mapping.tsv.gz"),
  list(dataset_id = "GSE100927", expr_file = "GSE100927_expr.tsv.gz", mapping_file = "GPL17077_mapping.tsv.gz"),
  list(dataset_id = "GSE28829", expr_file = "GSE28829_expr.tsv.gz", mapping_file = "GPL570_mapping.tsv.gz"),
  list(dataset_id = "GSE21545", expr_file = "GSE21545_expr.tsv.gz", mapping_file = "GPL570_mapping.tsv.gz")
)

summary_rows <- list()

for (dataset in datasets) {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk", dataset$expr_file))
  mapping <- read_tsv_auto(project_path("data", "processed", "annotation", dataset$mapping_file))
  gene_expr <- collapse_to_gene_level(expr, mapping)
  write_tsv_gz(gene_expr, project_path("data", "processed", "bulk_gene", paste0(dataset$dataset_id, "_gene_expr.tsv.gz")))
  summary_rows[[dataset$dataset_id]] <- build_summary_row(dataset$dataset_id, expr, gene_expr)
}

summary_table <- do.call(rbind, summary_rows)
write_tsv(summary_table, project_path("res", "qc", "bulk", "bulk_gene_level_summary.tsv"))

cat("Gene-level bulk matrices written to data/processed/bulk_gene\n")
print(summary_table)
