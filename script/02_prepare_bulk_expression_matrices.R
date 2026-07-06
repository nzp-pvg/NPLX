source("script/R/00_project_config.R")

ensure_project_dirs()

reorder_matrix <- function(expr, sample_ids) {
  missing_ids <- setdiff(sample_ids, colnames(expr))
  if (length(missing_ids) > 0) {
    stop("Missing samples in matrix: ", paste(missing_ids, collapse = ", "))
  }
  expr[, c("feature_id", sample_ids), drop = FALSE]
}

build_summary_row <- function(dataset_id, expr, pheno, comparison_label) {
  data.frame(
    dataset_id = dataset_id,
    comparison_label = comparison_label,
    n_features = nrow(expr),
    n_samples = ncol(expr) - 1,
    n_groups = length(unique(pheno$group_std)),
    groups = paste(sort(unique(pheno$group_std)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

bulk_dir <- project_path("data", "processed", "bulk")
summary_rows <- list()

gse43292_expr <- read_geo_series_matrix(project_path("data", "raw", "series_matrix", "GSE43292_series_matrix.txt.gz"))
gse43292_pheno <- load_pheno("GSE43292")
gse43292_pheno$group_std <- ifelse(gse43292_pheno$group_label == "plaque", "Plaque", "Intact")
gse43292_pheno$keep_for_primary_contrast <- "yes"
gse43292_expr <- reorder_matrix(gse43292_expr, gse43292_pheno$sample_id)
write_tsv_gz(gse43292_expr, file.path(bulk_dir, "GSE43292_expr.tsv.gz"))
write_tsv(gse43292_pheno, file.path(bulk_dir, "GSE43292_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE43292", gse43292_expr, gse43292_pheno, "plaque_vs_intact_paired")

gse100927_expr <- read_geo_series_matrix(project_path("data", "raw", "series_matrix", "GSE100927_series_matrix.txt.gz"))
gse100927_pheno <- load_pheno("GSE100927")
gse100927_pheno$group_std <- ifelse(gse100927_pheno$group_label == "atherosclerotic", "Atherosclerotic", "Control")
gse100927_pheno$vascular_bed <- gse100927_pheno$subgroup
gse100927_pheno$keep_for_primary_contrast <- ifelse(gse100927_pheno$group_std %in% c("Atherosclerotic", "Control"), "yes", "no")
gse100927_expr <- reorder_matrix(gse100927_expr, gse100927_pheno$sample_id)
write_tsv_gz(gse100927_expr, file.path(bulk_dir, "GSE100927_expr.tsv.gz"))
write_tsv(gse100927_pheno, file.path(bulk_dir, "GSE100927_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE100927", gse100927_expr, gse100927_pheno, "atherosclerotic_vs_control_adjust_bed")

gse28829_expr <- read_geo_series_matrix(project_path("data", "raw", "series_matrix", "GSE28829_series_matrix.txt.gz"))
gse28829_pheno <- load_pheno("GSE28829")
gse28829_pheno$group_std <- ifelse(gse28829_pheno$group_label == "advanced", "Advanced", "Early")
gse28829_pheno$keep_for_primary_contrast <- "yes"
gse28829_expr <- reorder_matrix(gse28829_expr, gse28829_pheno$sample_id)
write_tsv_gz(gse28829_expr, file.path(bulk_dir, "GSE28829_expr.tsv.gz"))
write_tsv(gse28829_pheno, file.path(bulk_dir, "GSE28829_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE28829", gse28829_expr, gse28829_pheno, "advanced_vs_early")

gse21545_expr <- read_geo_series_matrix(project_path("data", "raw", "series_matrix", "GSE21545_series_matrix.txt.gz"))
gse21545_pheno <- load_pheno("GSE21545")
gse21545_pheno$group_std <- ifelse(gse21545_pheno$subgroup == "TRUE", "IschemicEvent", "NoEvent")
gse21545_expr <- reorder_matrix(gse21545_expr, gse21545_pheno$sample_id)
write_tsv_gz(gse21545_expr, file.path(bulk_dir, "GSE21545_expr.tsv.gz"))
write_tsv(gse21545_pheno, file.path(bulk_dir, "GSE21545_pheno.tsv"))
summary_rows[[length(summary_rows) + 1]] <- build_summary_row("GSE21545", gse21545_expr, gse21545_pheno, "pbmc_ischemic_event_support")

summary_table <- do.call(rbind, summary_rows)
write_tsv(summary_table, project_path("res", "qc", "bulk", "bulk_matrix_prep_summary.tsv"))

cat("Bulk expression matrices written to data/processed/bulk\n")
print(summary_table)
