source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

ensure_project_dirs()

expr <- read_geo_series_matrix(project_path("data", "raw", "series_matrix", "GSE46394_series_matrix.txt.gz"))
pheno <- load_pheno("GSE46394")
pheno$group_std <- ifelse(pheno$group_label == "diseased", "Diseased", "Healthy")
expr <- expr[, c("feature_id", pheno$sample_id), drop = FALSE]

mapping <- read_tsv_auto(project_path("data", "processed", "annotation", "GPL13534_mapping.tsv.gz"))
expr_annot <- merge(expr, mapping[, c("feature_id", "gene_symbol", "entrez_id", "gene_title")], by = "feature_id", all.x = TRUE, sort = FALSE)
expr_annot <- expr_annot[, c("feature_id", "gene_symbol", "entrez_id", "gene_title", pheno$sample_id), drop = FALSE]

write_tsv_gz(expr_annot, project_path("data", "processed", "methylation", "GSE46394_beta.tsv.gz"))
write_tsv(pheno, project_path("data", "processed", "methylation", "GSE46394_pheno.tsv"))

summary_table <- data.frame(
  dataset_id = "GSE46394",
  n_cpg = nrow(expr_annot),
  n_samples = length(pheno$sample_id),
  n_diseased = sum(pheno$group_std == "Diseased"),
  n_healthy = sum(pheno$group_std == "Healthy"),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, project_path("res", "qc", "mechanism", "GSE46394_methylation_matrix_summary.tsv"))

cat("GSE46394 methylation matrix written to data/processed/methylation/GSE46394_beta.tsv.gz\n")
print(summary_table)
