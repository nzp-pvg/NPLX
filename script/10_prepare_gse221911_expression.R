source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(readxl))

ensure_project_dirs()

raw_count_file <- project_path("data", "raw", "supplementary", "GSE221911", "GSE221911_Illumina_raw_read_count_177_pt.xlsx")
tpm_file <- project_path("data", "raw", "supplementary", "GSE221911", "GSE221911_Illumina_Merged_TPM_177_pt.xlsx")
tpm_df <- as.data.frame(read_excel(tpm_file), stringsAsFactors = FALSE)
names(tpm_df)[1:2] <- c("gene_symbol", "gene_length")

pheno <- load_pheno("GSE221911")
sample_ids <- pheno$title
tpm_df <- tpm_df[, c("gene_symbol", "gene_length", sample_ids), drop = FALSE]

write_tsv_gz(tpm_df, project_path("data", "processed", "bulk", "GSE221911_tpm.tsv.gz"))

summary_table <- data.frame(
  dataset_id = "GSE221911",
  n_genes = nrow(tpm_df),
  n_samples = length(sample_ids),
  n_low = sum(pheno$group_label == "LOW"),
  n_mid = sum(pheno$group_label == "MID"),
  n_cad = sum(pheno$group_label == "CAD"),
  stringsAsFactors = FALSE
)
write_tsv(summary_table, project_path("res", "qc", "bulk", "GSE221911_expression_prep_summary.tsv"))

cat("GSE221911 processed expression files written to data/processed/bulk\n")
print(summary_table)
