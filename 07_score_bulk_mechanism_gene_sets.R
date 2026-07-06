source("script/R/00_project_config.R")

ensure_project_dirs()

score_gene_sets <- function(dataset_id, expr_file, pheno_file, case_label, control_label) {
  expr <- read_tsv_auto(expr_file)
  pheno <- read_tsv_auto(pheno_file)
  gene_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "athero_mechanism_gene_sets.tsv"))

  sample_ids <- pheno$sample_id
  mat <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- expr$gene_symbol
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0

  rows <- list()
  for (set_name in unique(gene_sets$gene_set)) {
    members <- intersect(gene_sets$gene_symbol[gene_sets$gene_set == set_name], rownames(mat_z))
    if (length(members) == 0) {
      next
    }
    scores <- colMeans(mat_z[members, , drop = FALSE])
    score_df <- data.frame(
      dataset_id = dataset_id,
      gene_set = set_name,
      sample_id = names(scores),
      score = as.numeric(scores),
      stringsAsFactors = FALSE
    )
    score_df <- merge(score_df, pheno[, c("sample_id", "group_std")], by = "sample_id", all.x = TRUE)
    case_scores <- score_df$score[score_df$group_std == case_label]
    control_scores <- score_df$score[score_df$group_std == control_label]
    test <- wilcox.test(case_scores, control_scores, exact = FALSE)
    summary_df <- data.frame(
      dataset_id = dataset_id,
      gene_set = set_name,
      case_group = case_label,
      control_group = control_label,
      n_case = length(case_scores),
      n_control = length(control_scores),
      case_mean = mean(case_scores, na.rm = TRUE),
      control_mean = mean(control_scores, na.rm = TRUE),
      delta_case_minus_control = mean(case_scores, na.rm = TRUE) - mean(control_scores, na.rm = TRUE),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- list(sample = score_df, summary = summary_df)
  }
  rows
}

specs <- list(
  list(dataset_id = "GSE43292", expr_file = project_path("data", "processed", "bulk_gene", "GSE43292_gene_expr.tsv.gz"), pheno_file = project_path("data", "processed", "bulk", "GSE43292_pheno.tsv"), case_label = "Plaque", control_label = "Intact"),
  list(dataset_id = "GSE100927", expr_file = project_path("data", "processed", "bulk_gene", "GSE100927_gene_expr.tsv.gz"), pheno_file = project_path("data", "processed", "bulk", "GSE100927_pheno.tsv"), case_label = "Atherosclerotic", control_label = "Control"),
  list(dataset_id = "GSE28829", expr_file = project_path("data", "processed", "bulk_gene", "GSE28829_gene_expr.tsv.gz"), pheno_file = project_path("data", "processed", "bulk", "GSE28829_pheno.tsv"), case_label = "Advanced", control_label = "Early"),
  list(dataset_id = "GSE21545_PBMC", expr_file = project_path("data", "processed", "bulk_gene", "GSE21545_gene_expr.tsv.gz"), pheno_file = project_path("data", "processed", "bulk", "GSE21545_pheno.tsv"), case_label = "IschemicEvent", control_label = "NoEvent")
)

sample_rows <- list()
summary_rows <- list()

for (spec in specs) {
  scored <- score_gene_sets(spec$dataset_id, spec$expr_file, spec$pheno_file, spec$case_label, spec$control_label)
  for (item in scored) {
    sample_rows[[length(sample_rows) + 1]] <- item$sample
    summary_rows[[length(summary_rows) + 1]] <- item$summary
  }
}

sample_table <- do.call(rbind, sample_rows)
summary_table <- do.call(rbind, summary_rows)
summary_table$adj_p_value <- p.adjust(summary_table$p_value, method = "BH")

write_tsv_gz(sample_table, project_path("res", "tables", "mechanism", "bulk_mechanism_scores.tsv.gz"))
write_tsv(summary_table, project_path("res", "tables", "mechanism", "bulk_mechanism_score_summary.tsv"))

cat("Bulk mechanism scores written to res/tables/mechanism\n")
print(summary_table)
