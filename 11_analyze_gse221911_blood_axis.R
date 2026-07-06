source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(limma)
})

ensure_project_dirs()

expr <- read_tsv_auto(project_path("data", "processed", "bulk", "GSE221911_tpm.tsv.gz"))
pheno <- load_pheno("GSE221911")
pheno <- pheno[pheno$group_label %in% c("LOW", "MID", "CAD"), , drop = FALSE]

sample_ids <- as.character(pheno$title)
expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
storage.mode(expr_matrix) <- "double"
rownames(expr_matrix) <- make.unique(expr$gene_symbol)
expr_matrix <- log2(expr_matrix + 1)

group <- factor(pheno$group_label, levels = c("LOW", "MID", "CAD"))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("LOW", "MID", "CAD")

fit <- lmFit(expr_matrix, design)
contrast_matrix <- makeContrasts(
  MID_vs_LOW = MID - LOW,
  CAD_vs_LOW = CAD - LOW,
  CAD_vs_MID = CAD - MID,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

result_rows <- list()
for (coef_name in colnames(contrast_matrix)) {
  tt <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt$contrast <- coef_name
  tt$gene_symbol <- rownames(tt)
  tt$dataset_id <- "GSE221911"
  result_rows[[coef_name]] <- tt[, c("dataset_id", "contrast", "gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}

limma_table <- do.call(rbind, result_rows)
write_tsv_gz(limma_table, project_path("res", "tables", "bulk", "GSE221911_blood_axis_limma.tsv.gz"))

summary_table <- do.call(
  rbind,
  lapply(split(limma_table, limma_table$contrast), function(df) {
    data.frame(
      contrast = unique(df$contrast),
      n_tested = nrow(df),
      n_sig_fdr = sum(df$adj.P.Val < 0.05, na.rm = TRUE),
      top_gene = df$gene_symbol[[1]],
      top_logFC = df$logFC[[1]],
      top_adj_p = df$adj.P.Val[[1]],
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(summary_table, project_path("res", "qc", "bulk", "GSE221911_blood_axis_limma_summary.tsv"))

cat("GSE221911 blood-axis limma written to res/tables/bulk/GSE221911_blood_axis_limma.tsv.gz\n")
print(summary_table)
