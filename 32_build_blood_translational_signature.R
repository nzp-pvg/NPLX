source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(ggplot2)
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

plaque_core_signature <- c("NPL", "FABP5", "GPNMB", "APOC1", "PLA2G7", "SPP1", "CD36", "APOE", "LGALS3")
blood_translational_signature <- c("NPL", "FABP5", "PLA2G7", "SPP1", "CD36", "LGALS3", "HMOX1", "CCL2", "VWF")
oxidative_signature <- c("HMOX1", "GPX1", "NQO1", "FTL", "FTH1")

score_signature <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  colMeans(mat[genes, , drop = FALSE])
}

score_df <- data.frame(
  sample_id = sample_ids,
  group_label = pheno$group_label,
  plaque_core_signature = score_signature(expr_matrix, plaque_core_signature),
  blood_translational_signature = score_signature(expr_matrix, blood_translational_signature),
  oxidative_signature = score_signature(expr_matrix, oxidative_signature),
  stringsAsFactors = FALSE
)
score_df$group_num <- c(LOW = 0, MID = 1, CAD = 2)[score_df$group_label]

sig_rows <- list()
for (score_name in c("plaque_core_signature", "blood_translational_signature", "oxidative_signature")) {
  vals <- score_df[[score_name]]
  sig_rows[[score_name]] <- data.frame(
    score_name = score_name,
    spearman_group_trend = suppressWarnings(cor(vals, score_df$group_num, method = "spearman")),
    mean_low = mean(vals[score_df$group_label == "LOW"], na.rm = TRUE),
    mean_mid = mean(vals[score_df$group_label == "MID"], na.rm = TRUE),
    mean_cad = mean(vals[score_df$group_label == "CAD"], na.rm = TRUE),
    low_vs_cad_p = wilcox.test(vals[score_df$group_label == "LOW"], vals[score_df$group_label == "CAD"], exact = FALSE)$p.value,
    stringsAsFactors = FALSE
  )
}
signature_summary <- do.call(rbind, sig_rows)

cad_low <- score_df[score_df$group_label %in% c("LOW", "CAD"), , drop = FALSE]
cad_low$y <- ifelse(cad_low$group_label == "CAD", 1, 0)
candidate_genes <- intersect(blood_translational_signature, rownames(expr_matrix))
x <- t(expr_matrix[candidate_genes, cad_low$sample_id, drop = FALSE])
y <- cad_low$y

set.seed(42)
cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 0.5, type.measure = "auc", nfolds = 5)
pred_prob <- as.numeric(predict(cvfit, newx = x, s = "lambda.1se", type = "response"))
roc_obj <- pROC::roc(response = y, predictor = pred_prob, quiet = TRUE)
coef_mat <- as.matrix(coef(cvfit, s = "lambda.1se"))
coef_df <- data.frame(
  feature = rownames(coef_mat),
  coefficient = coef_mat[, 1],
  stringsAsFactors = FALSE
)
coef_df <- coef_df[coef_df$coefficient != 0, , drop = FALSE]

cad_low$translational_classifier_prob <- pred_prob
cad_low$blood_translational_signature_z <- scale(cad_low$blood_translational_signature)[, 1]

write_tsv(score_df, project_path("res", "tables", "bulk", "gse221911_signature_scores.tsv"))
write_tsv(signature_summary, project_path("res", "qc", "bulk", "gse221911_signature_summary.tsv"))
write_tsv(coef_df, project_path("res", "tables", "bulk", "gse221911_translational_classifier_coefficients.tsv"))
write_tsv(cad_low, project_path("res", "tables", "bulk", "gse221911_translational_classifier_predictions.tsv"))
write_tsv(
  data.frame(
    auc = as.numeric(pROC::auc(roc_obj)),
    n_low = sum(cad_low$group_label == "LOW"),
    n_cad = sum(cad_low$group_label == "CAD"),
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    stringsAsFactors = FALSE
  ),
  project_path("res", "qc", "bulk", "gse221911_translational_classifier_summary.tsv")
)

p_score <- ggplot(score_df, aes(group_label, blood_translational_signature, fill = group_label)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = "Blood translational signature in GSE221911", x = NULL, y = "Mean log2(TPM+1)")

p_prob <- ggplot(cad_low, aes(group_label, translational_classifier_prob, fill = group_label)) +
  geom_boxplot(outlier.size = 0.3) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = paste0("CAD vs LOW classifier (AUC=", sprintf("%.3f", as.numeric(pROC::auc(roc_obj))), ")"),
       x = NULL, y = "Predicted CAD probability")

ggsave(project_path("figure", "export", "bulk", "gse221911_blood_translational_signature.pdf"), p_score, width = 5.5, height = 4.2)
ggsave(project_path("figure", "export", "bulk", "gse221911_blood_translational_classifier.pdf"), p_prob, width = 5.5, height = 4.2)

cat("Blood translational signature outputs written to res/tables/bulk and figure/export/bulk\n")
print(signature_summary)
print(coef_df)
cat("AUC: ", as.numeric(pROC::auc(roc_obj)), "\n", sep = "")
