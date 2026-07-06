source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
})

ensure_project_dirs()

expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", "GSE21545_gene_expr.tsv.gz"))
pheno <- read_tsv_auto(project_path("data", "processed", "bulk", "GSE21545_pheno.tsv"))
modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))

expr_mat <- as.matrix(expr[, pheno$sample_id, drop = FALSE])
rownames(expr_mat) <- expr$gene_symbol

score_gene_set_z <- function(expr_mat, genes) {
  genes <- intersect(unique(genes), rownames(expr_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(expr_mat)))
  }
  sub <- expr_mat[genes, , drop = FALSE]
  sub_z <- t(scale(t(sub)))
  sub_z[!is.finite(sub_z)] <- 0
  colMeans(sub_z)
}

compare_groups <- function(values, groups, case_label, control_label) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- groups[keep]
  case_vals <- values[groups == case_label]
  control_vals <- values[groups == control_label]
  if (length(case_vals) < 3 || length(control_vals) < 3) {
    return(data.frame(
      n_case = length(case_vals),
      n_control = length(control_vals),
      mean_case = mean(case_vals),
      mean_control = mean(control_vals),
      delta_case_minus_control = mean(case_vals) - mean(control_vals),
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  wt <- wilcox.test(case_vals, control_vals, exact = FALSE)
  data.frame(
    n_case = length(case_vals),
    n_control = length(control_vals),
    mean_case = mean(case_vals),
    mean_control = mean(control_vals),
    delta_case_minus_control = mean(case_vals) - mean(control_vals),
    p_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}

build_signature <- function(signature_name, genes) {
  data.frame(
    sample_id = pheno$sample_id,
    signature_name = signature_name,
    score = score_gene_set_z(expr_mat, genes)[pheno$sample_id],
    stringsAsFactors = FALSE
  )
}

sig_npl <- build_signature("NPL_FOAM_MACROPHAGE_COMPACT", modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"])
sig_ox <- build_signature("OXIDATIVE_RESPONSE", c("HMOX1", "GPX1", "NQO1", "FTL", "FTH1"))
sig_endo <- build_signature("ENDOTHELIAL_ACTIVATION", c("VCAM1", "CCL2", "VWF"))

score_table <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), list(sig_npl, sig_ox, sig_endo))
colnames(score_table) <- c("sample_id", "signature_name.x", "npl_score", "signature_name.y", "oxidative_score", "signature_name", "endothelial_score")
score_table <- merge(pheno, score_table[, c("sample_id", "npl_score", "oxidative_score", "endothelial_score")], by = "sample_id", all.x = TRUE)

long_scores <- rbind(
  data.frame(score_table[, c("sample_id", "tissue_context", "group_std")], signature_name = "NPL_FOAM_MACROPHAGE_COMPACT", score = score_table$npl_score, stringsAsFactors = FALSE),
  data.frame(score_table[, c("sample_id", "tissue_context", "group_std")], signature_name = "OXIDATIVE_RESPONSE", score = score_table$oxidative_score, stringsAsFactors = FALSE),
  data.frame(score_table[, c("sample_id", "tissue_context", "group_std")], signature_name = "ENDOTHELIAL_ACTIVATION", score = score_table$endothelial_score, stringsAsFactors = FALSE)
)

ctx_rows <- list()
for (ctx in sort(unique(long_scores$tissue_context))) {
  for (sig in unique(long_scores$signature_name)) {
    sub <- long_scores[long_scores$tissue_context == ctx & long_scores$signature_name == sig, ]
    cmp <- compare_groups(sub$score, sub$group_std, case_label = "IschemicEvent", control_label = "NoEvent")
    cmp$tissue_context <- ctx
    cmp$signature_name <- sig
    ctx_rows[[paste(ctx, sig, sep = "::")]] <- cmp
  }
}
ctx_summary <- do.call(rbind, ctx_rows)
ctx_summary$adj_p_value <- p.adjust(ctx_summary$p_value, method = "BH")

tissue_rows <- list()
for (sig in unique(long_scores$signature_name)) {
  sub <- long_scores[long_scores$signature_name == sig, ]
  cmp <- compare_groups(sub$score, sub$tissue_context, case_label = "carotid plaque", control_label = "peripheral blood mononuclear cells")
  cmp$contrast <- "plaque_vs_pbmc"
  cmp$signature_name <- sig
  tissue_rows[[sig]] <- cmp
}
tissue_summary <- do.call(rbind, tissue_rows)
tissue_summary$adj_p_value <- p.adjust(tissue_summary$p_value, method = "BH")

pbmc <- subset(score_table, tissue_context == "peripheral blood mononuclear cells")
pbmc$event01 <- ifelse(pbmc$group_std == "IschemicEvent", 1, 0)

fit_npl <- glm(event01 ~ npl_score, data = pbmc, family = binomial())
fit_npl_ox <- glm(event01 ~ npl_score + oxidative_score, data = pbmc, family = binomial())
fit_all <- glm(event01 ~ npl_score + oxidative_score + endothelial_score, data = pbmc, family = binomial())

extract_model <- function(fit, model_name) {
  cf <- coef(summary(fit))
  out <- data.frame(
    model_name = model_name,
    term = rownames(cf),
    estimate = cf[, "Estimate"],
    std_error = cf[, "Std. Error"],
    z_value = cf[, "z value"],
    p_value = cf[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  preds <- predict(fit, type = "response")
  auc <- as.numeric(pROC::roc(pbmc$event01, preds, quiet = TRUE)$auc)
  list(coefficients = out, auc = auc, preds = preds)
}

mod_npl <- extract_model(fit_npl, "PBMC_event_NPL")
mod_npl_ox <- extract_model(fit_npl_ox, "PBMC_event_NPL_OX")
mod_all <- extract_model(fit_all, "PBMC_event_NPL_OX_ENDO")

model_coeffs <- rbind(mod_npl$coefficients, mod_npl_ox$coefficients, mod_all$coefficients)
model_summary <- data.frame(
  model_name = c("PBMC_event_NPL", "PBMC_event_NPL_OX", "PBMC_event_NPL_OX_ENDO"),
  auc = c(mod_npl$auc, mod_npl_ox$auc, mod_all$auc),
  stringsAsFactors = FALSE
)

model_predictions <- data.frame(
  sample_id = pbmc$sample_id,
  group_std = pbmc$group_std,
  npl_only = mod_npl$preds,
  npl_oxidative = mod_npl_ox$preds,
  full_model = mod_all$preds,
  stringsAsFactors = FALSE
)

write_tsv(score_table, project_path("res", "tables", "bulk", "gse21545_translational_signature_scores.tsv"))
write_tsv(ctx_summary, project_path("res", "qc", "bulk", "gse21545_signature_context_summary.tsv"))
write_tsv(tissue_summary, project_path("res", "qc", "bulk", "gse21545_signature_tissue_summary.tsv"))
write_tsv(model_coeffs, project_path("res", "tables", "bulk", "gse21545_pbmc_event_model_coefficients.tsv"))
write_tsv(model_summary, project_path("res", "qc", "bulk", "gse21545_pbmc_event_model_summary.tsv"))
write_tsv(model_predictions, project_path("res", "tables", "bulk", "gse21545_pbmc_event_model_predictions.tsv"))

p_ctx <- ggplot(long_scores, aes(group_std, score, fill = group_std)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.size = 0.25) +
  facet_grid(signature_name ~ tissue_context, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(
    title = "GSE21545 predefined signatures by tissue context and ischemic event status",
    x = NULL,
    y = "Signature score (gene-wise z mean)"
  )

roc_rows <- rbind(
  data.frame(fpr = 1 - pROC::roc(pbmc$event01, mod_npl$preds, quiet = TRUE)$specificities,
             tpr = pROC::roc(pbmc$event01, mod_npl$preds, quiet = TRUE)$sensitivities,
             model_name = "PBMC_event_NPL"),
  data.frame(fpr = 1 - pROC::roc(pbmc$event01, mod_npl_ox$preds, quiet = TRUE)$specificities,
             tpr = pROC::roc(pbmc$event01, mod_npl_ox$preds, quiet = TRUE)$sensitivities,
             model_name = "PBMC_event_NPL_OX"),
  data.frame(fpr = 1 - pROC::roc(pbmc$event01, mod_all$preds, quiet = TRUE)$specificities,
             tpr = pROC::roc(pbmc$event01, mod_all$preds, quiet = TRUE)$sensitivities,
             model_name = "PBMC_event_NPL_OX_ENDO")
)

p_roc <- ggplot(roc_rows, aes(fpr, tpr, color = model_name)) +
  geom_path(linewidth = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey50") +
  theme_bw(base_size = 10) +
  labs(
    title = "GSE21545 PBMC ischemic-event support",
    subtitle = paste(
      sprintf("%s AUC=%.3f", model_summary$model_name, model_summary$auc),
      collapse = " | "
    ),
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Model"
  )

ggsave(project_path("figure", "export", "bulk", "gse21545_translational_signature_contexts.pdf"), p_ctx, width = 9.2, height = 6.6)
ggsave(project_path("figure", "export", "bulk", "gse21545_pbmc_event_roc.pdf"), p_roc, width = 7.2, height = 5.2)

cat("GSE21545 translational support outputs written to res/tables/bulk, res/qc/bulk, and figure/export/bulk\n")
print(ctx_summary)
print(tissue_summary)
print(model_summary)
