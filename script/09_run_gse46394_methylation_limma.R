source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(limma))

ensure_project_dirs()

meth <- read_tsv_auto(project_path("data", "processed", "methylation", "GSE46394_beta.tsv.gz"))
pheno <- read_tsv_auto(project_path("data", "processed", "methylation", "GSE46394_pheno.tsv"))

sample_ids <- pheno$sample_id
beta <- as.matrix(meth[, sample_ids, drop = FALSE])
storage.mode(beta) <- "double"
rownames(beta) <- meth$feature_id

group <- factor(pheno$group_std, levels = c("Healthy", "Diseased"))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("Healthy", "Diseased")

fit <- lmFit(beta, design)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Diseased - Healthy, levels = design)))
tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$feature_id <- rownames(tt)
tt$gene_symbol <- meth$gene_symbol[match(tt$feature_id, meth$feature_id)]
tt$entrez_id <- meth$entrez_id[match(tt$feature_id, meth$feature_id)]
tt$gene_title <- meth$gene_title[match(tt$feature_id, meth$feature_id)]
tt$delta_beta <- tt$logFC
tt$direction <- ifelse(tt$delta_beta > 0, "Hyper_in_Diseased", "Hypo_in_Diseased")
tt$sig_fdr <- ifelse(tt$adj.P.Val < 0.05, "yes", "no")
tt$sig_fdr_beta <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$delta_beta) >= 0.10, "yes", "no")
tt <- tt[, c("feature_id", "gene_symbol", "entrez_id", "gene_title", "delta_beta", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "direction", "sig_fdr", "sig_fdr_beta")]

write_tsv_gz(tt, project_path("res", "tables", "mechanism", "GSE46394_methylation_limma.tsv.gz"))

summary_table <- data.frame(
  dataset_id = "GSE46394",
  n_tested = nrow(tt),
  n_sig_fdr = sum(tt$sig_fdr == "yes", na.rm = TRUE),
  n_sig_fdr_beta = sum(tt$sig_fdr_beta == "yes", na.rm = TRUE),
  top_feature = tt$feature_id[[1]],
  top_gene = tt$gene_symbol[[1]],
  top_delta_beta = tt$delta_beta[[1]],
  top_adj_p = tt$adj.P.Val[[1]],
  stringsAsFactors = FALSE
)
write_tsv(summary_table, project_path("res", "qc", "mechanism", "GSE46394_methylation_limma_summary.tsv"))

cat("GSE46394 methylation limma written to res/tables/mechanism/GSE46394_methylation_limma.tsv.gz\n")
print(summary_table)
