source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(limma))

ensure_project_dirs()

run_limma_gse43292 <- function() {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", "GSE43292_gene_expr.tsv.gz"))
  pheno <- read_tsv_auto(project_path("data", "processed", "bulk", "GSE43292_pheno.tsv"))

  sample_ids <- pheno$sample_id
  expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- make.unique(expr$gene_symbol)

  group <- factor(pheno$group_std, levels = c("Intact", "Plaque"))
  patient <- factor(pheno$paired_id)
  design <- model.matrix(~ 0 + group + patient)
  colnames(design)[1:2] <- c("Intact", "Plaque")

  fit <- lmFit(expr_matrix, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Plaque - Intact, levels = design)))
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene_symbol <- expr$gene_symbol[match(rownames(tt), rownames(expr_matrix))]
  tt$feature_id <- expr$feature_id[match(rownames(tt), rownames(expr_matrix))]
  tt$dataset_id <- "GSE43292"
  tt$contrast <- "Plaque_vs_Intact"
  tt$sign_direction <- ifelse(tt$logFC > 0, "Up_in_Plaque", "Down_in_Plaque")
  tt$sig_fdr <- ifelse(tt$adj.P.Val < 0.05, "yes", "no")
  tt$sig_fdr_fc <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5, "yes", "no")
  tt[, c("dataset_id", "contrast", "gene_symbol", "feature_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "sign_direction", "sig_fdr", "sig_fdr_fc")]
}

run_limma_gse100927 <- function() {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", "GSE100927_gene_expr.tsv.gz"))
  pheno <- read_tsv_auto(project_path("data", "processed", "bulk", "GSE100927_pheno.tsv"))
  pheno <- pheno[pheno$group_std %in% c("Atherosclerotic", "Control"), , drop = FALSE]

  sample_ids <- pheno$sample_id
  expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- make.unique(expr$gene_symbol)

  group <- factor(pheno$group_std, levels = c("Control", "Atherosclerotic"))
  bed <- factor(pheno$vascular_bed)
  design <- model.matrix(~ 0 + group + bed)
  colnames(design) <- make.names(colnames(design))
  colnames(design)[1:2] <- c("Control", "Atherosclerotic")

  fit <- lmFit(expr_matrix, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Atherosclerotic - Control, levels = design)))
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene_symbol <- expr$gene_symbol[match(rownames(tt), rownames(expr_matrix))]
  tt$feature_id <- expr$feature_id[match(rownames(tt), rownames(expr_matrix))]
  tt$dataset_id <- "GSE100927"
  tt$contrast <- "Atherosclerotic_vs_Control_adjusted_vascular_bed"
  tt$sign_direction <- ifelse(tt$logFC > 0, "Up_in_Atherosclerotic", "Down_in_Atherosclerotic")
  tt$sig_fdr <- ifelse(tt$adj.P.Val < 0.05, "yes", "no")
  tt$sig_fdr_fc <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5, "yes", "no")
  tt[, c("dataset_id", "contrast", "gene_symbol", "feature_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "sign_direction", "sig_fdr", "sig_fdr_fc")]
}

run_limma_gse28829 <- function() {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", "GSE28829_gene_expr.tsv.gz"))
  pheno <- read_tsv_auto(project_path("data", "processed", "bulk", "GSE28829_pheno.tsv"))

  sample_ids <- pheno$sample_id
  expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- make.unique(expr$gene_symbol)

  group <- factor(pheno$group_std, levels = c("Early", "Advanced"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("Early", "Advanced")

  fit <- lmFit(expr_matrix, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(Advanced - Early, levels = design)))
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene_symbol <- expr$gene_symbol[match(rownames(tt), rownames(expr_matrix))]
  tt$feature_id <- expr$feature_id[match(rownames(tt), rownames(expr_matrix))]
  tt$dataset_id <- "GSE28829"
  tt$contrast <- "Advanced_vs_Early"
  tt$sign_direction <- ifelse(tt$logFC > 0, "Up_in_Advanced", "Down_in_Advanced")
  tt$sig_fdr <- ifelse(tt$adj.P.Val < 0.05, "yes", "no")
  tt$sig_fdr_fc <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5, "yes", "no")
  tt[, c("dataset_id", "contrast", "gene_symbol", "feature_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "sign_direction", "sig_fdr", "sig_fdr_fc")]
}

run_limma_gse21545 <- function() {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", "GSE21545_gene_expr.tsv.gz"))
  pheno <- read_tsv_auto(project_path("data", "processed", "bulk", "GSE21545_pheno.tsv"))
  pheno <- pheno[pheno$group_label == "peripheral blood mononuclear cells", , drop = FALSE]

  sample_ids <- pheno$sample_id
  expr_matrix <- as.matrix(expr[, sample_ids, drop = FALSE])
  storage.mode(expr_matrix) <- "double"
  rownames(expr_matrix) <- make.unique(expr$gene_symbol)

  group <- factor(pheno$group_std, levels = c("NoEvent", "IschemicEvent"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("NoEvent", "IschemicEvent")

  fit <- lmFit(expr_matrix, design)
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(IschemicEvent - NoEvent, levels = design)))
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$gene_symbol <- expr$gene_symbol[match(rownames(tt), rownames(expr_matrix))]
  tt$feature_id <- expr$feature_id[match(rownames(tt), rownames(expr_matrix))]
  tt$dataset_id <- "GSE21545_PBMC"
  tt$contrast <- "IschemicEvent_vs_NoEvent"
  tt$sign_direction <- ifelse(tt$logFC > 0, "Up_in_IschemicEvent", "Down_in_IschemicEvent")
  tt$sig_fdr <- ifelse(tt$adj.P.Val < 0.05, "yes", "no")
  tt$sig_fdr_fc <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5, "yes", "no")
  tt[, c("dataset_id", "contrast", "gene_symbol", "feature_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "sign_direction", "sig_fdr", "sig_fdr_fc")]
}

results <- list(
  GSE43292 = run_limma_gse43292(),
  GSE100927 = run_limma_gse100927(),
  GSE28829 = run_limma_gse28829(),
  GSE21545_PBMC = run_limma_gse21545()
)

summary_rows <- list()

for (nm in names(results)) {
  df <- results[[nm]]
  write_tsv_gz(df, project_path("res", "tables", "bulk", paste0(nm, "_limma.tsv.gz")))
  summary_rows[[nm]] <- data.frame(
    dataset_id = nm,
    n_tested = nrow(df),
    n_sig_fdr = sum(df$sig_fdr == "yes"),
    n_sig_fdr_fc = sum(df$sig_fdr_fc == "yes"),
    top_gene = df$gene_symbol[[1]],
    top_logFC = df$logFC[[1]],
    top_adj_p = df$adj.P.Val[[1]],
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
write_tsv(summary_table, project_path("res", "qc", "bulk", "bulk_limma_summary.tsv"))

cat("Bulk limma results written to res/tables/bulk\n")
print(summary_table)
