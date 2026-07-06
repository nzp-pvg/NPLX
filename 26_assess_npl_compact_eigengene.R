source("script/R/00_project_config.R")

ensure_project_dirs()

modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
compact_genes <- modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"]

prepare_bulk_dataset <- function(dataset_id) {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", paste0(dataset_id, "_gene_expr.tsv.gz")))
  pheno <- load_pheno(dataset_id)
  rownames(expr) <- expr$gene_symbol
  expr_mat <- as.matrix(expr[, !(colnames(expr) %in% c("gene_symbol", "feature_id")), drop = FALSE])
  storage.mode(expr_mat) <- "double"
  list(expr = expr_mat, pheno = pheno)
}

compute_pca_module <- function(expr_mat, genes) {
  genes <- intersect(genes, rownames(expr_mat))
  dat <- t(expr_mat[genes, , drop = FALSE])
  dat <- scale(dat)
  dat[is.na(dat)] <- 0
  pca <- prcomp(dat, center = FALSE, scale. = FALSE)
  loadings <- pca$rotation[, 1]
  scores <- pca$x[, 1]
  if (!is.na(loadings["NPL"]) && loadings["NPL"] < 0) {
    loadings <- -loadings
    scores <- -scores
  }
  list(
    genes = genes,
    scores = scores,
    loadings = loadings,
    var_explained = (pca$sdev[1]^2) / sum(pca$sdev^2)
  )
}

project_module_scores <- function(expr_mat, loadings) {
  genes <- intersect(names(loadings), rownames(expr_mat))
  dat <- t(expr_mat[genes, , drop = FALSE])
  dat <- scale(dat)
  dat[is.na(dat)] <- 0
  as.numeric(dat %*% loadings[genes])
}

group_effect <- function(scores, pheno, group_a, group_b) {
  df <- data.frame(sample_id = names(scores), score = as.numeric(scores), stringsAsFactors = FALSE)
  df <- merge(df, pheno, by = "sample_id", all.x = TRUE)
  a <- df$score[df$group_label == group_a]
  b <- df$score[df$group_label == group_b]
  delta <- mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE)
  pooled_sd <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
  cohen_d <- ifelse(is.na(pooled_sd) || pooled_sd == 0, NA_real_, delta / pooled_sd)
  data.frame(
    mean_a = mean(a, na.rm = TRUE),
    mean_b = mean(b, na.rm = TRUE),
    delta_b_minus_a = delta,
    cohen_d = cohen_d,
    wilcox_p = tryCatch(wilcox.test(b, a, exact = FALSE)$p.value, error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )
}

datasets <- list(
  GSE43292 = c("intact", "plaque"),
  GSE100927 = c("control", "atherosclerotic"),
  GSE28829 = c("early", "advanced")
)

prepared <- lapply(names(datasets), prepare_bulk_dataset)
names(prepared) <- names(datasets)

observed_rows <- list()
loading_rows <- list()
loo_rows <- list()

for (dataset_id in names(prepared)) {
  obj <- prepared[[dataset_id]]
  pca_res <- compute_pca_module(obj$expr, compact_genes)
  score_vec <- pca_res$scores
  names(score_vec) <- colnames(obj$expr)
  effect <- group_effect(score_vec, obj$pheno, datasets[[dataset_id]][1], datasets[[dataset_id]][2])
  observed_rows[[dataset_id]] <- cbind(
    data.frame(
      dataset_id = dataset_id,
      eigengene_source = paste0(dataset_id, "_native"),
      var_explained = pca_res$var_explained,
      stringsAsFactors = FALSE
    ),
    effect
  )

  loading_rows[[dataset_id]] <- data.frame(
    dataset_id = dataset_id,
    gene_symbol = names(pca_res$loadings),
    loading = as.numeric(pca_res$loadings),
    stringsAsFactors = FALSE
  )

  for (drop_gene in compact_genes) {
    sub_genes <- setdiff(compact_genes, drop_gene)
    loo_res <- compute_pca_module(obj$expr, sub_genes)
    loo_scores <- loo_res$scores
    if (cor(loo_scores, pca_res$scores) < 0) {
      loo_scores <- -loo_scores
    }
    loo_rows[[length(loo_rows) + 1]] <- data.frame(
      dataset_id = dataset_id,
      dropped_gene = drop_gene,
      score_cor = cor(loo_scores, pca_res$scores),
      stringsAsFactors = FALSE
    )
  }
}

native_summary <- do.call(rbind, observed_rows)
loading_table <- do.call(rbind, loading_rows)
loo_table <- do.call(rbind, loo_rows)
loo_summary <- aggregate(score_cor ~ dataset_id, data = loo_table, FUN = function(x) c(mean = mean(x), min = min(x)))
loo_summary <- data.frame(
  dataset_id = loo_summary$dataset_id,
  mean_loo_score_cor = loo_summary$score_cor[, "mean"],
  min_loo_score_cor = loo_summary$score_cor[, "min"],
  stringsAsFactors = FALSE
)

## Cross-dataset projection using GSE43292 as discovery reference
ref <- compute_pca_module(prepared[["GSE43292"]]$expr, compact_genes)
projection_rows <- list()
for (dataset_id in names(prepared)) {
  obj <- prepared[[dataset_id]]
  proj_scores <- project_module_scores(obj$expr, ref$loadings)
  names(proj_scores) <- colnames(obj$expr)
  effect <- group_effect(proj_scores, obj$pheno, datasets[[dataset_id]][1], datasets[[dataset_id]][2])
  projection_rows[[dataset_id]] <- cbind(
    data.frame(
      dataset_id = dataset_id,
      eigengene_source = "GSE43292_projected",
      var_explained = NA_real_,
      stringsAsFactors = FALSE
    ),
    effect
  )
}
projection_summary <- do.call(rbind, projection_rows)

## Loading consistency
loading_wide <- reshape(
  loading_table,
  idvar = "gene_symbol",
  timevar = "dataset_id",
  direction = "wide"
)
dataset_ids <- names(prepared)
consistency_rows <- list()
for (i in seq_along(dataset_ids)) {
  for (j in seq_along(dataset_ids)) {
    if (j <= i) {
      next
    }
    a <- loading_wide[[paste0("loading.", dataset_ids[i])]]
    b <- loading_wide[[paste0("loading.", dataset_ids[j])]]
    consistency_rows[[length(consistency_rows) + 1]] <- data.frame(
      dataset_a = dataset_ids[i],
      dataset_b = dataset_ids[j],
      loading_cor = cor(a, b),
      stringsAsFactors = FALSE
    )
  }
}
consistency_table <- do.call(rbind, consistency_rows)

write_tsv(native_summary, project_path("res", "tables", "mechanism", "npl_compact_native_eigengene_summary.tsv"))
write_tsv(projection_summary, project_path("res", "tables", "mechanism", "npl_compact_projected_eigengene_summary.tsv"))
write_tsv(loading_table, project_path("res", "tables", "mechanism", "npl_compact_loading_table.tsv"))
write_tsv(consistency_table, project_path("res", "tables", "mechanism", "npl_compact_loading_consistency.tsv"))
write_tsv(loo_table, project_path("res", "tables", "mechanism", "npl_compact_loo_gene_scores.tsv"))
write_tsv(loo_summary, project_path("res", "tables", "mechanism", "npl_compact_loo_summary.tsv"))

cat("NPL compact eigengene tables written to res/tables/mechanism\n")
print(native_summary)
cat("\nProjected summary:\n")
print(projection_summary)
cat("\nLoading consistency:\n")
print(consistency_table)
cat("\nLeave-one-gene-out summary:\n")
print(loo_summary)
