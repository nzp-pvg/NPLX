source("script/R/00_project_config.R")

ensure_project_dirs()

set.seed(20260330)

modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
observed_genes <- modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"]

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
  list(loadings = loadings, scores = scores)
}

project_module_scores <- function(expr_mat, loadings) {
  genes <- intersect(names(loadings), rownames(expr_mat))
  dat <- t(expr_mat[genes, , drop = FALSE])
  dat <- scale(dat)
  dat[is.na(dat)] <- 0
  as.numeric(dat %*% loadings[genes])
}

group_effect_size <- function(scores, pheno, group_a, group_b) {
  df <- data.frame(sample_id = names(scores), score = as.numeric(scores), stringsAsFactors = FALSE)
  df <- merge(df, pheno, by = "sample_id", all.x = TRUE)
  a <- df$score[df$group_label == group_a]
  b <- df$score[df$group_label == group_b]
  delta <- mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE)
  pooled_sd <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  delta / pooled_sd
}

datasets <- list(
  GSE43292 = c("intact", "plaque"),
  GSE100927 = c("control", "atherosclerotic"),
  GSE28829 = c("early", "advanced")
)

prepared <- lapply(names(datasets), prepare_bulk_dataset)
names(prepared) <- names(datasets)
common_genes <- Reduce(intersect, lapply(prepared, function(x) rownames(x$expr)))
common_genes <- setdiff(common_genes, observed_genes)
module_size <- length(observed_genes)
n_iter <- 1000

score_module_transfer <- function(genes) {
  ref <- compute_pca_module(prepared[["GSE43292"]]$expr, genes)
  ref_scores <- ref$scores
  names(ref_scores) <- colnames(prepared[["GSE43292"]]$expr)
  ref_effect <- group_effect_size(ref_scores, prepared[["GSE43292"]]$pheno, datasets[["GSE43292"]][1], datasets[["GSE43292"]][2])
  if (!is.na(ref_effect) && ref_effect < 0) {
    ref$loadings <- -ref$loadings
    ref_scores <- -ref_scores
    ref_effect <- -ref_effect
  }

  effects <- c(GSE43292 = ref_effect)
  for (dataset_id in c("GSE100927", "GSE28829")) {
    proj_scores <- project_module_scores(prepared[[dataset_id]]$expr, ref$loadings)
    names(proj_scores) <- colnames(prepared[[dataset_id]]$expr)
    effects[dataset_id] <- group_effect_size(proj_scores, prepared[[dataset_id]]$pheno, datasets[[dataset_id]][1], datasets[[dataset_id]][2])
  }

  data.frame(
    mean_transfer_d = mean(effects, na.rm = TRUE),
    min_transfer_d = min(effects, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

observed_score <- score_module_transfer(observed_genes)

null_rows <- vector("list", n_iter)
for (iter in seq_len(n_iter)) {
  genes <- sample(common_genes, size = module_size, replace = FALSE)
  score <- score_module_transfer(genes)
  null_rows[[iter]] <- data.frame(
    iter = iter,
    mean_transfer_d = score$mean_transfer_d,
    min_transfer_d = score$min_transfer_d,
    stringsAsFactors = FALSE
  )
}

null_table <- do.call(rbind, null_rows)
benchmark_summary <- data.frame(
  module_name = "NPL_FOAM_MACROPHAGE_COMPACT",
  observed_mean_transfer_d = observed_score$mean_transfer_d,
  observed_min_transfer_d = observed_score$min_transfer_d,
  null_mean_mean_transfer_d = mean(null_table$mean_transfer_d, na.rm = TRUE),
  null_sd_mean_transfer_d = sd(null_table$mean_transfer_d, na.rm = TRUE),
  empirical_p_mean_transfer = mean(null_table$mean_transfer_d >= observed_score$mean_transfer_d, na.rm = TRUE),
  empirical_p_min_transfer = mean(null_table$min_transfer_d >= observed_score$min_transfer_d, na.rm = TRUE),
  mean_transfer_percentile = mean(null_table$mean_transfer_d <= observed_score$mean_transfer_d, na.rm = TRUE),
  min_transfer_percentile = mean(null_table$min_transfer_d <= observed_score$min_transfer_d, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write_tsv(null_table, project_path("res", "tables", "mechanism", "npl_compact_null_benchmark.tsv"))
write_tsv(benchmark_summary, project_path("res", "tables", "mechanism", "npl_compact_null_benchmark_summary.tsv"))

cat("NPL compact null benchmark written to res/tables/mechanism\n")
print(benchmark_summary)
