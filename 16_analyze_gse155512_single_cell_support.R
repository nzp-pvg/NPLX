source("script/R/00_project_config.R")

ensure_project_dirs()

pheno <- load_pheno("GSE155512")
markers <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_marker_sets.tsv"))

normalize_log1p_dense <- function(mat) {
  lib <- colSums(mat)
  lib[lib == 0] <- 1
  log1p(t(t(mat) / lib * 10000))
}

score_gene_set_dense <- function(norm_mat, genes) {
  genes <- intersect(genes, rownames(norm_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(norm_mat)))
  }
  colMeans(norm_mat[genes, , drop = FALSE])
}

sample_files <- list.files(project_path("data", "raw", "single_cell", "GSE155512"), pattern = "_matrix\\.txt\\.gz$", full.names = TRUE)

cell_rows <- list()
sample_summary_rows <- list()

for (file in sample_files) {
  gsm <- sub("_.*$", "", basename(file))
  sample_pheno <- pheno[pheno$sample_id == gsm, , drop = FALSE]
  symptom <- if (nrow(sample_pheno) == 1) sample_pheno$group_label else NA_character_

  mat <- read.delim(gzfile(file), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gene_symbol <- mat[[1]]
  expr <- as.matrix(mat[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- make.unique(gene_symbol)

  n_count <- colSums(expr)
  n_feature <- colSums(expr > 0)
  keep <- n_feature >= 200 & n_count >= 500
  expr <- expr[, keep, drop = FALSE]
  n_count <- n_count[keep]
  n_feature <- n_feature[keep]
  norm <- normalize_log1p_dense(expr)

  score_mat <- sapply(split(markers$gene_symbol, markers$cell_type), function(genes) score_gene_set_dense(norm, genes))
  if (is.vector(score_mat)) {
    score_mat <- matrix(score_mat, ncol = 1, dimnames = list(colnames(norm), unique(markers$cell_type)))
  }
  pred_type <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
  top_score <- apply(score_mat, 1, max, na.rm = TRUE)
  pred_type[top_score < 0.15 | is.na(top_score)] <- "UNRESOLVED"

  cell_rows[[gsm]] <- data.frame(
    sample_id = gsm,
    symptom = symptom,
    barcode = colnames(norm),
    n_count = as.numeric(n_count),
    n_feature = as.numeric(n_feature),
    predicted_cell_type = pred_type,
    top_marker_score = as.numeric(top_score),
    stringsAsFactors = FALSE
  )

  tab <- table(pred_type)
  sample_summary_rows[[gsm]] <- data.frame(
    sample_id = gsm,
    symptom = symptom,
    cell_type = names(tab),
    n_cells = as.integer(tab),
    fraction = as.numeric(tab) / sum(tab),
    stringsAsFactors = FALSE
  )
}

cell_table <- do.call(rbind, cell_rows)
sample_summary <- do.call(rbind, sample_summary_rows)
group_summary <- aggregate(cbind(n_cells, fraction) ~ symptom + cell_type, data = sample_summary, FUN = mean)

write_tsv_gz(cell_table, project_path("res", "tables", "mechanism", "gse155512_cell_level_typing.tsv.gz"))
write_tsv(sample_summary, project_path("res", "tables", "mechanism", "gse155512_sample_celltype_summary.tsv"))
write_tsv(group_summary, project_path("res", "tables", "mechanism", "gse155512_group_celltype_summary.tsv"))

cat("GSE155512 single-cell support summaries written to res/tables/mechanism\n")
print(group_summary)
