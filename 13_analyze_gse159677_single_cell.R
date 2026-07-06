source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(Matrix))

ensure_project_dirs()

untar_dir <- project_path("data", "raw", "single_cell", "GSE159677", "per_sample")
pheno <- load_pheno("GSE159677")
markers <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_marker_sets.tsv"))
mechanism_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "athero_mechanism_gene_sets.tsv"))

find_triplet_dir <- function(root_dir) {
  candidates <- list.files(root_dir, pattern = "matrix.mtx.gz$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) != 1) {
    stop("Expected exactly one matrix.mtx.gz under ", root_dir, ", found ", length(candidates))
  }
  dirname(candidates[[1]])
}

read_10x_triplet <- function(triplet_dir) {
  matrix_path <- file.path(triplet_dir, "matrix.mtx.gz")
  feature_path <- file.path(triplet_dir, "features.tsv.gz")
  barcode_path <- file.path(triplet_dir, "barcodes.tsv.gz")

  mat <- readMM(matrix_path)
  features <- read.delim(gzfile(feature_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(barcode_path), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  rownames(mat) <- make.unique(features$V2)
  colnames(mat) <- barcodes$V1
  mat
}

normalize_log1p <- function(mat) {
  lib <- Matrix::colSums(mat)
  lib[lib == 0] <- 1
  norm <- t(t(mat) / lib * 10000)
  log1p(norm)
}

score_gene_set <- function(norm_mat, genes) {
  genes <- intersect(genes, rownames(norm_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(norm_mat)))
  }
  Matrix::colMeans(norm_mat[genes, , drop = FALSE])
}

sample_dirs <- list.dirs(untar_dir, recursive = FALSE, full.names = TRUE)
if (length(sample_dirs) == 0) {
  stop("No extracted per-sample directories found under ", untar_dir)
}

cell_rows <- list()
sample_summary_rows <- list()
mechanism_sample_rows <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  gsm <- sub("_.*$", "", sample_name)
  triplet_dir <- find_triplet_dir(sample_dir)
  counts <- read_10x_triplet(triplet_dir)

  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)
  mt_genes <- grepl("^MT-", rownames(counts))
  mt_percent <- Matrix::colSums(counts[mt_genes, , drop = FALSE]) / pmax(n_count, 1) * 100
  keep <- n_feature >= 200 & n_feature <= 7500 & n_count >= 500 & mt_percent <= 20

  counts <- counts[, keep, drop = FALSE]
  n_count <- n_count[keep]
  n_feature <- n_feature[keep]
  mt_percent <- mt_percent[keep]
  norm <- normalize_log1p(counts)

  score_mat <- sapply(split(markers$gene_symbol, markers$cell_type), function(genes) score_gene_set(norm, genes))
  if (is.vector(score_mat)) {
    score_mat <- matrix(score_mat, ncol = 1, dimnames = list(colnames(norm), unique(markers$cell_type)))
  }
  pred_type <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
  top_score <- apply(score_mat, 1, max, na.rm = TRUE)
  pred_type[top_score < 0.15 | is.na(top_score)] <- "UNRESOLVED"

  sample_pheno <- pheno[pheno$sample_id == gsm, , drop = FALSE]
  location <- if (nrow(sample_pheno) == 1) sample_pheno$group_label else NA_character_

  mechanism_scores <- sapply(split(mechanism_sets$gene_symbol, mechanism_sets$gene_set), function(genes) score_gene_set(norm, genes))
  if (is.vector(mechanism_scores)) {
    mechanism_scores <- matrix(mechanism_scores, ncol = 1, dimnames = list(colnames(norm), unique(mechanism_sets$gene_set)))
  }

  cell_rows[[sample_name]] <- data.frame(
    sample_name = sample_name,
    sample_id = gsm,
    location = location,
    barcode = colnames(norm),
    n_count = as.numeric(n_count),
    n_feature = as.numeric(n_feature),
    mt_percent = as.numeric(mt_percent),
    predicted_cell_type = pred_type,
    top_marker_score = as.numeric(top_score),
    stringsAsFactors = FALSE
  )

  tab <- table(pred_type)
  sample_summary_rows[[sample_name]] <- data.frame(
    sample_name = sample_name,
    sample_id = gsm,
    location = location,
    cell_type = names(tab),
    n_cells = as.integer(tab),
    fraction = as.numeric(tab) / sum(tab),
    stringsAsFactors = FALSE
  )

  mechanism_sample_rows[[sample_name]] <- data.frame(
    sample_name = sample_name,
    sample_id = gsm,
    location = location,
    gene_set = colnames(mechanism_scores),
    mean_score = colMeans(mechanism_scores, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cell_table <- do.call(rbind, cell_rows)
sample_summary <- do.call(rbind, sample_summary_rows)
mechanism_summary <- do.call(rbind, mechanism_sample_rows)

group_summary <- aggregate(
  cbind(n_cells, fraction) ~ location + cell_type,
  data = sample_summary,
  FUN = mean
)
mechanism_group_summary <- aggregate(
  mean_score ~ location + gene_set,
  data = mechanism_summary,
  FUN = mean
)

write_tsv_gz(cell_table, project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
write_tsv(sample_summary, project_path("res", "tables", "mechanism", "gse159677_sample_celltype_summary.tsv"))
write_tsv(group_summary, project_path("res", "tables", "mechanism", "gse159677_group_celltype_summary.tsv"))
write_tsv(mechanism_summary, project_path("res", "tables", "mechanism", "gse159677_sample_mechanism_scores.tsv"))
write_tsv(mechanism_group_summary, project_path("res", "tables", "mechanism", "gse159677_group_mechanism_scores.tsv"))

qc_table <- data.frame(
  n_total_cells = nrow(cell_table),
  n_samples = length(unique(cell_table$sample_id)),
  stringsAsFactors = FALSE
)
write_tsv(qc_table, project_path("res", "qc", "mechanism", "gse159677_qc_summary.tsv"))

cat("GSE159677 single-cell summaries written to res/tables/mechanism\n")
print(group_summary)
print(mechanism_group_summary)
