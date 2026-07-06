source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(Matrix))

ensure_project_dirs()

find_triplet_dir <- function(root_dir) {
  candidates <- list.files(root_dir, pattern = "matrix.mtx.gz$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) != 1) {
    stop("Expected exactly one matrix.mtx.gz under ", root_dir, ", found ", length(candidates))
  }
  dirname(candidates[[1]])
}

read_10x_triplet <- function(triplet_dir) {
  mat <- readMM(file.path(triplet_dir, "matrix.mtx.gz"))
  features <- read.delim(gzfile(file.path(triplet_dir, "features.tsv.gz")), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(file.path(triplet_dir, "barcodes.tsv.gz")), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
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

normalize_log1p_dense <- function(mat) {
  lib <- colSums(mat)
  lib[lib == 0] <- 1
  log1p(t(t(mat) / lib * 10000))
}

score_gene_set <- function(norm_mat, genes) {
  genes <- intersect(genes, rownames(norm_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(norm_mat)))
  }
  Matrix::colMeans(norm_mat[genes, , drop = FALSE])
}

infer_patient_id <- function(sample_id) {
  ifelse(
    sample_id %in% c("GSM4837523", "GSM4837524"),
    "P1",
    ifelse(sample_id %in% c("GSM4837525", "GSM4837526"), "P2", "P3")
  )
}

modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))
cell159677 <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
cell155512 <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse155512_cell_level_typing.tsv.gz"))
meth <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE46394_methylation_limma.tsv.gz"))

## GSE159677 macrophage layer
sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
mats <- list()
meta_rows <- list()
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  if (!grepl("^GSM", sample_name)) {
    next
  }
  counts <- read_10x_triplet(find_triplet_dir(sample_dir))
  keep <- cell159677$barcode[cell159677$sample_name == sample_name & cell159677$predicted_cell_type == "MACROPHAGE"]
  keep <- intersect(colnames(counts), keep)
  if (length(keep) == 0) {
    next
  }
  counts <- counts[, keep, drop = FALSE]
  colnames(counts) <- paste(sample_name, keep, sep = "|")
  mats[[sample_name]] <- counts
  sample_id <- unique(cell159677$sample_id[cell159677$sample_name == sample_name])[1]
  meta_rows[[sample_name]] <- data.frame(
    cell_id = colnames(counts),
    sample_name = sample_name,
    sample_id = sample_id,
    patient_id = infer_patient_id(sample_id),
    location = unique(cell159677$location[cell159677$sample_name == sample_name])[1],
    stringsAsFactors = FALSE
  )
}

counts159 <- do.call(cbind, mats)
meta159 <- do.call(rbind, meta_rows)
norm159 <- normalize_log1p(counts159)

state_names <- c("MACROPHAGE_INFLAMMATORY", "MACROPHAGE_FOAM_TREM2", "MACROPHAGE_C1Q", "MACROPHAGE_IFN")
state_scores <- sapply(state_names, function(state_name) {
  score_gene_set(norm159, state_sets$gene_symbol[state_sets$state_set == state_name])
})
if (is.vector(state_scores)) {
  state_scores <- matrix(state_scores, ncol = 1, dimnames = list(colnames(norm159), state_names))
}
meta159$dominant_state <- colnames(state_scores)[max.col(state_scores, ties.method = "first")]

cell_module_rows <- list()
for (module_name in unique(modules$module_name)) {
  genes <- modules$gene_symbol[modules$module_name == module_name]
  cell_module_rows[[module_name]] <- data.frame(
    cell_id = colnames(norm159),
    module_name = module_name,
    module_score = score_gene_set(norm159, genes),
    stringsAsFactors = FALSE
  )
}
cell_module_scores <- do.call(rbind, cell_module_rows)
cell_module_scores <- merge(cell_module_scores, meta159, by = "cell_id", all.x = TRUE)

state_summary <- aggregate(
  module_score ~ module_name + location + dominant_state,
  data = cell_module_scores,
  FUN = mean
)

sample_summary <- aggregate(
  module_score ~ module_name + sample_id + patient_id + location,
  data = cell_module_scores,
  FUN = mean
)

pair_rows <- list()
for (module_name in unique(sample_summary$module_name)) {
  sub <- sample_summary[sample_summary$module_name == module_name, ]
  paired <- merge(
    sub[sub$location == "core", c("patient_id", "module_score")],
    sub[sub$location == "adjacent", c("patient_id", "module_score")],
    by = "patient_id",
    suffixes = c("_core", "_adjacent")
  )
  paired$module_name <- module_name
  paired$core_minus_adjacent <- paired$module_score_core - paired$module_score_adjacent
  pair_rows[[module_name]] <- paired
}
pair_summary <- do.call(rbind, pair_rows)

## GSE155512 support layer
sample_suffix <- c(GSM4705589 = "RPE004", GSM4705590 = "RPE005", GSM4705591 = "RPE006")
support_rows <- list()
for (gsm in names(sample_suffix)) {
  file <- project_path("data", "raw", "single_cell", "GSE155512", paste0(gsm, "_", sample_suffix[[gsm]], "_matrix.txt.gz"))
  mat <- read.delim(gzfile(file), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gene_symbol <- make.unique(mat[[1]])
  expr <- as.matrix(mat[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- gene_symbol
  norm <- normalize_log1p_dense(expr)

  meta <- cell155512[cell155512$sample_id == gsm, , drop = FALSE]
  for (cell_type in unique(meta$predicted_cell_type)) {
    keep <- intersect(colnames(norm), meta$barcode[meta$predicted_cell_type == cell_type])
    if (length(keep) < 10) {
      next
    }
    for (module_name in unique(modules$module_name)) {
      genes <- modules$gene_symbol[modules$module_name == module_name]
      support_rows[[length(support_rows) + 1]] <- data.frame(
        sample_id = gsm,
        symptom = unique(meta$symptom)[1],
        cell_type = cell_type,
        module_name = module_name,
        module_score = mean(score_gene_set(norm[, keep, drop = FALSE], genes), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}
support_scores <- do.call(rbind, support_rows)
support_summary <- aggregate(module_score ~ module_name + symptom + cell_type, data = support_scores, FUN = mean)

## Bulk tissue layer
score_bulk_dataset <- function(dataset_id, comparison_a, comparison_b) {
  expr <- read_tsv_auto(project_path("data", "processed", "bulk_gene", paste0(dataset_id, "_gene_expr.tsv.gz")))
  pheno <- load_pheno(dataset_id)
  rownames(expr) <- expr$gene_symbol
  expr_mat <- as.matrix(expr[, !(colnames(expr) %in% c("gene_symbol", "feature_id")), drop = FALSE])
  storage.mode(expr_mat) <- "double"

  rows <- list()
  for (module_name in unique(modules$module_name)) {
    genes <- intersect(modules$gene_symbol[modules$module_name == module_name], rownames(expr_mat))
    if (length(genes) == 0) {
      next
    }
    sample_scores <- data.frame(
      sample_id = colnames(expr_mat),
      module_name = module_name,
      module_score = colMeans(expr_mat[genes, , drop = FALSE]),
      stringsAsFactors = FALSE
    )
    sample_scores <- merge(sample_scores, pheno, by = "sample_id", all.x = TRUE)
    rows[[length(rows) + 1]] <- sample_scores
  }
  out <- do.call(rbind, rows)
  out$dataset_id <- dataset_id
  out$comparison_a <- comparison_a
  out$comparison_b <- comparison_b
  out
}

bulk_scores <- rbind(
  score_bulk_dataset("GSE43292", "intact", "plaque"),
  score_bulk_dataset("GSE100927", "control", "atherosclerotic"),
  score_bulk_dataset("GSE28829", "early", "advanced")
)

bulk_group_summary <- aggregate(
  module_score ~ dataset_id + module_name + group_label,
  data = bulk_scores,
  FUN = mean
)

bulk_contrast_rows <- list()
for (dataset_id in unique(bulk_scores$dataset_id)) {
  sub_dataset <- bulk_scores[bulk_scores$dataset_id == dataset_id, ]
  comparison_a <- unique(sub_dataset$comparison_a)[1]
  comparison_b <- unique(sub_dataset$comparison_b)[1]
  for (module_name in unique(sub_dataset$module_name)) {
    sub <- sub_dataset[sub_dataset$module_name == module_name, ]
    group_a <- sub$module_score[sub$group_label == comparison_a]
    group_b <- sub$module_score[sub$group_label == comparison_b]
    delta <- mean(group_b, na.rm = TRUE) - mean(group_a, na.rm = TRUE)
    p_val <- tryCatch(wilcox.test(group_b, group_a, exact = FALSE)$p.value, error = function(e) NA_real_)
    bulk_contrast_rows[[length(bulk_contrast_rows) + 1]] <- data.frame(
      dataset_id = dataset_id,
      module_name = module_name,
      comparison_a = comparison_a,
      comparison_b = comparison_b,
      mean_a = mean(group_a, na.rm = TRUE),
      mean_b = mean(group_b, na.rm = TRUE),
      delta_b_minus_a = delta,
      wilcox_p = p_val,
      stringsAsFactors = FALSE
    )
  }
}
bulk_contrast_summary <- do.call(rbind, bulk_contrast_rows)

## Methylation module support
meth_sub <- meth[!(is.na(meth$gene_symbol) | meth$gene_symbol == ""), ]
meth_sub <- meth_sub[order(meth_sub$adj.P.Val), ]
meth_sub <- meth_sub[!duplicated(meth_sub$gene_symbol), c("gene_symbol", "delta_beta", "adj.P.Val")]
names(meth_sub)[2:3] <- c("delta_beta_meth", "adj_p_meth")

meth_rows <- list()
for (module_name in unique(modules$module_name)) {
  sub <- meth_sub[meth_sub$gene_symbol %in% modules$gene_symbol[modules$module_name == module_name], ]
  meth_rows[[length(meth_rows) + 1]] <- data.frame(
    module_name = module_name,
    n_genes_with_meth = nrow(sub),
    mean_delta_beta = mean(sub$delta_beta_meth, na.rm = TRUE),
    median_delta_beta = median(sub$delta_beta_meth, na.rm = TRUE),
    n_fdr_sig = sum(sub$adj_p_meth < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
meth_summary <- do.call(rbind, meth_rows)

write_tsv(state_summary, project_path("res", "tables", "mechanism", "gse159677_npl_module_state_scores.tsv"))
write_tsv(sample_summary, project_path("res", "tables", "mechanism", "gse159677_npl_module_sample_scores.tsv"))
write_tsv(pair_summary, project_path("res", "tables", "mechanism", "gse159677_npl_module_pair_summary.tsv"))
write_tsv(support_summary, project_path("res", "tables", "mechanism", "gse155512_npl_module_celltype_scores.tsv"))
write_tsv(bulk_scores, project_path("res", "tables", "mechanism", "npl_module_bulk_scores.tsv"))
write_tsv(bulk_group_summary, project_path("res", "tables", "mechanism", "npl_module_bulk_group_summary.tsv"))
write_tsv(bulk_contrast_summary, project_path("res", "tables", "mechanism", "npl_module_bulk_contrast_summary.tsv"))
write_tsv(meth_summary, project_path("res", "tables", "mechanism", "npl_module_methylation_support.tsv"))

cat("NPL module scores written to res/tables/mechanism\n")
print(state_summary)
cat("\nBulk contrast summary:\n")
print(bulk_contrast_summary)
cat("\nMethylation summary:\n")
print(meth_summary)
