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

state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))
cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
core_integrated <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_core_vs_adjacent_integrated.tsv"))

sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
mats <- list()
meta_rows <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  if (!grepl("^GSM", sample_name)) {
    next
  }

  counts <- read_10x_triplet(find_triplet_dir(sample_dir))
  keep <- cell_typing$barcode[cell_typing$sample_name == sample_name & cell_typing$predicted_cell_type == "MACROPHAGE"]
  keep <- intersect(colnames(counts), keep)
  if (length(keep) == 0) {
    next
  }

  counts <- counts[, keep, drop = FALSE]
  colnames(counts) <- paste(sample_name, keep, sep = "|")
  mats[[sample_name]] <- counts
  meta_rows[[sample_name]] <- data.frame(
    cell_id = colnames(counts),
    sample_name = sample_name,
    location = unique(cell_typing$location[cell_typing$sample_name == sample_name])[1],
    stringsAsFactors = FALSE
  )
}

counts <- do.call(cbind, mats)
meta <- do.call(rbind, meta_rows)
norm <- normalize_log1p(counts)

foam_score <- score_gene_set(norm, state_sets$gene_symbol[state_sets$state_set == "MACROPHAGE_FOAM_TREM2"])
c1q_score <- score_gene_set(norm, state_sets$gene_symbol[state_sets$state_set == "MACROPHAGE_C1Q"])
transition_axis <- foam_score - c1q_score

core_keep <- meta$location == "core"
core_axis <- transition_axis[core_keep]
candidate_genes <- intersect(core_integrated$gene_symbol, rownames(norm))

axis_rows <- lapply(candidate_genes, function(gene_symbol) {
  data.frame(
    gene_symbol = gene_symbol,
    rho_core_axis = suppressWarnings(cor(as.numeric(norm[gene_symbol, core_keep]), core_axis, method = "spearman")),
    stringsAsFactors = FALSE
  )
})

axis_table <- do.call(rbind, axis_rows)
axis_table <- axis_table[order(-axis_table$rho_core_axis), ]
axis_table$rank_core_axis <- seq_len(nrow(axis_table))

candidate_subset <- axis_table[axis_table$gene_symbol %in% c(
  "NPL", "RNASET2", "CAPG", "VAMP8", "NPC2", "FABP5", "SPP1",
  "LGALS3", "HMOX1", "GPNMB", "PLA2G7", "C15orf48", "APOE", "APOC1"
), ]

write_tsv(axis_table, project_path("res", "tables", "mechanism", "gse159677_macrophage_core_axis_allgenes.tsv"))
write_tsv(candidate_subset, project_path("res", "qc", "mechanism", "gse159677_macrophage_core_axis_candidate_cor.tsv"))

cat("Macrophage transition-axis correlation tables written to res/tables/mechanism and res/qc/mechanism\n")
print(head(axis_table, 40))
cat("\nCandidate subset:\n")
print(candidate_subset)
