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

cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))

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

core_keep <- meta$location == "core"
core_expr <- norm[, core_keep, drop = FALSE]
core_npl <- as.numeric(core_expr["NPL", ])

lipid_panel <- c("SPP1", "FABP5", "APOE", "APOC1", "GPNMB", "LGALS3", "CYP27A1", "LPL", "CD36", "PLA2G7")
inflammatory_panel <- c("HMOX1", "IL1RN", "CCL2", "NFKBIA", "FOS", "JUNB", "PTGS2", "TNF")

mean_expr_rank <- order(Matrix::rowMeans(core_expr), decreasing = TRUE)
candidate_genes <- unique(c(
  rownames(core_expr)[mean_expr_rank[seq_len(min(500, length(mean_expr_rank)))]],
  lipid_panel,
  inflammatory_panel,
  "NPL",
  "RNASET2",
  "CAPG",
  "VAMP8"
))
candidate_genes <- intersect(candidate_genes, rownames(core_expr))

cor_rows <- lapply(candidate_genes, function(gene_symbol) {
  data.frame(
    gene_symbol = gene_symbol,
    rho_with_npl = suppressWarnings(cor(core_npl, as.numeric(core_expr[gene_symbol, ]), method = "spearman")),
    mean_core_expr = mean(core_expr[gene_symbol, ]),
    stringsAsFactors = FALSE
  )
})

cor_table <- do.call(rbind, cor_rows)
cor_table <- cor_table[order(-cor_table$rho_with_npl), ]
cor_table$rank_with_npl <- seq_len(nrow(cor_table))

panel_rows <- list(
  data.frame(
    panel = "lipid_handling_foam",
    mean_rho = mean(cor_table$rho_with_npl[cor_table$gene_symbol %in% lipid_panel], na.rm = TRUE),
    median_rho = median(cor_table$rho_with_npl[cor_table$gene_symbol %in% lipid_panel], na.rm = TRUE),
    n_genes = sum(cor_table$gene_symbol %in% lipid_panel),
    stringsAsFactors = FALSE
  ),
  data.frame(
    panel = "inflammatory_stress",
    mean_rho = mean(cor_table$rho_with_npl[cor_table$gene_symbol %in% inflammatory_panel], na.rm = TRUE),
    median_rho = median(cor_table$rho_with_npl[cor_table$gene_symbol %in% inflammatory_panel], na.rm = TRUE),
    n_genes = sum(cor_table$gene_symbol %in% inflammatory_panel),
    stringsAsFactors = FALSE
  )
)
panel_summary <- do.call(rbind, panel_rows)

candidate_subset <- cor_table[cor_table$gene_symbol %in% unique(c(lipid_panel, inflammatory_panel, "NPL", "RNASET2", "CAPG", "VAMP8")), ]

write_tsv(cor_table, project_path("res", "tables", "mechanism", "gse159677_npl_core_macrophage_correlations.tsv"))
write_tsv(head(cor_table, 100), project_path("res", "qc", "mechanism", "gse159677_npl_core_macrophage_top100.tsv"))
write_tsv(candidate_subset, project_path("res", "qc", "mechanism", "gse159677_npl_candidate_correlations.tsv"))
write_tsv(panel_summary, project_path("res", "qc", "mechanism", "gse159677_npl_panel_summary.tsv"))

cat("NPL core-macrophage correlation outputs written to res/tables/mechanism and res/qc/mechanism\n")
print(head(cor_table, 40))
cat("\nPanel summary:\n")
print(panel_summary)
