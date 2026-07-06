source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
})

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

pheno <- load_pheno("GSE159677")
cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
mechanism_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "athero_mechanism_gene_sets.tsv"))
modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))

sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
rows <- list()

for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  if (!grepl("^GSM", sample_name)) {
    next
  }
  sample_id <- sub("_.*$", "", sample_name)
  counts <- read_10x_triplet(find_triplet_dir(sample_dir))
  norm <- normalize_log1p(counts)
  sample_meta <- cell_typing[cell_typing$sample_name == sample_name, , drop = FALSE]
  location <- pheno$group_label[match(sample_id, pheno$sample_id)]

  for (cell_type in sort(unique(sample_meta$predicted_cell_type))) {
    keep <- intersect(colnames(norm), sample_meta$barcode[sample_meta$predicted_cell_type == cell_type])
    if (length(keep) < 20) {
      next
    }

    sub <- norm[, keep, drop = FALSE]
    all_sets <- unique(c(mechanism_sets$gene_set, modules$module_name))
    for (set_name in all_sets) {
      genes <- if (set_name %in% mechanism_sets$gene_set) {
        mechanism_sets$gene_symbol[mechanism_sets$gene_set == set_name]
      } else {
        modules$gene_symbol[modules$module_name == set_name]
      }
      score <- score_gene_set(sub, genes)
      rows[[length(rows) + 1]] <- data.frame(
        sample_name = sample_name,
        sample_id = sample_id,
        location = location,
        cell_type = cell_type,
        gene_set = set_name,
        mean_score = mean(score, na.rm = TRUE),
        median_score = median(score, na.rm = TRUE),
        n_cells = length(keep),
        stringsAsFactors = FALSE
      )
    }
  }
}

score_table <- do.call(rbind, rows)
group_table <- aggregate(
  cbind(mean_score, median_score, n_cells) ~ location + cell_type + gene_set,
  data = score_table,
  FUN = mean
)

write_tsv(score_table, project_path("res", "tables", "mechanism", "gse159677_celltype_restricted_program_scores.tsv"))
write_tsv(group_table, project_path("res", "tables", "mechanism", "gse159677_celltype_restricted_program_group_scores.tsv"))

focus_sets <- c("KEAP1_NRF2_RESPONSE", "ENDOTHELIAL_ACTIVATION_STRESS", "SMC_PHENOTYPE_SWITCH_STRESS", "NPL_FOAM_MACROPHAGE_COMPACT")
plot_df <- score_table[score_table$gene_set %in% focus_sets, , drop = FALSE]

p <- ggplot(plot_df, aes(location, mean_score, fill = location)) +
  geom_boxplot(outlier.size = 0.25) +
  facet_grid(gene_set ~ cell_type, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Cell type-restricted program scores in GSE159677", x = NULL, y = "Mean per-cell score")

ggsave(project_path("figure", "export", "mechanism", "gse159677_celltype_restricted_program_scores.pdf"), p, width = 12, height = 7.5)

cat("Cell type-restricted program scores written to res/tables/mechanism and figure/export/mechanism\n")
print(group_table[group_table$gene_set %in% focus_sets, ])
