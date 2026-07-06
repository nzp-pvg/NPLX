source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(limma)
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

infer_patient_id <- function(sample_id) {
  ifelse(
    sample_id %in% c("GSM4837523", "GSM4837524"),
    "P1",
    ifelse(sample_id %in% c("GSM4837525", "GSM4837526"), "P2", "P3")
  )
}

build_pseudobulk <- function() {
  pheno <- load_pheno("GSE159677")
  cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
  sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)

  mats <- list()
  meta_rows <- list()

  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    if (!grepl("^GSM", sample_name)) {
      next
    }
    sample_id <- sub("_.*$", "", sample_name)
    counts <- read_10x_triplet(find_triplet_dir(sample_dir))
    keep <- cell_typing$barcode[cell_typing$sample_name == sample_name & cell_typing$predicted_cell_type == "MACROPHAGE"]
    keep <- intersect(colnames(counts), keep)
    if (length(keep) == 0) {
      next
    }
    counts <- counts[, keep, drop = FALSE]
    mats[[sample_id]] <- Matrix::rowSums(counts)
    meta_rows[[sample_id]] <- data.frame(
      sample_id = sample_id,
      patient_id = infer_patient_id(sample_id),
      location = pheno$group_label[match(sample_id, pheno$sample_id)],
      n_cells = length(keep),
      stringsAsFactors = FALSE
    )
  }

  mat <- do.call(cbind, mats)
  meta <- do.call(rbind, meta_rows)
  colnames(mat) <- meta$sample_id
  list(counts = mat, meta = meta)
}

obj <- build_pseudobulk()
counts <- obj$counts
meta <- obj$meta[match(colnames(counts), obj$meta$sample_id), , drop = FALSE]

lib_size <- colSums(counts)
logcpm <- log2(t(t(counts) / pmax(lib_size, 1) * 1e6) + 1)

patient <- factor(meta$patient_id, levels = c("P1", "P2", "P3"))
location <- factor(meta$location, levels = c("adjacent", "core"))
design <- model.matrix(~ patient + location)

fit <- eBayes(lmFit(logcpm, design))
tt <- topTable(fit, coef = "locationcore", number = Inf, sort.by = "P")
tt$gene_symbol <- rownames(tt)
tt$n_pairs_core_gt_adj <- NA_integer_
tt$mean_pair_delta <- NA_real_
tt$median_pair_delta <- NA_real_

paired_delta_rows <- list()
for (gene_symbol in rownames(logcpm)) {
  vals <- as.numeric(logcpm[gene_symbol, ])
  tmp <- data.frame(
    patient_id = meta$patient_id,
    location = meta$location,
    logcpm = vals,
    stringsAsFactors = FALSE
  )
  wide <- reshape(tmp, idvar = "patient_id", timevar = "location", direction = "wide")
  wide$delta_core_minus_adj <- wide$logcpm.core - wide$logcpm.adjacent
  tt[tt$gene_symbol == gene_symbol, "n_pairs_core_gt_adj"] <- sum(wide$delta_core_minus_adj > 0, na.rm = TRUE)
  tt[tt$gene_symbol == gene_symbol, "mean_pair_delta"] <- mean(wide$delta_core_minus_adj, na.rm = TRUE)
  tt[tt$gene_symbol == gene_symbol, "median_pair_delta"] <- median(wide$delta_core_minus_adj, na.rm = TRUE)
  paired_delta_rows[[gene_symbol]] <- data.frame(gene_symbol = gene_symbol, wide, stringsAsFactors = FALSE)
}

paired_delta_table <- do.call(rbind, paired_delta_rows)
tt <- tt[order(tt$P.Value, -tt$mean_pair_delta), ]

modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
module_rows <- list()
for (module_name in unique(modules$module_name)) {
  genes <- intersect(modules$gene_symbol[modules$module_name == module_name], rownames(logcpm))
  if (length(genes) == 0) {
    next
  }
  score <- colMeans(logcpm[genes, , drop = FALSE])
  module_rows[[module_name]] <- data.frame(
    sample_id = names(score),
    module_name = module_name,
    module_score = as.numeric(score),
    stringsAsFactors = FALSE
  )
}
module_scores <- do.call(rbind, module_rows)
module_scores <- merge(module_scores, meta, by = "sample_id", all.x = TRUE)

module_pair <- reshape(
  module_scores[, c("patient_id", "location", "module_name", "module_score")],
  idvar = c("patient_id", "module_name"),
  timevar = "location",
  direction = "wide"
)
module_pair$delta_core_minus_adjacent <- module_pair$module_score.core - module_pair$module_score.adjacent

module_summary <- do.call(
  rbind,
  lapply(split(module_pair$delta_core_minus_adjacent, module_pair$module_name), function(x) {
    data.frame(
      mean_delta = mean(x, na.rm = TRUE),
      median_delta = median(x, na.rm = TRUE),
      n_pairs_core_gt_adj = sum(x > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
module_summary$module_name <- rownames(module_summary)
rownames(module_summary) <- NULL
module_summary <- module_summary[, c("module_name", "mean_delta", "median_delta", "n_pairs_core_gt_adj")]

top_genes <- c("NPL", "FABP5", "GPNMB", "APOC1", "SPP1", "CD36", "APOE", "PLA2G7", "LGALS3", "HMOX1", "CYP27A1")
gene_plot_df <- data.frame(
  sample_id = rep(colnames(logcpm), each = length(top_genes)),
  gene_symbol = rep(top_genes, times = ncol(logcpm)),
  logcpm = as.numeric(logcpm[top_genes, , drop = FALSE]),
  stringsAsFactors = FALSE
)
gene_plot_df <- merge(gene_plot_df, meta, by = "sample_id", all.x = TRUE)

write_tsv(tt, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_de.tsv"))
write_tsv(paired_delta_table, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_pair_deltas.tsv"))
write_tsv(meta, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_sample_meta.tsv"))
write_tsv(module_scores, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_module_scores.tsv"))
write_tsv(module_pair, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudobulk_module_pairs.tsv"))
write_tsv(module_summary, project_path("res", "qc", "mechanism", "gse159677_macrophage_pseudobulk_module_summary.tsv"))
write_tsv(head(tt[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "mean_pair_delta", "n_pairs_core_gt_adj")], 40),
          project_path("res", "qc", "mechanism", "gse159677_macrophage_pseudobulk_top40.tsv"))

p_module <- ggplot(module_scores, aes(location, module_score, group = patient_id, color = patient_id)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ module_name, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(title = "Macrophage pseudobulk module shifts", x = NULL, y = "Mean logCPM")

p_gene <- ggplot(gene_plot_df, aes(location, logcpm, group = patient_id, color = patient_id)) +
  geom_point(size = 1.5) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ gene_symbol, scales = "free_y") +
  theme_bw(base_size = 9) +
  labs(title = "Selected macrophage genes in pseudobulk pairs", x = NULL, y = "logCPM")

ggsave(project_path("figure", "export", "mechanism", "gse159677_macrophage_pseudobulk_modules.pdf"), p_module, width = 8.5, height = 5.5)
ggsave(project_path("figure", "export", "mechanism", "gse159677_macrophage_pseudobulk_selected_genes.pdf"), p_gene, width = 9.5, height = 7.0)

cat("Macrophage pseudobulk outputs written to res/tables/mechanism and figure/export/mechanism\n")
print(head(tt[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val", "mean_pair_delta", "n_pairs_core_gt_adj")], 20))
print(module_summary)
