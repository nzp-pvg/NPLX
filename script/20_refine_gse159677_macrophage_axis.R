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

normalize_log1p_dense <- function(mat) {
  lib <- colSums(mat)
  lib[lib == 0] <- 1
  log1p(t(t(mat) / lib * 10000))
}

infer_patient_id <- function(sample_id) {
  ifelse(
    sample_id %in% c("GSM4837523", "GSM4837524"),
    "P1",
    ifelse(sample_id %in% c("GSM4837525", "GSM4837526"), "P2", "P3")
  )
}

score_gse159677_macrophage_samples <- function(candidate_genes, cell_typing) {
  sample_meta <- unique(cell_typing[cell_typing$predicted_cell_type == "MACROPHAGE", c("sample_name", "sample_id", "location")])
  rows <- list()

  for (idx in seq_len(nrow(sample_meta))) {
    sample_name <- sample_meta$sample_name[[idx]]
    sample_dir <- project_path("data", "raw", "single_cell", "GSE159677", "per_sample", sample_name)
    counts <- read_10x_triplet(find_triplet_dir(sample_dir))
    keep <- cell_typing$barcode[cell_typing$sample_name == sample_name & cell_typing$predicted_cell_type == "MACROPHAGE"]
    keep <- intersect(colnames(counts), keep)
    if (length(keep) == 0) {
      next
    }

    norm <- normalize_log1p(counts[, keep, drop = FALSE])
    gene_keep <- intersect(candidate_genes, rownames(norm))
    if (length(gene_keep) == 0) {
      next
    }

    rows[[length(rows) + 1]] <- data.frame(
      sample_name = sample_name,
      sample_id = sample_meta$sample_id[[idx]],
      patient_id = infer_patient_id(sample_meta$sample_id[[idx]]),
      location = sample_meta$location[[idx]],
      gene_symbol = gene_keep,
      mean_expr = apply(norm[gene_keep, , drop = FALSE], 1, mean),
      n_cells = length(keep),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

score_gse155512_celltype_support <- function(candidate_genes, cell_typing) {
  sample_suffix <- c(
    GSM4705589 = "RPE004",
    GSM4705590 = "RPE005",
    GSM4705591 = "RPE006"
  )
  rows <- list()

  for (gsm in names(sample_suffix)) {
    file <- project_path("data", "raw", "single_cell", "GSE155512", paste0(gsm, "_", sample_suffix[[gsm]], "_matrix.txt.gz"))
    mat <- read.delim(gzfile(file), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    gene_symbol <- make.unique(mat[[1]])
    expr <- as.matrix(mat[, -1, drop = FALSE])
    storage.mode(expr) <- "double"
    rownames(expr) <- gene_symbol
    norm <- normalize_log1p_dense(expr)

    meta <- cell_typing[cell_typing$sample_id == gsm, , drop = FALSE]
    for (cell_type in unique(meta$predicted_cell_type)) {
      keep <- intersect(colnames(norm), meta$barcode[meta$predicted_cell_type == cell_type])
      if (length(keep) < 10) {
        next
      }
      gene_keep <- intersect(candidate_genes, rownames(norm))
      if (length(gene_keep) == 0) {
        next
      }
      rows[[length(rows) + 1]] <- data.frame(
        sample_id = gsm,
        symptom = unique(meta$symptom)[1],
        cell_type = cell_type,
        gene_symbol = gene_keep,
        mean_expr = apply(norm[gene_keep, keep, drop = FALSE], 1, mean),
        n_cells = length(keep),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

core_integrated <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_core_vs_adjacent_integrated.tsv"))
state_markers <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_state_markers.tsv"))
state_sample <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_state_sample_summary.tsv"))
cell159677 <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
cell155512 <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse155512_cell_level_typing.tsv.gz"))

foam_markers <- state_markers[state_markers$dominant_state == "MACROPHAGE_FOAM_TREM2", ]
foam_markers <- foam_markers[order(-foam_markers$delta_state_vs_other), ]
core_integrated <- core_integrated[order(-core_integrated$priority_score), ]

candidate_genes <- unique(c(
  head(core_integrated$gene_symbol, 120),
  head(foam_markers$gene_symbol, 120),
  c("NPL", "RNASET2", "CAPG", "VAMP8", "NPC2", "FABP5", "SPP1", "LGALS3", "HMOX1", "GPNMB", "PLA2G7", "C15orf48")
))
candidate_genes <- candidate_genes[!(is.na(candidate_genes) | candidate_genes == "")]

sample_gene_expr <- score_gse159677_macrophage_samples(candidate_genes, cell159677)
sample_gene_expr$location <- factor(sample_gene_expr$location, levels = c("adjacent", "core"))
sample_gene_expr <- sample_gene_expr[order(sample_gene_expr$gene_symbol, sample_gene_expr$patient_id, sample_gene_expr$location), ]

foam_fraction <- state_sample[state_sample$dominant_state == "MACROPHAGE_FOAM_TREM2", c("sample_id", "location", "fraction")]
names(foam_fraction)[3] <- "foam_fraction"
sample_meta <- unique(sample_gene_expr[, c("sample_id", "patient_id", "location")])
sample_meta <- merge(sample_meta, foam_fraction, by = c("sample_id", "location"), all.x = TRUE)
sample_meta$foam_fraction[is.na(sample_meta$foam_fraction)] <- 0
sample_gene_expr <- merge(sample_gene_expr, sample_meta, by = c("sample_id", "patient_id", "location"), all.x = TRUE)

gene_rows <- list()
pair_rows <- list()

for (gene in unique(sample_gene_expr$gene_symbol)) {
  sub <- sample_gene_expr[sample_gene_expr$gene_symbol == gene, ]
  core_vals <- sub$mean_expr[sub$location == "core"]
  adj_vals <- sub$mean_expr[sub$location == "adjacent"]
  paired <- merge(
    sub[sub$location == "core", c("patient_id", "mean_expr", "foam_fraction")],
    sub[sub$location == "adjacent", c("patient_id", "mean_expr", "foam_fraction")],
    by = "patient_id",
    suffixes = c("_core", "_adj")
  )
  paired$delta_core_minus_adjacent <- paired$mean_expr_core - paired$mean_expr_adj
  pair_rows[[length(pair_rows) + 1]] <- data.frame(
    gene_symbol = gene,
    paired,
    stringsAsFactors = FALSE
  )

  foam_rho <- suppressWarnings(cor(sub$mean_expr, sub$foam_fraction, method = "spearman"))
  gene_rows[[length(gene_rows) + 1]] <- data.frame(
    gene_symbol = gene,
    sample_mean_core = mean(core_vals, na.rm = TRUE),
    sample_mean_adjacent = mean(adj_vals, na.rm = TRUE),
    paired_mean_delta = mean(paired$delta_core_minus_adjacent, na.rm = TRUE),
    paired_median_delta = median(paired$delta_core_minus_adjacent, na.rm = TRUE),
    n_patients_core_gt_adj = sum(paired$delta_core_minus_adjacent > 0, na.rm = TRUE),
    min_patient_delta = min(paired$delta_core_minus_adjacent, na.rm = TRUE),
    foam_fraction_spearman = foam_rho,
    stringsAsFactors = FALSE
  )
}

gene_summary <- do.call(rbind, gene_rows)
pair_summary <- do.call(rbind, pair_rows)

support_celltype <- score_gse155512_celltype_support(candidate_genes, cell155512)
support_asym <- support_celltype[support_celltype$symptom == "asymptomatic", ]
support_group <- aggregate(mean_expr ~ cell_type + gene_symbol, data = support_asym, FUN = mean)

support_gene_rows <- list()
for (gene in unique(support_group$gene_symbol)) {
  sub <- support_group[support_group$gene_symbol == gene, ]
  sub <- sub[order(-sub$mean_expr), ]
  macrophage_mean <- sub$mean_expr[sub$cell_type == "MACROPHAGE"]
  macrophage_mean <- if (length(macrophage_mean) == 0) NA_real_ else macrophage_mean[[1]]
  macrophage_rank <- match("MACROPHAGE", sub$cell_type)
  non_mac <- sub$mean_expr[sub$cell_type != "MACROPHAGE"]
  second_mean <- if (length(non_mac) == 0) NA_real_ else max(non_mac)
  support_gene_rows[[length(support_gene_rows) + 1]] <- data.frame(
    gene_symbol = gene,
    external_top_cell_type = sub$cell_type[[1]],
    external_macrophage_mean = macrophage_mean,
    external_second_mean = second_mean,
    external_macrophage_rank = ifelse(is.na(macrophage_rank), NA_integer_, macrophage_rank),
    external_macrophage_specificity_ratio = ifelse(
      is.na(macrophage_mean) || is.na(second_mean),
      NA_real_,
      (macrophage_mean + 1e-06) / (second_mean + 1e-06)
    ),
    stringsAsFactors = FALSE
  )
}

support_gene_summary <- do.call(rbind, support_gene_rows)

foam_sub <- foam_markers[, c("gene_symbol", "mean_in_state", "mean_out_state", "delta_state_vs_other")]
names(foam_sub)[2:4] <- c("foam_mean_expr", "nonfoam_mean_expr", "foam_delta_vs_other")

refined <- merge(core_integrated, foam_sub, by = "gene_symbol", all.x = TRUE)
refined <- merge(refined, gene_summary, by = "gene_symbol", all.x = TRUE)
refined <- merge(refined, support_gene_summary, by = "gene_symbol", all.x = TRUE)
refined$bulk_integrated_priority_score <- refined$integrated_priority_score

refined$patient_consistency_fraction <- refined$n_patients_core_gt_adj / 3
refined$bulk_direction_consistency <- rowSums(as.matrix(cbind(
  refined$logFC_GSE43292 > 0,
  refined$logFC_GSE100927 > 0,
  refined$logFC_GSE28829 > 0
)), na.rm = TRUE)
refined$external_macrophage_specificity_ratio[!is.finite(refined$external_macrophage_specificity_ratio)] <- NA_real_

core_weight <- (-log10(pmax(refined$adj_p_value, 1e-300))) * pmax(refined$delta_core_minus_adjacent, 0.05)
foam_weight <- pmax(refined$foam_delta_vs_other, 0.05)
bulk_weight <- log10(pmax(refined$bulk_integrated_priority_score, 1) + 1)
meth_weight <- ifelse(
  is.na(refined$adj_p_meth),
  0.5,
  (-log10(pmax(refined$adj_p_meth, 1e-300))) * pmax(abs(refined$delta_beta_meth), 0.02)
)
patient_weight <- pmax(refined$patient_consistency_fraction, 1 / 3)
external_weight <- ifelse(
  is.na(refined$external_macrophage_rank),
  1,
  ifelse(refined$external_macrophage_rank == 1, 1.5, ifelse(refined$external_macrophage_rank == 2, 1.2, 0.9))
)

refined$refined_macrophage_priority_score <- core_weight * foam_weight * bulk_weight * meth_weight * patient_weight * external_weight
refined <- refined[order(-refined$refined_macrophage_priority_score, -refined$foam_delta_vs_other, -refined$paired_mean_delta), ]

top_cols <- c(
  "gene_symbol",
  "delta_core_minus_adjacent",
  "adj_p_value",
  "foam_delta_vs_other",
  "paired_mean_delta",
  "n_patients_core_gt_adj",
  "foam_fraction_spearman",
  "external_top_cell_type",
  "external_macrophage_specificity_ratio",
  "bulk_integrated_priority_score",
  "adj_p_meth",
  "refined_macrophage_priority_score"
)
top_table <- head(refined[, top_cols], 40)

write_tsv(sample_gene_expr, project_path("res", "tables", "mechanism", "gse159677_macrophage_sample_gene_expression.tsv"))
write_tsv(pair_summary, project_path("res", "tables", "mechanism", "gse159677_macrophage_patient_pair_gene_summary.tsv"))
write_tsv(support_celltype, project_path("res", "tables", "mechanism", "gse155512_macrophage_candidate_celltype_support.tsv"))
write_tsv(support_gene_summary, project_path("res", "tables", "mechanism", "gse155512_macrophage_candidate_gene_support_summary.tsv"))
write_tsv(refined, project_path("res", "tables", "mechanism", "gse159677_macrophage_refined_candidates.tsv"))
write_tsv(top_table, project_path("res", "qc", "mechanism", "gse159677_macrophage_refined_top40.tsv"))

cat("Refined macrophage-axis outputs written to res/tables/mechanism and res/qc/mechanism\n")
print(top_table)
