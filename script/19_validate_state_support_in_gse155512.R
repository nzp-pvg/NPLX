source("script/R/00_project_config.R")

ensure_project_dirs()

state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))
cell_table <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse155512_cell_level_typing.tsv.gz"))
pheno <- load_pheno("GSE155512")

score_gene_set_dense <- function(norm_mat, genes) {
  genes <- intersect(genes, rownames(norm_mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(norm_mat)))
  }
  colMeans(norm_mat[genes, , drop = FALSE])
}

normalize_log1p_dense <- function(mat) {
  lib <- colSums(mat)
  lib[lib == 0] <- 1
  log1p(t(t(mat) / lib * 10000))
}

rows <- list()

for (gsm in unique(pheno$sample_id)) {
  file <- project_path("data", "raw", "single_cell", "GSE155512", paste0(gsm, "_", ifelse(gsm == "GSM4705589", "RPE004", ifelse(gsm == "GSM4705590", "RPE005", "RPE006")), "_matrix.txt.gz"))
  mat <- read.delim(gzfile(file), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gene_symbol <- mat[[1]]
  expr <- as.matrix(mat[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- make.unique(gene_symbol)
  norm <- normalize_log1p_dense(expr)

  cell_sub <- cell_table[cell_table$sample_id == gsm, , drop = FALSE]
  symptom <- unique(cell_sub$symptom)[1]

  for (target_cell_type in c("MACROPHAGE", "SMC")) {
    keep <- intersect(colnames(norm), cell_sub$barcode[cell_sub$predicted_cell_type == target_cell_type])
    if (length(keep) < 20) {
      next
    }
    target_norm <- norm[, keep, drop = FALSE]
    relevant <- if (target_cell_type == "MACROPHAGE") {
      c("MACROPHAGE_INFLAMMATORY", "MACROPHAGE_FOAM_TREM2", "MACROPHAGE_C1Q", "MACROPHAGE_IFN")
    } else {
      c("SMC_CONTRACTILE", "SMC_FIBROMYOCYTE", "SMC_OSTEO_STRESS", "SMC_INFLAMMATORY")
    }
    for (state_name in relevant) {
      genes <- state_sets$gene_symbol[state_sets$state_set == state_name]
      rows[[length(rows) + 1]] <- data.frame(
        sample_id = gsm,
        symptom = symptom,
        target_cell_type = target_cell_type,
        state_set = state_name,
        mean_score = mean(score_gene_set_dense(target_norm, genes), na.rm = TRUE),
        n_cells = length(keep),
        stringsAsFactors = FALSE
      )
    }
  }
}

summary_table <- do.call(rbind, rows)
group_summary <- aggregate(mean_score ~ symptom + target_cell_type + state_set, data = summary_table, FUN = mean)

write_tsv(summary_table, project_path("res", "tables", "mechanism", "gse155512_state_support_sample_scores.tsv"))
write_tsv(group_summary, project_path("res", "tables", "mechanism", "gse155512_state_support_group_scores.tsv"))

cat("GSE155512 state-support summaries written to res/tables/mechanism\n")
print(group_summary)
