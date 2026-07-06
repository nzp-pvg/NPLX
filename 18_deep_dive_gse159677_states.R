source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(Matrix))

ensure_project_dirs()

args <- commandArgs(trailingOnly = TRUE)
target_cell_type <- if (length(args) >= 1) toupper(args[[1]]) else "MACROPHAGE"
if (!target_cell_type %in% c("MACROPHAGE", "SMC")) {
  stop("Target cell type must be MACROPHAGE or SMC")
}

pheno <- load_pheno("GSE159677")
cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))
bulk_cons <- read_tsv_auto(project_path("res", "tables", "bulk", "athero_bulk_discovery_consensus.tsv.gz"))
meth <- read_tsv_auto(project_path("res", "tables", "mechanism", "GSE46394_methylation_limma.tsv.gz"))

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

build_celltype_object <- function(target_cell_type) {
  sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
  mats <- list()
  meta_rows <- list()

  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    if (!grepl("^GSM", sample_name)) {
      next
    }
    gsm <- sub("_.*$", "", sample_name)
    triplet_dir <- find_triplet_dir(sample_dir)
    counts <- read_10x_triplet(triplet_dir)
    typing_sub <- cell_typing[cell_typing$sample_name == sample_name & cell_typing$predicted_cell_type == target_cell_type, , drop = FALSE]
    if (nrow(typing_sub) == 0) {
      next
    }
    keep_barcodes <- intersect(colnames(counts), typing_sub$barcode)
    if (length(keep_barcodes) == 0) {
      next
    }
    counts <- counts[, keep_barcodes, drop = FALSE]
    colnames(counts) <- paste(sample_name, keep_barcodes, sep = "|")
    mats[[sample_name]] <- counts

    location <- pheno$group_label[match(gsm, pheno$sample_id)]
    meta_rows[[sample_name]] <- data.frame(
      cell_id = colnames(counts),
      sample_name = sample_name,
      sample_id = gsm,
      location = location,
      stringsAsFactors = FALSE
    )
  }

  list(counts = do.call(cbind, mats), meta = do.call(rbind, meta_rows))
}

analyze_celltype_states <- function(target_cell_type) {
  obj <- build_celltype_object(target_cell_type)
  counts <- obj$counts
  meta <- obj$meta
  norm <- normalize_log1p(counts)

  relevant_states <- if (target_cell_type == "MACROPHAGE") {
    c("MACROPHAGE_INFLAMMATORY", "MACROPHAGE_FOAM_TREM2", "MACROPHAGE_C1Q", "MACROPHAGE_IFN")
  } else {
    c("SMC_CONTRACTILE", "SMC_FIBROMYOCYTE", "SMC_OSTEO_STRESS", "SMC_INFLAMMATORY")
  }

  score_mat <- sapply(relevant_states, function(set_name) {
    genes <- state_sets$gene_symbol[state_sets$state_set == set_name]
    score_gene_set(norm, genes)
  })
  if (is.vector(score_mat)) {
    score_mat <- matrix(score_mat, ncol = 1, dimnames = list(colnames(norm), relevant_states))
  }

  meta$dominant_state <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
  meta$dominant_score <- apply(score_mat, 1, max, na.rm = TRUE)
  meta$state_margin <- apply(score_mat, 1, function(x) {
    x <- sort(as.numeric(x), decreasing = TRUE)
    if (length(x) < 2) return(NA_real_)
    x[1] - x[2]
  })

  state_sample_summary <- aggregate(
    rep(1, nrow(meta)),
    by = list(sample_id = meta$sample_id, location = meta$location, dominant_state = meta$dominant_state),
    FUN = sum
  )
  names(state_sample_summary)[4] <- "n_cells"
  totals <- aggregate(n_cells ~ sample_id + location, data = state_sample_summary, FUN = sum)
  state_sample_summary <- merge(state_sample_summary, totals, by = c("sample_id", "location"), suffixes = c("", "_total"))
  state_sample_summary$fraction <- state_sample_summary$n_cells / state_sample_summary$n_cells_total

  state_group_summary <- aggregate(
    cbind(n_cells, fraction) ~ location + dominant_state,
    data = state_sample_summary,
    FUN = mean
  )

  state_score_summary <- aggregate(
    score_mat,
    by = list(location = meta$location, dominant_state = meta$dominant_state),
    FUN = mean
  )

  selected_genes <- unique(c(
    state_sets$gene_symbol,
    bulk_cons$gene_symbol[1:min(800, nrow(bulk_cons))]
  ))
  selected_genes <- intersect(selected_genes, rownames(norm))

  location_factor <- factor(meta$location, levels = c("adjacent", "core"))
  de_rows <- list()
  for (gene in selected_genes) {
    vals <- as.numeric(norm[gene, ])
    if (var(vals) == 0) {
      next
    }
    wt <- tryCatch(wilcox.test(vals[location_factor == "core"], vals[location_factor == "adjacent"], exact = FALSE), error = function(e) NULL)
    if (is.null(wt)) {
      next
    }
    de_rows[[length(de_rows) + 1]] <- data.frame(
      gene_symbol = gene,
      mean_core = mean(vals[location_factor == "core"]),
      mean_adjacent = mean(vals[location_factor == "adjacent"]),
      delta_core_minus_adjacent = mean(vals[location_factor == "core"]) - mean(vals[location_factor == "adjacent"]),
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
  de_table <- do.call(rbind, de_rows)
  de_table$adj_p_value <- p.adjust(de_table$p_value, method = "BH")

  # State-specific markers using dominant-state partitions.
  marker_rows <- list()
  for (state_name in unique(meta$dominant_state)) {
    in_state <- meta$dominant_state == state_name
    out_state <- meta$dominant_state != state_name
    for (gene in selected_genes) {
      vals <- as.numeric(norm[gene, ])
      marker_rows[[length(marker_rows) + 1]] <- data.frame(
        dominant_state = state_name,
        gene_symbol = gene,
        mean_in_state = mean(vals[in_state]),
        mean_out_state = mean(vals[out_state]),
        delta_state_vs_other = mean(vals[in_state]) - mean(vals[out_state]),
        stringsAsFactors = FALSE
      )
    }
  }
  marker_table <- do.call(rbind, marker_rows)
  marker_table <- marker_table[order(marker_table$dominant_state, -marker_table$delta_state_vs_other), ]

  meth_sub <- meth[!(is.na(meth$gene_symbol) | meth$gene_symbol == ""), ]
  meth_sub <- meth_sub[order(meth_sub$adj.P.Val), ]
  meth_sub <- meth_sub[!duplicated(meth_sub$gene_symbol), c("gene_symbol", "delta_beta", "adj.P.Val")]
  names(meth_sub)[-1] <- c("delta_beta_meth", "adj_p_meth")
  bulk_sub <- bulk_cons[, c("gene_symbol", "integrated_priority_score", "discovery_fisher_fdr", "logFC_GSE43292", "logFC_GSE100927", "logFC_GSE28829")]

  integrated <- merge(de_table, bulk_sub, by = "gene_symbol", all.x = TRUE)
  integrated <- merge(integrated, meth_sub, by = "gene_symbol", all.x = TRUE)
  integrated$priority_score <- with(
    integrated,
    (-log10(pmax(adj_p_value, 1e-300))) *
      pmax(abs(delta_core_minus_adjacent), 0.01) *
      ifelse(!is.na(integrated_priority_score), pmax(integrated_priority_score, 0.1), 0.1)
  )
  integrated <- integrated[order(-integrated$priority_score), ]

  list(
    meta = meta,
    state_sample_summary = state_sample_summary,
    state_group_summary = state_group_summary,
    state_score_summary = state_score_summary,
    marker_table = marker_table,
    integrated = integrated
  )
}

res <- analyze_celltype_states(target_cell_type)
prefix <- tolower(target_cell_type)

write_tsv(res$state_sample_summary, project_path("res", "tables", "mechanism", paste0("gse159677_", prefix, "_state_sample_summary.tsv")))
write_tsv(res$state_group_summary, project_path("res", "tables", "mechanism", paste0("gse159677_", prefix, "_state_group_summary.tsv")))
write_tsv(res$state_score_summary, project_path("res", "tables", "mechanism", paste0("gse159677_", prefix, "_state_score_summary.tsv")))
write_tsv(res$marker_table, project_path("res", "tables", "mechanism", paste0("gse159677_", prefix, "_state_markers.tsv")))
write_tsv(res$integrated, project_path("res", "tables", "mechanism", paste0("gse159677_", prefix, "_core_vs_adjacent_integrated.tsv")))

top_tbl <- head(res$integrated[, c("gene_symbol", "delta_core_minus_adjacent", "adj_p_value", "integrated_priority_score", "adj_p_meth")], 30)
write_tsv(top_tbl, project_path("res", "qc", "mechanism", paste0("gse159677_", prefix, "_top30.tsv")))

cat("GSE159677 state analysis completed for ", target_cell_type, "\n", sep = "")
print(res$state_group_summary)
cat("\nTop candidates:\n")
print(top_tbl)
