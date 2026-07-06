source("script/R/00_project_config.R")
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scTenifoldKnk))

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

build_macrophage_counts <- function(cell_typing, keep_location = NULL) {
  sample_dirs <- list.dirs(project_path("data", "raw", "single_cell", "GSE159677", "per_sample"), recursive = FALSE, full.names = TRUE)
  mats <- list()
  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    if (!grepl("^GSM", sample_name)) {
      next
    }
    counts <- read_10x_triplet(find_triplet_dir(sample_dir))
    keep <- cell_typing$barcode[cell_typing$sample_name == sample_name & cell_typing$predicted_cell_type == "MACROPHAGE"]
    if (!is.null(keep_location)) {
      keep <- intersect(keep, cell_typing$barcode[cell_typing$sample_name == sample_name & cell_typing$location == keep_location])
    }
    keep <- intersect(colnames(counts), keep)
    if (length(keep) == 0) {
      next
    }
    counts <- counts[, keep, drop = FALSE]
    colnames(counts) <- paste(sample_name, keep, sep = "|")
    mats[[sample_name]] <- counts
  }
  do.call(cbind, mats)
}

prefilter_for_knockout <- function(counts, must_keep, min_cells = 50, top_n_var = 1000) {
  detected <- Matrix::rowSums(counts > 0)
  keep <- detected >= min_cells
  norm <- normalize_log1p(counts[keep, , drop = FALSE])
  vars <- apply(norm, 1, var)
  vars <- sort(vars, decreasing = TRUE)
  top_genes <- names(vars)[seq_len(min(top_n_var, length(vars)))]
  selected <- unique(c(top_genes, intersect(must_keep, rownames(counts))))
  counts[selected, , drop = FALSE]
}

gene_set_summary <- function(diff_tbl, gene_sets) {
  diff_tbl$rank_padj <- rank(diff_tbl$p.adj, ties.method = "average")
  diff_tbl$rank_z <- rank(-abs(diff_tbl$Z), ties.method = "average")
  rows <- list()
  top_sizes <- c(25, 50, 100, 200)

  for (set_name in names(gene_sets)) {
    genes <- intersect(gene_sets[[set_name]], diff_tbl$gene)
    if (length(genes) == 0) {
      next
    }
    sub <- diff_tbl[diff_tbl$gene %in% genes, , drop = FALSE]
    for (top_n in top_sizes) {
      top_genes <- head(diff_tbl$gene[order(diff_tbl$p.adj, -abs(diff_tbl$Z))], top_n)
      in_top <- sum(genes %in% top_genes)
      out_top <- top_n - in_top
      in_rest <- length(genes) - in_top
      out_rest <- nrow(diff_tbl) - top_n - in_rest
      p_val <- tryCatch(fisher.test(matrix(c(in_top, out_top, in_rest, out_rest), nrow = 2))$p.value, error = function(e) NA_real_)
      rows[[length(rows) + 1]] <- data.frame(
        gene_set = set_name,
        top_n = top_n,
        set_size = length(genes),
        in_top = in_top,
        fisher_p = p_val,
        mean_abs_z = NA_real_,
        median_padj = NA_real_,
        mean_fc = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    rows[[length(rows) + 1]] <- data.frame(
      gene_set = paste0(set_name, "_summary"),
      top_n = NA_integer_,
      set_size = length(genes),
      in_top = NA_integer_,
      fisher_p = NA_real_,
      mean_abs_z = mean(abs(sub$Z), na.rm = TRUE),
      median_padj = median(sub$p.adj, na.rm = TRUE),
      mean_fc = mean(sub$FC, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

run_knockout <- function(counts, label, gene_sets, n_cores = 4, nc_nNet = 2, nc_nCells = 150) {
  message("Running scTenifoldKnk for ", label, " (", nrow(counts), " genes x ", ncol(counts), " cells)")
  res <- scTenifoldKnk(
    countMatrix = as.matrix(counts),
    gKO = "NPL",
    qc = TRUE,
    qc_minLSize = 0,
    qc_minCells = 25,
    nc_nNet = nc_nNet,
    nc_nCells = nc_nCells,
    nc_nComp = 3,
    td_K = 3,
    ma_nDim = 2,
    nCores = n_cores
  )
  diff_tbl <- res$diffRegulation
  diff_tbl$context <- label
  diff_tbl <- diff_tbl[order(diff_tbl$p.adj, -abs(diff_tbl$Z)), ]
  enrichment <- gene_set_summary(diff_tbl, gene_sets)
  list(result = res, diff = diff_tbl, enrich = enrichment)
}

cell_typing <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_cell_level_typing.tsv.gz"))
modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))
state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))

compact_module <- modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"]
raw_module <- modules$gene_symbol[modules$module_name == "NPL_NEIGHBORHOOD_RAW"]
foam_set <- state_sets$gene_symbol[state_sets$state_set == "MACROPHAGE_FOAM_TREM2"]
inflam_set <- state_sets$gene_symbol[state_sets$state_set == "MACROPHAGE_INFLAMMATORY"]
c1q_set <- state_sets$gene_symbol[state_sets$state_set == "MACROPHAGE_C1Q"]

must_keep <- unique(c(compact_module, raw_module, foam_set, inflam_set, c1q_set, "NPL"))
gene_sets <- list(
  NPL_FOAM_MACROPHAGE_COMPACT = compact_module,
  NPL_NEIGHBORHOOD_RAW = raw_module,
  MACROPHAGE_FOAM_TREM2 = foam_set,
  MACROPHAGE_INFLAMMATORY = inflam_set,
  MACROPHAGE_C1Q = c1q_set
)

core_mac_counts <- build_macrophage_counts(cell_typing, keep_location = "core")

core_mac_counts <- prefilter_for_knockout(core_mac_counts, must_keep = must_keep)

core_res <- run_knockout(core_mac_counts, "core_macrophage", gene_sets, nc_nNet = 2, nc_nCells = 150)

write_tsv(core_res$diff, project_path("res", "tables", "mechanism", "npl_sctenifoldknk_core_macrophage_diff.tsv"))
write_tsv(core_res$enrich, project_path("res", "tables", "mechanism", "npl_sctenifoldknk_core_macrophage_enrichment.tsv"))
write_tsv(head(core_res$diff, 100), project_path("res", "qc", "mechanism", "npl_sctenifoldknk_core_top100.tsv"))

cat("scTenifoldKnk results written to res/tables/mechanism and res/qc/mechanism\n")
cat("\nCore macrophage enrichment summary:\n")
print(core_res$enrich)
cat("\nTop core perturbation hits:\n")
print(head(core_res$diff, 30))
