source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
  library(FNN)
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

infer_patient_id <- function(sample_id) {
  ifelse(
    sample_id %in% c("GSM4837523", "GSM4837524"),
    "P1",
    ifelse(sample_id %in% c("GSM4837525", "GSM4837526"), "P2", "P3")
  )
}

build_macrophage_object <- function() {
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
    colnames(counts) <- paste(sample_name, keep, sep = "|")
    mats[[sample_name]] <- counts

    meta_rows[[sample_name]] <- data.frame(
      cell_id = colnames(counts),
      sample_name = sample_name,
      sample_id = sample_id,
      patient_id = infer_patient_id(sample_id),
      location = pheno$group_label[match(sample_id, pheno$sample_id)],
      stringsAsFactors = FALSE
    )
  }

  list(
    counts = do.call(cbind, mats),
    meta = do.call(rbind, meta_rows)
  )
}

compute_hvg <- function(norm_mat, n_top = 1500) {
  detected <- Matrix::rowSums(norm_mat > 0) >= 50
  keep <- which(detected)
  if (length(keep) == 0) {
    stop("No genes passed detection filter")
  }
  sub <- as.matrix(norm_mat[keep, , drop = FALSE])
  gene_var <- apply(sub, 1, var)
  gene_mean <- rowMeans(sub)
  disp <- gene_var / pmax(gene_mean, 1e-3)
  ord <- order(disp, decreasing = TRUE)
  keep[ord[seq_len(min(n_top, length(ord)))]]
}

build_knn_graph <- function(pc_mat, k = 20) {
  knn <- FNN::get.knn(pc_mat, k = k)
  edge_from <- rep(seq_len(nrow(pc_mat)), each = k)
  edge_to <- as.vector(knn$nn.index)
  edge_weight <- as.vector(knn$nn.dist) + 1e-6
  edge_df <- data.frame(from = edge_from, to = edge_to, weight = edge_weight)
  edge_df <- edge_df[edge_df$from != edge_df$to, , drop = FALSE]
  graph <- igraph::graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = seq_len(nrow(pc_mat))))
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                            edge.attr.comb = list(weight = "min"))
  graph
}

compute_rooted_pseudotime <- function(graph, root_ids) {
  ext_graph <- igraph::add_vertices(graph, 1, name = "seed")
  seed_edges <- as.vector(rbind(rep("seed", length(root_ids)), as.character(root_ids)))
  ext_graph <- igraph::add_edges(ext_graph, seed_edges, attr = list(weight = rep(1e-9, length(root_ids))))
  d <- igraph::distances(ext_graph, v = "seed", to = as.character(seq_len(igraph::gorder(graph))), weights = igraph::E(ext_graph)$weight)
  pt <- as.numeric(d[1, ])
  pt <- pt - min(pt, na.rm = TRUE)
  pt / max(pt, na.rm = TRUE)
}

compute_projection_axis <- function(pc_mat, root_ids, terminal_ids) {
  root_center <- colMeans(pc_mat[root_ids, , drop = FALSE])
  terminal_center <- colMeans(pc_mat[terminal_ids, , drop = FALSE])
  axis_vec <- terminal_center - root_center
  axis_norm <- sqrt(sum(axis_vec ^ 2))
  if (!is.finite(axis_norm) || axis_norm == 0) {
    stop("Projection axis could not be defined")
  }
  axis_unit <- axis_vec / axis_norm
  proj <- as.numeric((pc_mat - matrix(root_center, nrow(pc_mat), ncol(pc_mat), byrow = TRUE)) %*% axis_unit)
  proj <- proj - min(proj, na.rm = TRUE)
  proj / max(proj, na.rm = TRUE)
}

orient_axis <- function(axis_vals, meta) {
  cor_val <- suppressWarnings(cor(axis_vals, meta$foam_minus_c1q, method = "spearman"))
  if (!is.na(cor_val) && cor_val < 0) {
    axis_vals <- 1 - axis_vals
  }
  core_mean <- mean(axis_vals[meta$location == "core"], na.rm = TRUE)
  adjacent_mean <- mean(axis_vals[meta$location == "adjacent"], na.rm = TRUE)
  if (is.finite(core_mean) && is.finite(adjacent_mean) && core_mean < adjacent_mean) {
    axis_vals <- 1 - axis_vals
  }
  axis_vals
}

obj <- build_macrophage_object()
norm <- normalize_log1p(obj$counts)
meta <- obj$meta

state_sets <- read_tsv_auto(project_path("res", "tables", "mechanism", "vascular_state_gene_sets.tsv"))
modules <- read_tsv_auto(project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))

mac_states <- c("MACROPHAGE_C1Q", "MACROPHAGE_FOAM_TREM2", "MACROPHAGE_INFLAMMATORY", "MACROPHAGE_IFN")
state_scores <- sapply(mac_states, function(set_name) {
  score_gene_set(norm, state_sets$gene_symbol[state_sets$state_set == set_name])
})
if (is.vector(state_scores)) {
  state_scores <- matrix(state_scores, ncol = 1, dimnames = list(colnames(norm), mac_states))
}

meta$dominant_state <- colnames(state_scores)[max.col(state_scores, ties.method = "first")]
meta$c1q_score <- state_scores[, "MACROPHAGE_C1Q"]
meta$foam_score <- state_scores[, "MACROPHAGE_FOAM_TREM2"]
meta$inflammatory_score <- state_scores[, "MACROPHAGE_INFLAMMATORY"]
meta$ifn_score <- state_scores[, "MACROPHAGE_IFN"]
meta$foam_minus_c1q <- meta$foam_score - meta$c1q_score
meta$npl_module_score <- score_gene_set(
  norm,
  modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"]
)

hvg_idx <- compute_hvg(norm, n_top = 1500)
hvg_mat <- as.matrix(norm[hvg_idx, , drop = FALSE])
hvg_scaled <- t(scale(t(hvg_mat)))
hvg_scaled[!is.finite(hvg_scaled)] <- 0

pca <- prcomp(t(hvg_scaled), center = FALSE, scale. = FALSE, rank. = 15)
pc_mat <- pca$x[, 1:15, drop = FALSE]

adjacent_root_cut <- quantile(meta$c1q_score[meta$location == "adjacent"], probs = 0.75, na.rm = TRUE)
foam_root_cut <- quantile(meta$foam_score[meta$location == "adjacent"], probs = 0.50, na.rm = TRUE)
root_ids <- which(meta$location == "adjacent" & meta$c1q_score >= adjacent_root_cut & meta$foam_score <= foam_root_cut)
if (length(root_ids) < 50) {
  root_ids <- which(meta$location == "adjacent" & meta$dominant_state == "MACROPHAGE_C1Q")
}

graph <- build_knn_graph(pc_mat, k = 20)
graph_pseudotime <- compute_rooted_pseudotime(graph, root_ids)

terminal_ids <- which(
  meta$location == "core" &
    meta$foam_score >= quantile(meta$foam_score[meta$location == "core"], probs = 0.75, na.rm = TRUE) &
    meta$c1q_score <= quantile(meta$c1q_score[meta$location == "core"], probs = 0.50, na.rm = TRUE)
)
if (length(terminal_ids) < 50) {
  terminal_ids <- which(meta$location == "core" & meta$dominant_state == "MACROPHAGE_FOAM_TREM2")
}
if (length(terminal_ids) < 50) {
  terminal_ids <- which(meta$dominant_state == "MACROPHAGE_FOAM_TREM2")
}

projection_pseudotime <- compute_projection_axis(pc_mat, root_ids, terminal_ids)
graph_pseudotime <- orient_axis(graph_pseudotime, meta)
projection_pseudotime <- orient_axis(projection_pseudotime, meta)

graph_cor <- suppressWarnings(cor(graph_pseudotime, meta$foam_minus_c1q, method = "spearman"))
proj_cor <- suppressWarnings(cor(projection_pseudotime, meta$foam_minus_c1q, method = "spearman"))
graph_sep <- mean(graph_pseudotime[meta$location == "core"], na.rm = TRUE) -
  mean(graph_pseudotime[meta$location == "adjacent"], na.rm = TRUE)
proj_sep <- mean(projection_pseudotime[meta$location == "core"], na.rm = TRUE) -
  mean(projection_pseudotime[meta$location == "adjacent"], na.rm = TRUE)

if (is.na(graph_cor)) graph_cor <- -Inf
if (is.na(proj_cor)) proj_cor <- -Inf

graph_score <- graph_cor + graph_sep
proj_score <- proj_cor + proj_sep

if (proj_score >= graph_score) {
  meta$pseudotime <- projection_pseudotime
  pseudotime_method <- "projection_axis"
} else {
  meta$pseudotime <- graph_pseudotime
  pseudotime_method <- "rooted_knn_graph"
}

meta$pc_1 <- pc_mat[, 1]
meta$pc_2 <- pc_mat[, 2]

core_integrated <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_core_vs_adjacent_integrated.tsv"))
foam_markers <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_state_markers.tsv"))
foam_markers <- foam_markers[foam_markers$dominant_state == "MACROPHAGE_FOAM_TREM2", ]

candidate_genes <- unique(c(
  "NPL", "FABP5", "GPNMB", "APOC1", "SPP1", "CD36", "APOE", "PLA2G7", "LGALS3", "HMOX1", "CYP27A1",
  head(core_integrated$gene_symbol, 120),
  head(foam_markers$gene_symbol[order(-foam_markers$delta_state_vs_other)], 120),
  modules$gene_symbol[modules$module_name == "NPL_FOAM_MACROPHAGE_COMPACT"]
))
candidate_genes <- intersect(unique(candidate_genes), rownames(norm))

gene_cor_rows <- lapply(candidate_genes, function(gene_symbol) {
  vals <- as.numeric(norm[gene_symbol, ])
  data.frame(
    gene_symbol = gene_symbol,
    rho_pseudotime = suppressWarnings(cor(vals, meta$pseudotime, method = "spearman")),
    rho_foam_minus_c1q = suppressWarnings(cor(vals, meta$foam_minus_c1q, method = "spearman")),
    stringsAsFactors = FALSE
  )
})
gene_cor <- do.call(rbind, gene_cor_rows)
gene_cor <- gene_cor[order(-gene_cor$rho_pseudotime), ]
gene_cor$rank_pseudotime <- seq_len(nrow(gene_cor))
candidate_cor <- gene_cor

meta$pt_bin <- cut(meta$pseudotime, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
trend_rows <- list()
for (gene_symbol in candidate_genes) {
  tmp <- data.frame(
    gene_symbol = gene_symbol,
    pt_bin = meta$pt_bin,
    expr = as.numeric(norm[gene_symbol, ]),
    stringsAsFactors = FALSE
  )
  trend_rows[[gene_symbol]] <- aggregate(expr ~ gene_symbol + pt_bin, data = tmp, FUN = mean)
}
trend_table <- do.call(rbind, trend_rows)

sample_pt <- aggregate(
  cbind(pseudotime, foam_minus_c1q, npl_module_score) ~ sample_id + patient_id + location + dominant_state,
  data = meta,
  FUN = mean
)

group_pt <- aggregate(
  cbind(pseudotime, foam_minus_c1q, npl_module_score) ~ location + dominant_state,
  data = meta,
  FUN = mean
)

paired_summary <- merge(
  aggregate(cbind(pseudotime, foam_minus_c1q, npl_module_score) ~ patient_id + location, data = meta, FUN = mean),
  aggregate(cbind(pseudotime, foam_minus_c1q, npl_module_score) ~ patient_id + location, data = meta, FUN = median),
  by = c("patient_id", "location"),
  suffixes = c("_mean", "_median")
)

global_summary <- data.frame(
  n_cells = nrow(meta),
  n_root_cells = length(root_ids),
  n_terminal_cells = length(terminal_ids),
  pseudotime_method = pseudotime_method,
  pseudotime_location_wilcox_p = wilcox.test(meta$pseudotime ~ meta$location, exact = FALSE)$p.value,
  pseudotime_core_mean = mean(meta$pseudotime[meta$location == "core"]),
  pseudotime_adjacent_mean = mean(meta$pseudotime[meta$location == "adjacent"]),
  pseudotime_foam_axis_spearman = suppressWarnings(cor(meta$pseudotime, meta$foam_minus_c1q, method = "spearman")),
  npl_module_pseudotime_spearman = suppressWarnings(cor(meta$npl_module_score, meta$pseudotime, method = "spearman")),
  graph_axis_spearman = graph_cor,
  projection_axis_spearman = proj_cor,
  graph_core_minus_adjacent = graph_sep,
  projection_core_minus_adjacent = proj_sep,
  stringsAsFactors = FALSE
)

write_tsv(meta, project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_cell_table.tsv"))
write_tsv(gene_cor, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudotime_gene_correlations.tsv"))
write_tsv(candidate_cor, project_path("res", "qc", "mechanism", "gse159677_macrophage_pseudotime_candidate_correlations.tsv"))
write_tsv(trend_table, project_path("res", "tables", "mechanism", "gse159677_macrophage_pseudotime_trends.tsv"))
write_tsv(sample_pt, project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_sample_summary.tsv"))
write_tsv(group_pt, project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_group_summary.tsv"))
write_tsv(paired_summary, project_path("res", "tables", "mechanism", "gse159677_macrophage_continuum_paired_summary.tsv"))
write_tsv(global_summary, project_path("res", "qc", "mechanism", "gse159677_macrophage_continuum_global_summary.tsv"))

p_state <- ggplot(meta, aes(pc_1, pc_2, color = dominant_state)) +
  geom_point(size = 0.25, alpha = 0.6) +
  theme_bw(base_size = 10) +
  labs(title = "GSE159677 macrophage continuum", subtitle = "Dominant macrophage state", x = "PC1", y = "PC2", color = "State")

p_pt <- ggplot(meta, aes(pc_1, pc_2, color = pseudotime)) +
  geom_point(size = 0.25, alpha = 0.7) +
  scale_color_viridis_c(option = "C") +
  theme_bw(base_size = 10) +
  labs(title = "GSE159677 macrophage continuum", subtitle = "Rooted graph pseudotime", x = "PC1", y = "PC2", color = "Pseudotime")

p_box <- ggplot(meta, aes(location, pseudotime, fill = location)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.size = 0.2) +
  facet_wrap(~ dominant_state, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = "Macrophage pseudotime by location and dominant state", x = NULL, y = "Pseudotime")

ggsave(project_path("figure", "export", "mechanism", "gse159677_macrophage_continuum_state_umap.pdf"), p_state, width = 7.5, height = 5.8)
ggsave(project_path("figure", "export", "mechanism", "gse159677_macrophage_continuum_pseudotime_umap.pdf"), p_pt, width = 7.5, height = 5.8)
ggsave(project_path("figure", "export", "mechanism", "gse159677_macrophage_continuum_boxplot.pdf"), p_box, width = 8.5, height = 5.2)

cat("Macrophage continuum outputs written to res/tables/mechanism and figure/export/mechanism\n")
print(global_summary)
print(candidate_cor)
