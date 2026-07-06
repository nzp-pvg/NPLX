source("script/R/00_project_config.R")

ensure_project_dirs()

npl_cor <- read_tsv_auto(project_path("res", "qc", "mechanism", "gse159677_npl_core_macrophage_top100.tsv"))
refined <- read_tsv_auto(project_path("res", "tables", "mechanism", "gse159677_macrophage_refined_candidates.tsv"))

raw_candidates <- merge(
  npl_cor,
  refined[, c(
    "gene_symbol",
    "foam_delta_vs_other",
    "n_patients_core_gt_adj",
    "bulk_direction_consistency",
    "external_top_cell_type",
    "external_macrophage_rank",
    "external_macrophage_specificity_ratio"
  )],
  by = "gene_symbol",
  all.x = TRUE
)

keep_raw <- !is.na(raw_candidates$gene_symbol) &
  !is.na(raw_candidates$foam_delta_vs_other) &
  !is.na(raw_candidates$n_patients_core_gt_adj) &
  !is.na(raw_candidates$bulk_direction_consistency) &
  raw_candidates$rank_with_npl <= 40 &
  raw_candidates$foam_delta_vs_other >= 0.3 &
  raw_candidates$n_patients_core_gt_adj == 3 &
  raw_candidates$bulk_direction_consistency == 3
raw_candidates <- raw_candidates[keep_raw, ]
raw_candidates <- raw_candidates[order(raw_candidates$rank_with_npl), ]

raw_module <- data.frame(
  module_name = "NPL_NEIGHBORHOOD_RAW",
  gene_symbol = raw_candidates$gene_symbol,
  selection_reason = "top_npl_neighbor_with_state_and_bulk_support",
  source_rank = raw_candidates$rank_with_npl,
  stringsAsFactors = FALSE
)

# Curated compact module removes broad lysosomal / housekeeping-heavy neighbors and
# keeps the disease-facing foam-remodeling genes most suitable for a mechanistic story.
compact_genes <- c(
  "NPL",
  "FABP5",
  "GPNMB",
  "APOC1",
  "PLA2G7",
  "SPP1",
  "CD36",
  "CYP27A1",
  "APOE"
)

compact_module <- data.frame(
  module_name = "NPL_FOAM_MACROPHAGE_COMPACT",
  gene_symbol = compact_genes,
  selection_reason = c(
    "underexplored_lead",
    rep("canonical_foam_neighbor", length(compact_genes) - 1)
  ),
  source_rank = match(compact_genes, npl_cor$gene_symbol),
  stringsAsFactors = FALSE
)

module_table <- rbind(raw_module, compact_module)

write_tsv(raw_candidates, project_path("res", "tables", "mechanism", "npl_raw_neighborhood_candidates.tsv"))
write_tsv(module_table, project_path("res", "tables", "mechanism", "npl_macrophage_modules.tsv"))

cat("NPL macrophage modules written to res/tables/mechanism\n")
print(module_table)
