source("script/R/00_project_config.R")
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

ensure_project_dirs()

first_non_empty <- function(values) {
  values <- trimws(values)
  values <- values[!(is.na(values) | values == "" | values == "---" | values == "NULL")]
  if (length(values) == 0) {
    return(NA_character_)
  }
  values[[1]]
}

collapse_unique_values <- function(values) {
  values <- trimws(values)
  values <- unique(values[!(is.na(values) | values == "" | values == "---" | values == "NULL")])
  if (length(values) == 0) {
    return(NA_character_)
  }
  paste(values, collapse = " /// ")
}

parse_assignment_field <- function(field_value) {
  if (is.na(field_value) || field_value == "" || field_value == "---") {
    return(c(
      gene_symbol = NA_character_,
      gene_symbol_all = NA_character_,
      entrez_id = NA_character_,
      entrez_id_all = NA_character_,
      gene_title = NA_character_,
      gene_title_all = NA_character_
    ))
  }

  entries <- strsplit(field_value, " /// ", fixed = TRUE)[[1]]
  parts <- lapply(entries, function(entry) trimws(strsplit(entry, " // ", fixed = TRUE)[[1]]))

  symbols <- vapply(parts, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1))
  titles <- vapply(parts, function(x) if (length(x) >= 3) x[[3]] else NA_character_, character(1))
  entrez <- vapply(parts, function(x) if (length(x) >= 5) x[[5]] else NA_character_, character(1))

  c(
    gene_symbol = first_non_empty(symbols),
    gene_symbol_all = collapse_unique_values(symbols),
    entrez_id = first_non_empty(entrez),
    entrez_id_all = collapse_unique_values(entrez),
    gene_title = first_non_empty(titles),
    gene_title_all = collapse_unique_values(titles)
  )
}

pick_first_token <- function(x, sep = " /// ") {
  ifelse(
    is.na(x) | x == "",
    NA_character_,
    trimws(vapply(strsplit(x, sep, fixed = TRUE), `[`, character(1), 1))
  )
}

pick_first_semicolon <- function(x) {
  ifelse(
    is.na(x) | x == "",
    NA_character_,
    trimws(vapply(strsplit(x, ";", fixed = TRUE), `[`, character(1), 1))
  )
}

gpl6244 <- read_geo_platform_table(project_path("data", "raw", "platform_tables", "GPL6244_platform_table.txt"))
gpl6244_parsed <- t(vapply(gpl6244$gene_assignment, parse_assignment_field, character(6)))
gpl6244_map <- data.frame(
  platform_id = "GPL6244",
  feature_id = gpl6244$ID,
  gene_symbol = clean_gene_symbol(gpl6244_parsed[, "gene_symbol"]),
  gene_symbol_all = clean_gene_symbol(gpl6244_parsed[, "gene_symbol_all"]),
  entrez_id = clean_gene_symbol(gpl6244_parsed[, "entrez_id"]),
  entrez_id_all = clean_gene_symbol(gpl6244_parsed[, "entrez_id_all"]),
  gene_title = clean_gene_symbol(gpl6244_parsed[, "gene_title"]),
  gene_title_all = clean_gene_symbol(gpl6244_parsed[, "gene_title_all"]),
  raw_annotation = gpl6244$gene_assignment,
  stringsAsFactors = FALSE
)

gpl17077 <- read_geo_platform_table(project_path("data", "raw", "platform_tables", "GPL17077_platform_table.txt"))
gpl17077_map <- data.frame(
  platform_id = "GPL17077",
  feature_id = gpl17077$ID,
  gene_symbol = clean_gene_symbol(gpl17077$GENE_SYMBOL),
  gene_symbol_all = clean_gene_symbol(gpl17077$GENE_SYMBOL),
  entrez_id = clean_gene_symbol(as.character(gpl17077$LOCUSLINK_ID)),
  entrez_id_all = clean_gene_symbol(as.character(gpl17077$LOCUSLINK_ID)),
  gene_title = clean_gene_symbol(gpl17077$GENE_NAME),
  gene_title_all = clean_gene_symbol(gpl17077$GENE_NAME),
  raw_annotation = gpl17077$DESCRIPTION,
  stringsAsFactors = FALSE
)
gpl17077_map <- gpl17077_map[!grepl("CONTROL|Corner", gpl17077_map$feature_id, ignore.case = TRUE), , drop = FALSE]

gpl570 <- read_geo_platform_table(project_path("data", "raw", "platform_tables", "GPL570_platform_table.txt"))
gpl570_map <- data.frame(
  platform_id = "GPL570",
  feature_id = gpl570$ID,
  gene_symbol = clean_gene_symbol(pick_first_token(gpl570$`Gene Symbol`)),
  gene_symbol_all = clean_gene_symbol(gpl570$`Gene Symbol`),
  entrez_id = clean_gene_symbol(pick_first_token(as.character(gpl570$ENTREZ_GENE_ID))),
  entrez_id_all = clean_gene_symbol(as.character(gpl570$ENTREZ_GENE_ID)),
  gene_title = clean_gene_symbol(pick_first_token(gpl570$`Gene Title`)),
  gene_title_all = clean_gene_symbol(gpl570$`Gene Title`),
  raw_annotation = gpl570$`Target Description`,
  stringsAsFactors = FALSE
)

gpl13534 <- read_geo_platform_table(project_path("data", "raw", "platform_tables", "GPL13534_platform_table.txt"))
gpl13534_symbols <- pick_first_semicolon(gpl13534$UCSC_RefGene_Name)
gpl13534_titles <- mapIds(
  org.Hs.eg.db,
  keys = gpl13534_symbols,
  column = "GENENAME",
  keytype = "SYMBOL",
  multiVals = "first"
)
gpl13534_entrez <- mapIds(
  org.Hs.eg.db,
  keys = gpl13534_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
gpl13534_map <- data.frame(
  platform_id = "GPL13534",
  feature_id = gpl13534$ID,
  gene_symbol = clean_gene_symbol(gpl13534_symbols),
  gene_symbol_all = clean_gene_symbol(gpl13534$UCSC_RefGene_Name),
  entrez_id = clean_gene_symbol(as.character(gpl13534_entrez)),
  entrez_id_all = clean_gene_symbol(as.character(gpl13534_entrez)),
  gene_title = clean_gene_symbol(as.character(gpl13534_titles)),
  gene_title_all = clean_gene_symbol(as.character(gpl13534_titles)),
  raw_annotation = gpl13534$UCSC_RefGene_Accession,
  stringsAsFactors = FALSE
)

mapping_list <- list(
  GPL6244 = gpl6244_map,
  GPL17077 = gpl17077_map,
  GPL570 = gpl570_map,
  GPL13534 = gpl13534_map
)

summary_rows <- list()

for (mapping_name in names(mapping_list)) {
  mapping <- mapping_list[[mapping_name]]
  out_file <- project_path("data", "processed", "annotation", paste0(mapping_name, "_mapping.tsv.gz"))
  write_tsv_gz(mapping, out_file)
  summary_rows[[mapping_name]] <- data.frame(
    platform_id = mapping_name,
    n_rows = nrow(mapping),
    n_nonempty_gene_symbol = sum(!(is.na(mapping$gene_symbol) | mapping$gene_symbol == "")),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
write_tsv(summary_table, project_path("res", "qc", "bulk", "platform_annotation_summary.tsv"))

cat("Platform annotation tables written to data/processed/annotation\n")
print(summary_table)
