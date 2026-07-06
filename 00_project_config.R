resolve_this_file <- function() {
  frame_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(frame_file) && nzchar(frame_file)) {
    return(normalizePath(frame_file, winslash = "/", mustWork = TRUE))
  }

  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(cmd_file) > 0) {
    return(normalizePath(sub("^--file=", "", cmd_file[[1]]), winslash = "/", mustWork = TRUE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

this_file <- resolve_this_file()
root_dir <- normalizePath(file.path(dirname(this_file), "..", ".."), winslash = "/", mustWork = TRUE)

project_path <- function(...) {
  file.path(root_dir, ...)
}

ensure_dir <- function(...) {
  dir.create(project_path(...), recursive = TRUE, showWarnings = FALSE)
}

ensure_project_dirs <- function() {
  dirs <- list(
    c("data", "processed", "annotation"),
    c("data", "processed", "bulk"),
    c("data", "processed", "bulk_gene"),
    c("data", "processed", "methylation"),
    c("res", "qc", "bulk"),
    c("res", "qc", "mechanism"),
    c("res", "tables", "bulk"),
    c("res", "tables", "mechanism"),
    c("figure", "export", "bulk"),
    c("figure", "export", "mechanism")
  )

  for (parts in dirs) {
    do.call(ensure_dir, as.list(parts))
  }
}

read_tsv_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    return(read.delim(gzfile(path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  }
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

read_geo_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  begin_idx <- which(lines == "!series_matrix_table_begin")
  end_idx <- which(lines == "!series_matrix_table_end")

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    stop("Could not find a valid series_matrix table block in: ", path)
  }

  block <- lines[(begin_idx + 1):(end_idx - 1)]
  handle <- textConnection(block)
  on.exit(close(handle), add = TRUE)

  mat <- read.delim(
    handle,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  colnames(mat)[1] <- "feature_id"
  mat
}

read_geo_platform_table <- function(path) {
  lines <- readLines(path, warn = FALSE)
  begin_idx <- which(lines == "!platform_table_begin")
  end_idx <- which(lines == "!platform_table_end")

  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    stop("Could not find a valid platform table block in: ", path)
  }

  block <- lines[(begin_idx + 1):(end_idx - 1)]
  handle <- textConnection(block)
  on.exit(close(handle), add = TRUE)

  read.delim(
    handle,
    sep = "\t",
    header = TRUE,
    quote = "\"",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
}

load_soft_samples <- function(dataset_id) {
  read_tsv_auto(project_path("data", "metadata", "soft_samples", paste0(dataset_id, "_samples.tsv")))
}

load_pheno <- function(dataset_id) {
  read_tsv_auto(project_path("data", "processed", "pheno", paste0(dataset_id, "_pheno.tsv")))
}

clean_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "---", "NA", "NULL")] <- NA_character_
  x
}
