#!/usr/bin/env Rscript

# Reproducible issue #27 acceptance benchmark for a prepared UKB-PPP NTRK3
# protein. Raw input, stores, logs, and the /usr/bin/time record are supplied
# through external paths and must remain outside the synced repository.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L || !args[[1L]] %in% c("full", "core", "core_plus") ||
    !args[[2L]] %in% c("compact", "none")) {
  stop("usage: benchmark_issue27_ntrk3.R full|core|core_plus compact|none INPUT OUTPUT SUMMARY.json",
       call. = FALSE)
}

selection <- args[[1L]]
qc_mode <- args[[2L]]
input <- args[[3L]]
output <- args[[4L]]
summary_path <- args[[5L]]
threads <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
commit <- Sys.getenv("COMPRESSOR_ISSUE27_COMMIT", unset = "unknown")
actual_commit <- Sys.getenv("COMPRESSOR_ISSUE27_ACTUAL_COMMIT", unset = commit)
native_commit <- Sys.getenv("COMPRESSOR_ISSUE27_NATIVE_LIBRARY_COMMIT", unset = "unknown")

if (!file.exists(input)) stop("input does not exist: ", input, call. = FALSE)
if (!grepl("^[0-9a-fA-F]{40}$", commit)) {
  stop("COMPRESSOR_ISSUE27_COMMIT must be a 40-character commit", call. = FALSE)
}
if (!grepl("^[0-9a-fA-F]{40}$", actual_commit) || !identical(actual_commit, commit)) {
  stop("COMPRESSOR_ISSUE27_ACTUAL_COMMIT must match COMPRESSOR_ISSUE27_COMMIT", call. = FALSE)
}

suppressPackageStartupMessages(library(CompreSSoR))

`%||%` <- function(x, y) if (is.null(x)) y else x
started <- unname(proc.time()[["elapsed"]])
store <- compress_sumstats(
  input, output,
  input_build = "GRCh38", store_build = "GRCh38",
  selection = selection, qc = qc_mode, row_policy = "report",
  threads = threads, overwrite = TRUE
)
elapsed <- unname(proc.time()[["elapsed"]]) - started
validation <- validate_compressor(store, full = FALSE)
manifest <- store$manifest
prep <- manifest$preparation$preparation %||% list()
structural <- prep$structural_qc %||% list()
selection_meta <- manifest$selection %||% list()
source_columns <- manifest$source_columns %||% character()
source_columns_read <- manifest$source_columns_read %||% character()
payload <- list.files(output, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)

source_input_rows <- prep$input_rows %||%
  structural$input_rows %||% selection_meta$source_input_rows %||% NA_integer_
accepted_rows <- structural$accepted_rows %||%
  structural$kept_rows %||% selection_meta$input_rows %||% NA_integer_
eaf <- manifest$eaf_observability %||% list()
writer <- manifest$writer %||% list()

result <- list(
  status = "PASS",
  benchmark = "issue27_ntrk3_prepared",
  commit = commit,
  actual_commit = actual_commit,
  native_library_commit = native_commit,
  package_version = as.character(packageVersion("CompreSSoR")),
  host = Sys.info()[["nodename"]],
  job_id = Sys.getenv("SLURM_JOB_ID", ""),
  source = list(
    path = normalizePath(input, mustWork = TRUE),
    basename = basename(input),
    bytes = as.numeric(file.info(input)$size),
    sha256 = manifest$source$sha256 %||% NA_character_,
    rows = as.integer(source_input_rows),
    columns_before = manifest$source$columns_before %||% length(source_columns),
    columns_read = manifest$source$columns_read %||% length(source_columns_read),
    columns = source_columns,
    columns_read_names = source_columns_read,
    projected = manifest$source$projected %||% NA
  ),
  build = list(input = "GRCh38", stored = "GRCh38"),
  parameters = list(
    selection = selection,
    qc = qc_mode,
    row_policy = manifest$row_policy %||% "report",
    threads_requested = as.integer(threads),
    threads_effective = manifest$threads$effective %||% writer$effective_workers %||% NA_integer_,
    writer_effective_workers = writer$effective_workers %||% NA_integer_
  ),
  result = list(
    elapsed_seconds = as.numeric(elapsed),
    stored_rows = as.integer(manifest$n_rows %||% manifest$rows),
    source_rows = as.integer(source_input_rows),
    accepted_rows_after_qc_or_identity = as.integer(accepted_rows),
    rejected_rows = as.integer(structural$rejected_rows %||%
                               structural$dropped_rows %||% NA_integer_),
    output_bytes = as.numeric(sum(file.info(payload)$size, na.rm = TRUE)),
    exception_rows = as.integer(manifest$semantic_codec$exception_rows %||% NA_integer_),
    eaf_coverage = eaf,
    selection = selection_meta,
    timings = manifest$timings %||% NULL,
    validation = validation
  ),
  provenance = list(
    format = manifest$format_version %||% NA_character_,
    profile = manifest$profile %||% NA_character_,
    codec = manifest$codec$name %||% NA_character_,
    external_reference_required = manifest$identity$external_reference_required %||% FALSE,
    input_build = manifest$input_build %||% "GRCh38",
    stored_build = manifest$genome_build %||% "GRCh38"
  )
)

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(result, summary_path, auto_unbox = TRUE, pretty = TRUE,
                     null = "null", digits = 17)
cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE,
                     null = "null", digits = 17), "\n")
