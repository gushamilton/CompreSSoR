#!/usr/bin/env Rscript

# Proportionate core-plus acceptance benchmark on the prepared UKB-PPP NTRK3
# input. Raw input, stores, logs, and summaries are supplied through external
# paths and must remain outside the synced repository.

suppressPackageStartupMessages(library(CompreSSoR))
`%||%` <- function(x, y) if (is.null(x)) y else x

source_root <- Sys.getenv("COMPRESSOR_CORE_PLUS_SOURCE_ROOT", unset = "")
if (nzchar(source_root)) {
  source_env <- new.env(parent = asNamespace("CompreSSoR"))
  source_files <- list.files(file.path(source_root, "R"), pattern = "[.]R$",
                             full.names = TRUE)
  for (source_file in sort(source_files)) sys.source(source_file, source_env)
  compress_sumstats <- source_env$compress_sumstats
} else {
  compress_sumstats <- CompreSSoR::compress_sumstats
}

input_path <- Sys.getenv("COMPRESSOR_CORE_PLUS_INPUT")
output_path <- Sys.getenv("COMPRESSOR_CORE_PLUS_STORE")
result_path <- Sys.getenv("COMPRESSOR_CORE_PLUS_RESULT")
commit <- Sys.getenv("COMPRESSOR_CORE_PLUS_COMMIT", unset = "unknown")
threads <- as.integer(Sys.getenv("COMPRESSOR_CORE_PLUS_THREADS", unset = "4"))
pvalue_threshold <- as.numeric(Sys.getenv(
  "COMPRESSOR_CORE_PLUS_PVALUE_THRESHOLD", unset = "1e-5"
))
region_padding <- as.integer(Sys.getenv(
  "COMPRESSOR_CORE_PLUS_REGION_PADDING", unset = "10000"
))

required <- c(input_path, output_path, result_path)
if (any(!nzchar(required))) stop("input, store, and result paths are required")
if (!file.exists(input_path)) stop("input not found: ", input_path)
if (region_padding != 10000L || pvalue_threshold != 1e-5) {
  stop("this acceptance harness is locked to p <= 1e-5 and a 10-kb window")
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
unlink(output_path, recursive = TRUE, force = TRUE)

elapsed <- system.time({
  store <- compress_sumstats(
    input_path, output_path,
    input_build = "GRCh38", store_build = "GRCh38",
    selection = "core_plus", pvalue_threshold = pvalue_threshold,
    region_padding = region_padding, threads = threads,
    overwrite = TRUE
  )
})[["elapsed"]]

files <- list.files(output_path, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
bytes <- sum(file.info(files)$size, na.rm = TRUE)
selection <- store$manifest$selection
result <- list(
  benchmark = "core_plus_ntrk3_prepared",
  commit = commit,
  source = list(
    path = normalizePath(input_path, mustWork = TRUE),
    basename = basename(input_path),
    bytes = unname(file.info(input_path)$size),
    rows = store$manifest$source$rows %||% NA_integer_,
    sha256 = store$manifest$source$sha256 %||% NA_character_
  ),
  build = list(input = "GRCh38", stored = "GRCh38"),
  parameters = list(
    selection = "core_plus", pvalue_threshold = pvalue_threshold,
    threshold_operator = selection$threshold_operator,
    region_padding = region_padding,
    window_bp_each_side = selection$window_bp_each_side,
    window_boundary = selection$window_boundary,
    threads = threads
  ),
  result = list(
    elapsed_seconds = as.numeric(elapsed),
    output_bytes = as.numeric(bytes),
    stored_rows = as.integer(store$manifest$n_rows),
    core_variant_rows = selection$core_variant_rows,
    seed_snps = selection$seed_snps,
    regions = selection$regions,
    source_input_rows = selection$source_input_rows %||% NA_integer_,
    timings = store$manifest$timings
  ),
  manifest = list(
    profile = store$manifest$profile,
    codec = store$manifest$codec,
    format_version = store$manifest$format_version
  )
)
jsonlite::write_json(result, result_path, auto_unbox = TRUE, pretty = TRUE,
                     digits = 17)
print(result)
