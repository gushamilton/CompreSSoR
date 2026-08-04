#!/usr/bin/env Rscript

# Repeated chr1 benchmark for the current native CompreSSoR Pcodec store.
# The input is a prefiltered FinnGen chr1 SNV TSV.gz. Large input/output
# stores stay on the external benchmark volume; only CSV/JSON evidence is
# written to the requested result directory.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

source_path <- Sys.getenv("COMPRESSOR_CHR1_SOURCE", unset = "")
bench_root <- Sys.getenv("COMPRESSOR_BENCH_ROOT", unset = "")
result_root <- Sys.getenv("COMPRESSOR_BENCH_RESULTS", file.path(getwd(), "outputs", "finngen-chr1"))
runs <- as.integer(Sys.getenv("COMPRESSOR_BENCH_RUNS", unset = "5"))
if (!nzchar(source_path) || !file.exists(source_path)) {
  stop("set COMPRESSOR_CHR1_SOURCE to an existing FinnGen chr1 SNV TSV.gz")
}
if (!nzchar(bench_root)) stop("set COMPRESSOR_BENCH_ROOT to external benchmark storage")
if (!is.finite(runs) || runs < 5L) stop("COMPRESSOR_BENCH_RUNS must be at least 5")

dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
store_path <- file.path(bench_root, "finngen-chr1-native.cpr")
run_path <- file.path(result_root, "finngen-chr1-native-runs.csv")
summary_path <- file.path(result_root, "finngen-chr1-native.csv")
json_path <- file.path(result_root, "finngen-chr1-native.json")

directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  if (!length(files)) return(0)
  sum(file.info(files)$size, na.rm = TRUE)
}

read_source <- function() {
  command <- paste("gzip -dc", shQuote(normalizePath(source_path)))
  x <- fread(cmd = command, data.table = TRUE, showProgress = FALSE)
  names(x)[names(x) == "#chrom"] <- "chrom"
  names(x)[names(x) == "sebeta"] <- "se"
  names(x)[names(x) == "maf"] <- "eaf"
  x
}

checksum <- function(x, derived = FALSE) {
  numeric_columns <- intersect(c("pos", "beta", "se", "eaf"), names(x))
  value <- sum(vapply(x[, numeric_columns, with = FALSE],
                      function(z) sum(as.numeric(z), na.rm = TRUE), numeric(1)))
  if (derived) {
    z <- x$beta / x$se
    value <- value + sum(z, na.rm = TRUE) + sum(2 * pnorm(-abs(z)), na.rm = TRUE)
  }
  value
}

checksum_store <- function(x, derived = FALSE) {
  numeric_columns <- intersect(c("base_pair_location", "z", "standard_error",
                                 "effect_allele_frequency", "beta", "p_value"), names(x))
  sum(vapply(x[numeric_columns], function(z) sum(as.numeric(z), na.rm = TRUE), numeric(1)))
}

source_probe <- read_source()
if (!nrow(source_probe)) stop("chr1 source has no rows")
source_keys <- compressor_variant_key(source_probe$chrom, source_probe$pos,
                                      source_probe$ref, source_probe$alt)
sparse_rows <- unique(as.integer(round(seq(1, nrow(source_probe), length.out = 1000L))))
sparse_keys <- source_keys[sparse_rows]
region <- "chr1:100000000-101000000"
core_columns <- c("chromosome", "base_pair_location", "z",
                  "standard_error", "effect_allele_frequency")
derived_columns <- c(core_columns, "beta", "p_value")
source_bytes <- unname(file.info(source_path)$size)
source_rows <- nrow(source_probe)

records <- list()
record <- function(format, workload, run, elapsed, rows, bytes, value) {
  records[[length(records) + 1L]] <<- data.table(
    format = format, workload = workload, run = as.integer(run),
    seconds = as.numeric(elapsed), rows = as.integer(rows),
    bytes = as.numeric(bytes), checksum = as.numeric(value)
  )
  fwrite(rbindlist(records), run_path)
  message(sprintf("[%s] %s %d/%d: %.3f s (%s rows)", format, workload, run,
                  runs, elapsed, format(rows, big.mark = ",")))
}

time_read <- function(expr) {
  gc(FALSE)
  value <- NULL
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(elapsed = elapsed, value = value)
}

# End-to-end compression is measured separately from access. Each repetition
# imports the gzip input and replaces the same external store.
for (i in seq_len(runs)) {
  measured <- time_read(compress_sumstats(
    source_path, store_path, reference = NULL, mode = "convert",
    backend = "pcodec", profile = "standard", overwrite = dir.exists(store_path),
    assume_grch38_ref_alt = TRUE
  ))
  record("CompreSSoR Pcodec", "compress_end_to_end", i, measured$elapsed,
         measured$value$manifest$n_rows, directory_bytes(store_path),
         measured$value$manifest$n_rows)
}

store <- open_compressor(store_path)
store_bytes <- directory_bytes(store_path)
for (i in seq_len(runs)) {
  measured <- time_read(read_sumstats(store, columns = core_columns))
  record("CompreSSoR Pcodec", "full_core", i, measured$elapsed,
         nrow(measured$value), store_bytes, checksum_store(measured$value))
  measured <- time_read(read_sumstats(store, columns = derived_columns))
  record("CompreSSoR Pcodec", "full_with_beta_p", i, measured$elapsed,
         nrow(measured$value), store_bytes, checksum_store(measured$value))
  measured <- time_read(read_sumstats(store, region = region, columns = core_columns))
  record("CompreSSoR Pcodec", "region_chr1_1mb", i, measured$elapsed,
         nrow(measured$value), store_bytes, checksum_store(measured$value))
  measured <- time_read(read_sumstats(store, variants = sparse_keys, columns = core_columns))
  record("CompreSSoR Pcodec", "sparse_1000_keys", i, measured$elapsed,
         nrow(measured$value), store_bytes, checksum_store(measured$value))
}

# TSV.gz is read end-to-end for every access workload. This is the fair
# unindexed baseline: regional and sparse requests must still inflate the file.
for (i in seq_len(runs)) {
  measured <- time_read({ x <- read_source(); list(x = x, derived = FALSE) })
  record("TSV gzip", "full_core", i, measured$elapsed,
         nrow(measured$value$x), source_bytes, checksum(measured$value$x))
  measured <- time_read({ x <- read_source(); list(x = x, derived = TRUE) })
  record("TSV gzip", "full_with_beta_p", i, measured$elapsed,
         nrow(measured$value$x), source_bytes, checksum(measured$value$x, derived = TRUE))
  measured <- time_read({ x <- read_source(); x <- x[pos >= 100000000 & pos <= 101000000]; x })
  record("TSV gzip", "region_chr1_1mb", i, measured$elapsed,
         nrow(measured$value), source_bytes, checksum(measured$value))
  measured <- time_read({ x <- read_source(); x[sparse_rows] })
  record("TSV gzip", "sparse_1000_keys", i, measured$elapsed,
         nrow(measured$value), source_bytes, checksum(measured$value))
}

runs_table <- rbindlist(records)
summary <- runs_table[, .(
  runs = .N, median_seconds = median(seconds), min_seconds = min(seconds),
  max_seconds = max(seconds), rows = max(rows), bytes = max(bytes),
  checksum_median = median(checksum)
), by = .(format, workload)]
fwrite(runs_table, run_path)
fwrite(summary, summary_path)
write_json(list(
  benchmark_id = "compressor_finngen_chr1_native",
  measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_version = as.character(utils::packageVersion("CompreSSoR")),
  git_sha = tryCatch(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
                     error = function(e) NA_character_),
  source = list(path = source_path, rows = source_rows, bytes = source_bytes,
                sha256 = digest::digest(source_path, algo = "sha256", file = TRUE),
                description = "FinnGen chr1 biallelic SNVs; external identity not counted"),
  store = list(path = store_path, rows = store$manifest$n_rows, bytes = store_bytes,
               format = store$manifest$format_version,
               codec = store$manifest$codec$name),
  runs = runs, summary = summary, records = runs_table
), json_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
message("Wrote ", summary_path)
