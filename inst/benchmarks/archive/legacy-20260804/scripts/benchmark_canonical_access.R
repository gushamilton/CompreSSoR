#!/usr/bin/env Rscript

# Five-repeat access benchmark for the canonical full-FinnGen production store.
# Large stores and bridges remain on the external volume; only compact records
# are written to inst/benchmarks.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

store_path <- Sys.getenv(
  "COMPRESSOR_CANONICAL_BENCH_STORE",
  "/Volumes/crucial_x9/CompreSSoR-benchmarks/finngen-full-pcodec-v03-release.cpr"
)
output_dir <- Sys.getenv(
  "COMPRESSOR_CANONICAL_BENCH_OUTPUT",
  file.path(getwd(), "inst", "benchmarks")
)
runs <- as.integer(Sys.getenv("COMPRESSOR_CANONICAL_BENCH_RUNS", "5"))
stopifnot(dir.exists(store_path), runs >= 5L)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

store <- open_compressor(store_path)
n <- store$manifest$n_rows
identity_columns <- c(
  "chromosome", "base_pair_location", "effect_allele", "other_allele"
)
row_ids <- unique(as.integer(floor(seq(0, n - 1, length.out = 25L))))
identity <- read_sumstats(store, variants = row_ids, columns = identity_columns)
keys <- compressor_variant_key(
  identity$chromosome, identity$base_pair_location,
  identity$other_allele, identity$effect_allele
)

core_columns <- c(identity_columns, "z", "standard_error", "effect_allele_frequency")
derived_columns <- c(core_columns, "beta", "p_value")
workloads <- list(
  canonical_25_beta_se = function() read_sumstats(
    store, variants = keys,
    columns = c(identity_columns, "beta", "standard_error")
  ),
  region_chr1_1mb = function() read_sumstats(
    store, region = "chr1:100000000-101000000", columns = core_columns
  ),
  full_core = function() read_sumstats(store, columns = core_columns),
  full_with_beta_p = function() read_sumstats(store, columns = derived_columns)
)

touch <- function(data) {
  numeric <- vapply(data, is.numeric, logical(1))
  sum(vapply(data[numeric], function(column) sum(column, na.rm = TRUE), numeric(1)))
}

records <- list()
for (workload in names(workloads)) {
  invisible(workloads[[workload]]())
  for (run in seq_len(runs)) {
    gc()
    value <- NULL
    elapsed <- system.time(value <- workloads[[workload]]())[["elapsed"]]
    records[[length(records) + 1L]] <- data.table(
      workload = workload,
      run = run,
      seconds = as.numeric(elapsed),
      rows = nrow(value),
      checksum = touch(value)
    )
    message(sprintf("%s run %d/%d: %.3f s", workload, run, runs, elapsed))
  }
}
records <- rbindlist(records)
stopifnot(records[, uniqueN(checksum), by = workload][, all(V1 == 1L)])
summary <- records[, .(
  runs = .N,
  median_seconds = median(seconds),
  min_seconds = min(seconds),
  max_seconds = max(seconds),
  rows = max(rows)
), by = workload]
files <- list.files(store_path, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, no.. = TRUE)
store_bytes <- sum(file.info(files)$size)
fwrite(records, file.path(output_dir, "pcodec-canonical-access-runs.csv"))
write_json(list(
  schema_version = "1.0.0",
  measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_version = as.character(utils::packageVersion("CompreSSoR")),
  source_revision = Sys.getenv(
    "COMPRESSOR_BENCH_COMMIT",
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
  ),
  machine = Sys.info()[c("nodename", "sysname", "release", "machine")],
  store = list(
    rows = n,
    bytes = unname(store_bytes),
    format_version = store$manifest$format_version,
    chunk_rows = store$manifest$wrapped_codec$chunk_rows,
    page_rows = store$manifest$wrapped_codec$page_rows,
    index_version = store$manifest$wrapped_codec$index_version,
    key_block_rows = store$manifest$key_block_rows,
    value_block_rows = store$manifest$value_block_rows,
    first_key = keys[[1L]]
  ),
  summary = summary,
  runs = records
), file.path(output_dir, "pcodec-canonical-access.json"),
auto_unbox = TRUE, pretty = TRUE)
print(summary)
