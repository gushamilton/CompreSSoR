#!/usr/bin/env Rscript

# Deterministic real-store stress suite for the wrapped Pcodec reader.  Large
# data remain external; this script writes only a compact JSON result record.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(jsonlite)
})

store_path <- Sys.getenv(
  "COMPRESSOR_STRESS_STORE",
  "/Volumes/crucial_x9/CompreSSoR-benchmarks/finngen-full-pcodec-wrapped.cpr"
)
output_path <- Sys.getenv(
  "COMPRESSOR_STRESS_OUTPUT",
  file.path(getwd(), "inst", "benchmarks", "pcodec-v03-stress.json")
)
stopifnot(dir.exists(store_path))
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

store <- open_compressor(store_path)
n <- as.integer(store$manifest$n_rows)
identity_columns <- c(
  "chromosome", "base_pair_location", "effect_allele", "other_allele"
)
core_columns <- c(
  identity_columns, "z", "standard_error", "effect_allele_frequency"
)

set.seed(3801)
row_ids <- sort(sample.int(n, 2500L) - 1L)
row_batches <- split(row_ids, rep(seq_len(100L), each = 25L))
row_values <- lapply(row_batches, function(rows) {
  read_sumstats(store, variants = rows, columns = core_columns)
})
key_batches <- lapply(row_values, function(value) {
  compressor_variant_key(
    value$chromosome, value$base_pair_location,
    value$other_allele, value$effect_allele
  )
})
invisible(read_sumstats(store, variants = key_batches[[1L]], columns = core_columns))
sparse_seconds <- vapply(seq_along(key_batches), function(index) {
  value <- NULL
  elapsed <- system.time(
    value <- read_sumstats(
      store, variants = key_batches[[index]], columns = core_columns
    )
  )[["elapsed"]]
  stopifnot(identical(value, row_values[[index]]))
  as.numeric(elapsed)
}, numeric(1))

chromosome_lengths <- unlist(
  store$manifest$identity$chromosome_lengths, use.names = TRUE
)
by_row <- do.call(rbind, row_values)
region_rows <- by_row[seq.int(1L, nrow(by_row), length.out = 50L), , drop = FALSE]
region_seconds <- vapply(seq_len(nrow(region_rows)), function(index) {
  chromosome <- as.character(region_rows$chromosome[[index]])
  position <- as.integer(region_rows$base_pair_location[[index]])
  start <- max(1L, position - 500000L)
  end <- min(as.integer(chromosome_lengths[[chromosome]]), position + 500000L)
  region <- sprintf("chr%s:%d-%d", chromosome, start, end)
  value <- NULL
  elapsed <- system.time(
    value <- read_sumstats(store, region = region, columns = core_columns)
  )[["elapsed"]]
  stopifnot(nrow(value) >= 1L)
  as.numeric(elapsed)
}, numeric(1))

missing <- read_sumstats(
  store, variants = "1:2:A:C", columns = core_columns
)
duplicated <- read_sumstats(
  store, variants = rep(key_batches[[1L]][[1L]], 3L), columns = core_columns
)
stopifnot(nrow(missing) == 0L, nrow(duplicated) == 1L)

summarise <- function(seconds) list(
  runs = length(seconds),
  median_seconds = unname(median(seconds)),
  p95_seconds = unname(quantile(seconds, 0.95)),
  max_seconds = unname(max(seconds))
)

write_json(list(
  schema_version = "1.0.0",
  measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_version = as.character(utils::packageVersion("CompreSSoR")),
  source_revision = Sys.getenv(
    "COMPRESSOR_BENCH_COMMIT",
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
  ),
  machine = Sys.info()[c("nodename", "sysname", "release", "machine")],
  store = list(rows = n, format_version = store$manifest$format_version),
  parity = list(
    sampled_rows = length(row_ids), exact_row_vs_key = TRUE
  ),
  sparse_25_random = summarise(sparse_seconds),
  random_1mb_regions = summarise(region_seconds),
  missing_key_rows = nrow(missing),
  duplicate_key_rows = nrow(duplicated)
), output_path, auto_unbox = TRUE, pretty = TRUE)

print(list(
  sparse_25_random = summarise(sparse_seconds),
  random_1mb_regions = summarise(region_seconds)
))
