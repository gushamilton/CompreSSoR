#!/usr/bin/env Rscript

# Benchmark optional parallel Pcodec reads without copying the GWAS into the
# repository. Set COMPRESSOR_THREAD_STORE to a .cpr directory and
# COMPRESSOR_THREAD_SOURCE to the matching TSV.gz used to choose keys.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
})

store_path <- Sys.getenv("COMPRESSOR_THREAD_STORE", unset = "")
source_path <- Sys.getenv("COMPRESSOR_THREAD_SOURCE", unset = "")
result_root <- Sys.getenv("COMPRESSOR_THREAD_RESULTS",
                          file.path(getwd(), "outputs", "pcodec-threads"))
runs <- as.integer(Sys.getenv("COMPRESSOR_THREAD_RUNS", unset = "5"))
thread_values <- as.integer(strsplit(
  Sys.getenv("COMPRESSOR_THREAD_VALUES", unset = "1,2,4"), ",", fixed = TRUE
)[[1]])

if (!dir.exists(store_path)) stop("set COMPRESSOR_THREAD_STORE to a .cpr directory")
if (!file.exists(source_path)) stop("set COMPRESSOR_THREAD_SOURCE to the matching TSV.gz")
if (!is.finite(runs) || runs < 5L) stop("COMPRESSOR_THREAD_RUNS must be at least 5")
if (any(!is.finite(thread_values) | thread_values < 1L)) stop("thread values must be positive integers")
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

read_source <- function() {
  x <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
             data.table = TRUE, showProgress = FALSE)
  names(x)[names(x) == "#chrom"] <- "chrom"
  x
}

source <- read_source()
keys <- compressor_variant_key(source$chrom, source$pos, source$ref, source$alt)
keys_25 <- keys[unique(as.integer(round(seq(1L, nrow(source), length.out = 25L))))]
keys_1000 <- keys[unique(as.integer(round(seq(1L, nrow(source), length.out = 1000L))))]
core_columns <- c("chromosome", "base_pair_location", "z",
                  "standard_error", "effect_allele_frequency")
stores_5 <- rep(store_path, 5L)
names(stores_5) <- paste0("gwas", seq_along(stores_5))

checksum <- function(value) {
  sum(vapply(value, function(x) sum(as.numeric(x), na.rm = TRUE), numeric(1)))
}

records <- list()
record <- function(threads, workload, run, elapsed, value) {
  records[[length(records) + 1L]] <<- data.table(
    threads = as.integer(threads), workload = workload, run = as.integer(run),
    seconds = as.numeric(elapsed), rows = as.integer(nrow(value)),
    checksum = checksum(value)
  )
}

time_read <- function(expr) {
  gc(FALSE)
  value <- NULL
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(elapsed = elapsed, value = value)
}

for (threads in thread_values) {
  for (run in seq_len(runs)) {
    measured <- time_read(read_sumstats(store_path, columns = core_columns,
                                        threads = threads))
    record(threads, "full_core", run, measured$elapsed, measured$value)

    measured <- time_read(read_sumstats(store_path,
                                        region = "chr1:100000000-101000000",
                                        columns = core_columns, threads = threads))
    record(threads, "region_chr1_1mb", run, measured$elapsed, measured$value)

    measured <- time_read(read_sumstats(store_path, variants = keys_25,
                                        columns = core_columns, threads = threads))
    record(threads, "sparse_25_keys", run, measured$elapsed, measured$value)

    measured <- time_read(read_sumstats(store_path, variants = keys_1000,
                                        columns = core_columns, threads = threads))
    record(threads, "sparse_1000_keys", run, measured$elapsed, measured$value)

    measured <- time_read(read_sumstats_batch(
      stores_5, variants = keys_25, columns = core_columns, threads = threads
    ))
    record(threads, "batch_5_stores_25_keys", run, measured$elapsed,
           do.call(rbind, measured$value))
  }
}

runs_table <- rbindlist(records)
summary <- runs_table[, .(
  runs = .N, median_seconds = median(seconds), min_seconds = min(seconds),
  max_seconds = max(seconds), rows = max(rows), checksum_median = median(checksum)
), by = .(threads, workload)]
fwrite(runs_table, file.path(result_root, "pcodec-threads-runs.csv"))
fwrite(summary, file.path(result_root, "pcodec-threads-summary.csv"))
message("Wrote ", file.path(result_root, "pcodec-threads-summary.csv"))
