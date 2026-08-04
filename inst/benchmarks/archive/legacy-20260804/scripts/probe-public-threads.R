#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(CompreSSoR)
})

store_path <- Sys.getenv("COMPRESSOR_PROBE_STORE")
if (!nzchar(store_path)) stop("COMPRESSOR_PROBE_STORE is required")
columns <- c("base_pair_location", "effect_allele", "other_allele", "z",
             "standard_error", "effect_allele_frequency")
checksum <- function(x) sum(as.numeric(x$base_pair_location),
                             as.numeric(x$z), as.numeric(x$standard_error),
                             as.numeric(x$effect_allele_frequency), na.rm = TRUE)
rows <- rbindlist(lapply(c(1L, 2L, 4L, 8L), function(n_threads) {
  times <- numeric(3L)
  checks <- numeric(3L)
  for (i in seq_along(times)) {
    started <- proc.time()[["elapsed"]]
    value <- read_sumstats(store_path, columns = columns, threads = n_threads)
    times[i] <- proc.time()[["elapsed"]] - started
    checks[i] <- checksum(value)
  }
  data.table(threads = n_threads, median_seconds = median(times),
             min_seconds = min(times), max_seconds = max(times),
             rows = nrow(value), checksum_min = min(checks),
             checksum_max = max(checks), checksum_match = length(unique(checks)) == 1L)
}))
fwrite(rows, Sys.getenv("COMPRESSOR_PROBE_OUTPUT", unset = "public-thread-probe.csv"))
print(rows)
