#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
})

source_path <- Sys.getenv("COMPRESSOR_FINAL_SOURCE", unset = "")
store_path <- Sys.getenv("COMPRESSOR_PCODEC_STORE", unset = "")
result_path <- Sys.getenv("COMPRESSOR_PCODEC_RESULT", unset = "")
runs <- as.integer(Sys.getenv("COMPRESSOR_PCODEC_RUNS", unset = "5"))
if (!file.exists(source_path) || !nzchar(store_path) || !nzchar(result_path)) {
  stop("source, store and result paths are required")
}
if (!is.finite(runs) || runs < 5L) stop("runs must be at least five")

source_data <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
                     showProgress = FALSE)
setnames(source_data, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p"))
source_data[, chrom := toupper(sub("^CHR", "", as.character(chrom)))]
source_data[, z := beta / se]
if (nrow(source_data) != 10000000L) stop("expected exactly 10m source rows")

dir.create(dirname(store_path), recursive = TRUE, showWarnings = FALSE)
store <- compress_sumstats(
  source_data[, .(chromosome = chrom, base_pair_location = pos,
                  effect_allele = alt, other_allele = ref, beta,
                  standard_error = se, effect_allele_frequency = eaf, z)],
  store_path, reference = NULL, mode = "convert", backend = "pcodec",
  profile = "standard", assume_grch38_ref_alt = TRUE, overwrite = TRUE
)
stopifnot(isTRUE(validate_compressor(store, full = TRUE)$valid))

columns <- c("global_position", "substitution", "z", "standard_error",
             "effect_allele_frequency")
checksum <- function(x) sum(vapply(x, function(value) sum(as.numeric(value), na.rm = TRUE), 0.0))
records <- list()
for (threads in c(1L, 4L)) {
  for (run in seq_len(runs)) {
    gc(FALSE)
    value <- NULL
    elapsed <- system.time(value <- read_sumstats(
      store_path, columns = columns, threads = threads))[['elapsed']]
    if (nrow(value) != nrow(source_data)) stop("row count mismatch")
    records[[length(records) + 1L]] <- data.table(
      threads = threads, run = run, seconds = as.numeric(elapsed),
      rows = nrow(value), checksum = checksum(value),
      storage_bytes = sum(file.info(list.files(store_path, full.names = TRUE,
                                               recursive = TRUE))$size, na.rm = TRUE))
  }
}
records <- rbindlist(records)
dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
fwrite(records, result_path)
summary <- records[, .(runs = .N, median_seconds = median(seconds),
                       min_seconds = min(seconds), max_seconds = max(seconds),
                       rows = max(rows), storage_bytes = max(storage_bytes)),
                   by = threads]
fwrite(summary, sub("\\.csv$", "-summary.csv", result_path))
cat("Wrote", result_path, "\n")
print(summary)
