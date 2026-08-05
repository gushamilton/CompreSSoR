#!/usr/bin/env Rscript

# Synthetic sizing/access benchmark for the experimental core-plus store.
# The generated GWAS and stores live outside the synced project; only compact
# CSV/Markdown summaries are written under outputs/.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

# For local development, source the working tree into an environment whose
# parent is the installed package namespace. This lets the harness exercise
# uninstalled R changes while reusing a compatible installed native library.
if (identical(Sys.getenv("COMPRESSOR_BENCH_SOURCE_TREE", unset = "0"), "1")) {
  source_root <- normalizePath(Sys.getenv("COMPRESSOR_SOURCE_ROOT", unset = "."),
                               mustWork = TRUE)
  source_env <- new.env(parent = asNamespace("CompreSSoR"))
  source_files <- list.files(file.path(source_root, "R"), pattern = "[.]R$",
                             full.names = TRUE)
  for (source_file in sort(source_files)) sys.source(source_file, source_env)
  compress_sumstats <- source_env$compress_sumstats
  open_compressor <- source_env$open_compressor
  read_sumstats <- source_env$read_sumstats
}

rows <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_ROWS", unset = "200000"))
signals <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_SIGNALS", unset = "10"))
repeats <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_RUNS", unset = "3"))
backend <- match.arg(Sys.getenv("COMPRESSOR_PVALUE_BACKEND", unset = "pcodec"),
                     c("pcodec", "parquet"))
if (backend == "parquet") options(CompreSSoR.native_decode = FALSE)
bench_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_BENCH_ROOT",
  unset = file.path(tempdir(), "compressor-pvalue-core-regions")
)
result_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_RESULTS",
  unset = file.path(getwd(), "outputs")
)
if (rows < 1000L || signals < 1L || repeats < 1L) {
  stop("rows must be >= 1000, signals >= 1, and repeats >= 1")
}
dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  as.numeric(sum(file.info(files)$size, na.rm = TRUE))
}

make_synthetic <- function(n, signal_count) {
  position <- seq.int(1000001L, by = 100L, length.out = n)
  signal_rows <- unique(round(seq.int(1L, n, length.out = signal_count)))
  z <- 0.8 * sin(seq_len(n) / 97)
  z[signal_rows] <- 5
  se <- 0.02 + (seq_len(n) %% 31) / 10000
  data.frame(
    chromosome = rep("1", n),
    base_pair_location = position,
    effect_allele = rep("A", n),
    other_allele = rep("C", n),
    beta = z * se,
    standard_error = se,
    effect_allele_frequency = 0.1 + (seq_len(n) %% 80) / 100,
    z = z,
    stringsAsFactors = FALSE
  )
}

input <- make_synthetic(rows, signals)
reference <- input[c("chromosome", "base_pair_location",
                     "effect_allele", "other_allele")]
core_rows <- unique(round(seq.int(max(1L, floor(rows / (signals * 2L))),
                                  rows, length.out = signals)))
core_panel <- input[core_rows, c("chromosome", "base_pair_location",
                                 "effect_allele", "other_allele")]
full_path <- file.path(bench_root, "synthetic-full.cpr")
core_path <- file.path(bench_root, "synthetic-core-pvalue-regions.cpr")
unlink(c(full_path, core_path), recursive = TRUE, force = TRUE)

write_time <- function(path, mode) {
  reference_arg <- if (mode == "core_plus") reference else NULL
  panel_arg <- if (mode == "core_plus") core_panel else NULL
  timer <- system.time({
    compress_sumstats(
      input, path, reference = reference_arg, mode = mode,
      variant_set = panel_arg,
      backend = backend, assume_grch38_ref_alt = TRUE,
      pvalue_threshold = 1e-5, region_padding = 50000L,
      overwrite = TRUE
    )
  })[["elapsed"]]
  as.numeric(timer)
}

full_write <- write_time(full_path, "convert")
core_write <- write_time(core_path, "core_plus")
full <- open_compressor(full_path)
core <- open_compressor(core_path)

read_runs <- function(store) {
  vapply(seq_len(repeats), function(i) {
    timer <- system.time({
      value <- read_sumstats(store, columns = c("z", "standard_error"))
    })[["elapsed"]]
    as.numeric(timer)
  }, numeric(1))
}

full_read <- read_runs(full)
core_read <- read_runs(core)
selection <- core$manifest$selection
full_bytes <- directory_bytes(full_path)
core_bytes <- directory_bytes(core_path)
records <- data.table(
  backend = backend,
  store = c("full", "core_plus_pvalue_regions"),
  source_rows = rows,
  stored_rows = c(full$manifest$n_rows, core$manifest$n_rows),
  stored_fraction = c(1, core$manifest$n_rows / rows),
  core_variant_rows = c(NA_integer_, selection$core_variant_rows),
  seed_snps = c(NA_integer_, selection$seed_snps),
  regions = c(NA_integer_, selection$regions),
  padding_bp = c(NA_integer_, selection$padding_bp),
  bytes = c(full_bytes, core_bytes),
  bytes_per_stored_row = c(full_bytes / full$manifest$n_rows,
                           core_bytes / core$manifest$n_rows),
  write_seconds = c(full_write, core_write),
  median_read_seconds = c(median(full_read), median(core_read)),
  core_as_fraction_of_full_bytes = c(1, core_bytes / full_bytes)
)

csv_path <- file.path(result_root, "pvalue-core-regions-benchmark.csv")
md_path <- file.path(result_root, "pvalue-core-regions-benchmark.md")
fwrite(records, csv_path)

lines <- c(
  "# Synthetic core-plus storage benchmark",
  "",
  sprintf("This is a synthetic chr1-only benchmark using the %s backend; it is intentionally not a FinnGen run.", backend),
  sprintf("Input rows: %s; injected p < 1e-5 signals: %s; padding: 50,000 bp.",
          format(rows, big.mark = ","), signals),
  "P-values are selected from the canonical Z score and are derived rather than stored.",
  "",
  "| Store | Rows | Row fraction | Regions | Bytes | Write (s) | Median full read (s) |",
  "|---|---:|---:|---:|---:|---:|---:|",
  sprintf("| full | %s | 100.0%% | — | %s | %.3f | %.3f |",
          format(full$manifest$n_rows, big.mark = ","),
          format(full_bytes, big.mark = ","), full_write, median(full_read)),
  sprintf("| core_plus_pvalue_regions | %s | %.2f%% | %s | %s | %.3f | %.3f |",
          format(core$manifest$n_rows, big.mark = ","),
          100 * core$manifest$n_rows / rows, selection$regions,
          format(core_bytes, big.mark = ","), core_write, median(core_read)),
  "",
  sprintf("Core-plus store size as a fraction of the full store: %.2f%%.",
          100 * core_bytes / full_bytes),
  sprintf("Region sidecar size: %s bytes.",
          file.info(file.path(core_path, selection$file))$size),
  "",
  "The benchmark measures storage and whole-store decode only; it does not measure regional reads."
)
writeLines(lines, md_path, useBytes = TRUE)
print(records)
cat("Wrote", csv_path, "and", md_path, "\n")
