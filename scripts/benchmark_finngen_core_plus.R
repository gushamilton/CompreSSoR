#!/usr/bin/env Rscript

# Real FinnGen chr1 storage benchmark for the prepared-input core_plus path.
# Any reference-backed harmonisation/liftover must happen before this script.
# Raw input and generated stores stay outside the synced project; only compact
# CSV/Markdown summaries are written under outputs/.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

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
  read_variant_set <- get("read_variant_set", source_env)
  variant_set_membership <- get("variant_set_membership", source_env)
} else {
  compress_sumstats <- CompreSSoR::compress_sumstats
  open_compressor <- CompreSSoR::open_compressor
  read_sumstats <- CompreSSoR::read_sumstats
  read_variant_set <- CompreSSoR:::read_variant_set
  variant_set_membership <- CompreSSoR:::variant_set_membership
}

input_path <- Sys.getenv(
  "COMPRESSOR_FINNGEN_INPUT",
  unset = "/Volumes/crucial_x9/CompreSSoR-benchmarks/chr1-native-final-20260804/finngen_chr1_snvs.tsv.gz"
)
if (!file.exists(input_path)) stop("FinnGen input not found: ", input_path)

core_panel_path <- Sys.getenv(
  "COMPRESSOR_FINNGEN_CORE_PANEL",
  unset = "/Volumes/crucial_x9/CompreSSoR-benchmarks/prepared-core-hm3-panels-20260805/core.tsv.gz"
)
if (!file.exists(core_panel_path)) {
  stop("prepared core panel not found: ", core_panel_path,
       "; run scripts/prepare_core_hm3_panels.R on the mini first", call. = FALSE)
}

bench_root <- Sys.getenv(
  "COMPRESSOR_FINNGEN_BENCH_ROOT",
  unset = "/Volumes/crucial_x9/CompreSSoR-benchmarks/finngen-core-plus-20260805"
)
result_root <- Sys.getenv(
  "COMPRESSOR_FINNGEN_RESULTS",
  unset = file.path(getwd(), "outputs")
)
repeats <- as.integer(Sys.getenv("COMPRESSOR_FINNGEN_RUNS", unset = "3"))
compression_threads <- as.integer(Sys.getenv("COMPRESSOR_FINNGEN_THREADS", unset = "4"))
pvalue_threshold <- as.numeric(Sys.getenv("COMPRESSOR_FINNGEN_PVALUE_THRESHOLD", unset = "1e-5"))
region_padding <- as.integer(Sys.getenv("COMPRESSOR_FINNGEN_REGION_PADDING", unset = "10000"))
if (repeats < 1L) stop("repeats must be positive")
if (compression_threads < 1L) stop("COMPRESSOR_FINNGEN_THREADS must be positive")

dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

preparation_time <- system.time({
  source_table <- fread(
    cmd = paste("gzip -dc", shQuote(normalizePath(input_path))),
    data.table = TRUE, showProgress = FALSE
  )
  column <- function(primary, alternatives = character()) {
    hit <- c(primary, alternatives)[c(primary, alternatives) %in% names(source_table)]
    if (!length(hit)) stop("FinnGen input is missing column: ", primary)
    source_table[[hit[[1L]]]]
  }
  prepared_effect <- toupper(as.character(column("alt", "effect_allele")))
  prepared_other <- toupper(as.character(column("ref", "other_allele")))
  sumstats <- data.frame(
    chromosome = as.character(column("#chrom", "chromosome")),
    base_pair_location = as.integer(column("pos", "base_pair_location")),
    reference_allele = prepared_other,
    alternate_allele = prepared_effect,
    effect_allele = prepared_effect,
    other_allele = prepared_other,
    beta = as.numeric(column("beta")),
    standard_error = as.numeric(column("sebeta", "standard_error")),
    effect_allele_frequency = as.numeric(column("maf", "eaf")),
    stringsAsFactors = FALSE
  )
  sumstats$z <- sumstats$beta / sumstats$standard_error
})[["elapsed"]]
preparation_time <- as.numeric(preparation_time)

panel_read_time <- system.time({
  benchmark_panel <- read_variant_set(
    core_panel_path, chromosomes = unique(sumstats$chromosome)
  )
})[["elapsed"]]
panel_read_time <- as.numeric(panel_read_time)
panel_match_time <- system.time({
  benchmark_panel_keep <- variant_set_membership(sumstats, benchmark_panel)
})[["elapsed"]]
panel_match_time <- as.numeric(panel_match_time)

full_path <- file.path(bench_root, "finngen-chr1-full.cpr")
core_plus_path <- file.path(bench_root, "finngen-chr1-core-plus.cpr")
unlink(c(full_path, core_plus_path), recursive = TRUE, force = TRUE)

write_time <- function(path, selection_scope) {
  timer <- system.time({
    if (identical(selection_scope, "full")) {
      compress_sumstats(
        sumstats, path, input_build = "GRCh38", store_build = "GRCh38",
        selection = "full", backend = "pcodec", threads = compression_threads,
        overwrite = TRUE
      )
    } else {
      compress_sumstats(
        sumstats, path, input_build = "GRCh38", store_build = "GRCh38",
        selection = "core_plus", variant_set = core_panel_path,
        backend = "pcodec", threads = compression_threads,
        pvalue_threshold = pvalue_threshold,
        region_padding = region_padding, overwrite = TRUE
      )
    }
  })[["elapsed"]]
  as.numeric(timer)
}

full_write <- write_time(full_path, "full")
core_plus_write <- write_time(core_plus_path, "core_plus")
full <- open_compressor(full_path)
core_plus <- open_compressor(core_plus_path)

directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  as.numeric(sum(file.info(files)$size, na.rm = TRUE))
}

read_runs <- function(store) {
  vapply(seq_len(repeats), function(i) {
    timer <- system.time({
      value <- read_sumstats(store, columns = c("z", "standard_error"))
    })[["elapsed"]]
    as.numeric(timer)
  }, numeric(1))
}

full_read <- read_runs(full)
core_plus_read <- read_runs(core_plus)
selection <- core_plus$manifest$selection
full_bytes <- directory_bytes(full_path)
core_plus_bytes <- directory_bytes(core_plus_path)

records <- data.table(
  backend = "pcodec",
  input = basename(input_path),
  source_rows = nrow(sumstats),
  source_bytes = as.numeric(file.info(input_path)$size),
  store = c("full", "core_plus"),
  stored_rows = c(full$manifest$n_rows, core_plus$manifest$n_rows),
  stored_fraction = c(1, core_plus$manifest$n_rows / nrow(sumstats)),
  core_panel_rows = c(NA_integer_, core_plus$manifest$selection$core_variant_rows),
  seed_snps = c(NA_integer_, selection$seed_snps),
  regions = c(NA_integer_, selection$regions),
  padding_bp = c(NA_integer_, selection$padding_bp),
  bytes = c(full_bytes, core_plus_bytes),
  write_seconds = c(full_write, core_plus_write),
  median_read_seconds = c(median(full_read), median(core_plus_read)),
  core_plus_as_fraction_of_full_bytes = c(1, core_plus_bytes / full_bytes),
  preparation_seconds = preparation_time,
  panel_read_seconds = panel_read_time,
  panel_match_seconds = panel_match_time,
  compression_threads = compression_threads,
  panel = basename(core_panel_path)
)

csv_path <- file.path(result_root, "finngen-core-plus-benchmark.csv")
md_path <- file.path(result_root, "finngen-core-plus-benchmark.md")
fwrite(records, csv_path)

sidecar <- file.info(file.path(core_plus_path, selection$file))$size
lines <- c(
  "# FinnGen chr1 core-plus storage benchmark",
  "",
  "This uses the real FinnGen chr1 biallelic-SNV input; it is not a full-genome FinnGen run.",
  sprintf("Input rows: %s; input gzip: %s bytes; staged core panel: %s rows.",
          format(nrow(sumstats), big.mark = ","),
          format(file.info(input_path)$size, big.mark = ","),
          format(core_plus$manifest$selection$core_variant_rows, big.mark = ",")),
  sprintf("P-value threshold: %s; padding: %s bp; p-values are derived from Z.",
          format(pvalue_threshold, scientific = TRUE), format(region_padding, big.mark = ",")),
  sprintf("Preparation: %.3f s; panel read: %.3f s; canonical membership: %.3f s; compression threads: %d.",
          preparation_time, panel_read_time, panel_match_time, compression_threads),
  paste0("Core panel: `", core_panel_path, "` (read against prepared identity columns; identity columns only)."),
  "Reference lookup, harmonisation, and liftover are external to this prepared-input benchmark.",
  "",
  "| Store | Rows | Row fraction | Regions | Bytes | Write (s) | Median full read (s) |",
  "|---|---:|---:|---:|---:|---:|---:|",
  sprintf("| full | %s | 100.0%% | — | %s | %.3f | %.3f |",
          format(full$manifest$n_rows, big.mark = ","),
          format(full_bytes, big.mark = ","), full_write, median(full_read)),
  sprintf("| core_plus | %s | %.2f%% | %s | %s | %.3f | %.3f |",
          format(core_plus$manifest$n_rows, big.mark = ","),
          100 * core_plus$manifest$n_rows / nrow(sumstats), selection$regions,
          format(core_plus_bytes, big.mark = ","), core_plus_write,
          median(core_plus_read)),
  "",
  sprintf("Core-plus store size as a fraction of full: %.2f%%.",
          100 * core_plus_bytes / full_bytes),
  sprintf("Region sidecar size: %s bytes.", format(sidecar, big.mark = ",")),
  "",
  "The benchmark measures storage and whole-store decode only; regional reads are not included."
)
writeLines(lines, md_path, useBytes = TRUE)
print(records)
cat("Wrote", csv_path, "and", md_path, "\n")
