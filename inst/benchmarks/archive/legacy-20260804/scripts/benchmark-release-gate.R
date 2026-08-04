#!/usr/bin/env Rscript

# Full release-gate benchmark. Run this on the Mac mini. Large inputs and
# temporary stores live on the external SSD; only the result CSV/Markdown are
# written into the project outputs directory.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
})

root <- Sys.getenv("COMPRESSOR_BENCH_ROOT", "/Volumes/crucial_x9/CompreSSoR-benchmarks")
project <- "/Users/fergushamilton/projects/CompreSSoR"
raw_dir <- file.path(root, "raw")
store_dir <- file.path(root, "stores")
cache_dir <- file.path(root, "cache", "reference")
report_dir <- file.path(root, "reports")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(store_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project, "outputs"), recursive = TRUE, showWarnings = FALSE)

source_url <- "https://storage.googleapis.com/finngen-public-data-r2/summary_stats/finngen_r2_ANTIDEPRESSANTS.gz"
source_path <- Sys.getenv("COMPRESSOR_SOURCE", file.path(raw_dir, basename(source_url)))
reference_path <- Sys.getenv(
  "COMPRESSOR_REFERENCE",
  "/Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz"
)
if (!file.exists(source_path)) {
  temporary <- tempfile("compressor-download-", tmpdir = raw_dir)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  download.file(source_url, temporary, mode = "wb", quiet = FALSE)
  if (!file.rename(temporary, source_path)) stop("could not install source GWAS: ", source_path)
}
if (!file.exists(reference_path)) stop("reference not found: ", reference_path)

reference <- list(
  id = "bp_shared_hm3",
  build = "GRCh38",
  source = "MR-atlas shared GRCh38 dictionary",
  variants = reference_path,
  cache_dir = cache_dir
)

read_head <- function(path, n = 100000L) {
  command <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    paste("gzip -dc", shQuote(normalizePath(path)))
  } else NULL
  if (is.null(command)) {
    data.table::fread(path, data.table = FALSE, nrows = n, showProgress = FALSE,
                      check.names = FALSE)
  } else {
    data.table::fread(cmd = command, data.table = FALSE, nrows = n,
                      showProgress = FALSE, check.names = FALSE)
  }
}

file_bytes <- function(path) {
  files <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE) else path
  files <- files[file.info(files)$isdir %in% FALSE & file.exists(files)]
  sum(file.info(files)$size, na.rm = TRUE)
}

records <- list()
record <- function(scenario, replicate, elapsed, rows = NA_integer_, bytes = NA_real_,
                   valid = NA, max_error = NA_real_, note = "") {
  records[[length(records) + 1L]] <<- data.frame(
    scenario = scenario, replicate = as.integer(replicate),
    elapsed_seconds = as.numeric(elapsed), rows = as.integer(rows),
    bytes = as.numeric(bytes), valid = as.logical(valid),
    max_error = as.numeric(max_error), note = note,
    stringsAsFactors = FALSE
  )
}

timed <- function(scenario, replicate, fun, rows = NA_integer_, path = NULL) {
  gc(FALSE)
  started <- proc.time()[["elapsed"]]
  value <- fun()
  elapsed <- proc.time()[["elapsed"]] - started
  valid <- if (inherits(value, "compressor_store")) validate_compressor(value)$valid else NA
  record(scenario, replicate, elapsed, rows = rows,
         bytes = if (!is.null(path)) file_bytes(path) else NA_real_, valid = valid)
  value
}

slice <- CompreSSoR:::normalise_sumstats_columns(read_head(source_path, 100000L))
slice_store_exact <- file.path(store_dir, "slice-exact.cpr")
slice_store_standard <- file.path(store_dir, "slice-standard.cpr")
slice_store_qc_serial <- file.path(store_dir, "slice-qc-serial.cpr")
slice_store_qc_parallel <- file.path(store_dir, "slice-qc-parallel.cpr")

for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_SLICE_RUNS", "5")))) {
  timed("slice_convert_exact_compress", i, function() compress_sumstats(
    slice, slice_store_exact, reference = NULL, mode = "convert",
    profile = "exact", overwrite = TRUE
  ), rows = nrow(slice), path = slice_store_exact)
}
exact_store <- open_compressor(slice_store_exact)
exact_decoded <- NULL
for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_SLICE_RUNS", "5")))) {
  exact_decoded <- timed("slice_exact_decompress", i, function() decompress_sumstats(exact_store), rows = nrow(slice))
  core <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "variant_id", "beta", "standard_error", "effect_allele_frequency", "p_value")
  same <- isTRUE(all.equal(exact_decoded[core], slice[core], check.attributes = FALSE))
  record("slice_exact_roundtrip_check", i, 0, rows = nrow(slice), valid = same)
}

for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_SLICE_RUNS", "5")))) {
  timed("slice_standard_compress", i, function() compress_sumstats(
    slice, slice_store_standard, reference = NULL, mode = "convert",
    profile = "standard", overwrite = TRUE
  ), rows = nrow(slice), path = slice_store_standard)
}
standard_store <- open_compressor(slice_store_standard)
standard_decoded <- NULL
for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_SLICE_RUNS", "5")))) {
  standard_decoded <- timed("slice_standard_decompress", i, function() decompress_sumstats(standard_store), rows = nrow(slice))
  max_error <- max(abs(standard_decoded$beta - slice$beta), na.rm = TRUE)
  record("slice_standard_roundtrip_check", i, 0, rows = nrow(slice),
         valid = nrow(standard_decoded) == nrow(slice), max_error = max_error)
}

build_cache(standard_store, overwrite = TRUE, block_rows = 4096L)
region_chromosome <- as.character(slice$chromosome[1L])
region_start <- min(slice$base_pair_location[slice$chromosome == region_chromosome], na.rm = TRUE)
region <- paste0("chr", region_chromosome, ":", region_start, "-", region_start + 1000000)
for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_REGION_RUNS", "5")))) {
  timed("slice_region_direct", i, function() read_sumstats(standard_store, region = region), rows = NA_integer_)
  timed("slice_region_q8_cache", i, function() read_sumstats(standard_store, region = region, use_cache = TRUE), rows = NA_integer_)
}

for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_SLICE_QC_RUNS", "3")))) {
  timed("slice_qc_serial_compress", i, function() compress_sumstats(
    slice, slice_store_qc_serial, reference = reference, mode = "qc",
    profile = "standard", chrom_threads = 1L, overwrite = TRUE
  ), rows = nrow(slice), path = slice_store_qc_serial)
  timed("slice_qc_chrom4_compress", i, function() compress_sumstats(
    slice, slice_store_qc_parallel, reference = reference, mode = "qc",
    profile = "standard", chrom_threads = 4L, overwrite = TRUE
  ), rows = nrow(slice), path = slice_store_qc_parallel)
}

full_store_convert <- file.path(store_dir, "finngen-convert.cpr")
full_store_qc <- file.path(store_dir, "finngen-qc.cpr")
full_runs <- as.integer(Sys.getenv("COMPRESSOR_FULL_RUNS", "5"))
for (i in seq_len(full_runs)) {
  timed("finngen_convert_standard_compress", i, function() compress_sumstats(
    source_path, full_store_convert, reference = NULL, mode = "convert",
    profile = "standard", overwrite = TRUE
  ), path = full_store_convert)
}
convert_store <- open_compressor(full_store_convert)
for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_FULL_DECOMP_RUNS", "5")))) {
  decoded <- timed("finngen_convert_standard_decompress", i, function() decompress_sumstats(convert_store),
                   rows = convert_store$manifest$n_rows, path = full_store_convert)
  core <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "variant_id", "p_value")
  same_head <- isTRUE(all.equal(decoded[seq_len(nrow(slice)), core], slice[core],
                                check.attributes = FALSE))
  record("finngen_convert_roundtrip_shape", i, 0, rows = nrow(decoded),
         valid = nrow(decoded) == convert_store$manifest$n_rows)
  record("finngen_convert_roundtrip_head", i, 0, rows = nrow(slice), valid = same_head)
}

qc_runs <- as.integer(Sys.getenv("COMPRESSOR_QC_RUNS", "3"))
for (i in seq_len(qc_runs)) {
  timed("finngen_qc_serial_compress", i, function() compress_sumstats(
    source_path, full_store_qc, reference = reference, mode = "qc",
    profile = "standard", chrom_threads = 1L, overwrite = TRUE
  ), path = full_store_qc)
}
qc_store <- open_compressor(full_store_qc)
for (i in seq_len(as.integer(Sys.getenv("COMPRESSOR_QC_DECOMP_RUNS", "3")))) {
  decoded <- timed("finngen_qc_decompress", i, function() decompress_sumstats(qc_store),
                   rows = qc_store$manifest$n_rows, path = full_store_qc)
  record("finngen_qc_roundtrip_shape", i, 0, rows = nrow(decoded),
         valid = nrow(decoded) == qc_store$manifest$n_rows)
}

result <- do.call(rbind, records)
result_stem <- Sys.getenv("COMPRESSOR_RESULT_STEM", "release-gate-benchmark")
result_path <- file.path(project, "outputs", paste0(result_stem, ".csv"))
write.csv(result, result_path, row.names = FALSE)
writeLines(c(
  "# CompreSSoR release-gate benchmark",
  "",
  paste0("Source: ", source_path),
  paste0("Reference: ", reference_path),
  paste0("External working root: ", root),
  "",
  "This run covers exact and lossy round trips, direct and q8 regional reads,", 
  "repeated conversion/decompression, and real FinnGen conversion plus QC.",
  "Large inputs and stores remain outside the synced project.",
  "",
  "```text",
  capture.output(print(result)),
  "```",
  "",
  paste0("R version: ", R.version.string),
  paste0("data.table threads: ", data.table::getDTthreads())
), file.path(project, "outputs", paste0(result_stem, ".md")))
saveRDS(list(source = source_path, reference = reference_path, root = root,
             results = result, session = sessionInfo()),
        file.path(report_dir, paste0(result_stem, ".rds")))
message("Wrote ", result_path)
