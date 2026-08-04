#!/usr/bin/env Rscript

# Final release benchmark for the exported CompreSSoR API. Large inputs,
# temporary bridges, and the durable store remain on the caller-supplied
# external volume; only compact CSV/JSON records are written to the project.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

source_path <- Sys.getenv("COMPRESSOR_BENCH_SOURCE", unset = "")
bench_root <- Sys.getenv("COMPRESSOR_BENCH_ROOT", unset = "")
result_root <- Sys.getenv(
  "COMPRESSOR_BENCH_RESULTS",
  file.path(getwd(), "outputs", "release-api")
)
python <- Sys.getenv("COMPRESSOR_PYTHON", unset = "")
runs <- as.integer(Sys.getenv("COMPRESSOR_BENCH_RUNS", unset = "5"))
public_labels <- !identical(Sys.getenv("COMPRESSOR_BENCH_PUBLIC", unset = "1"), "0")

if (!nzchar(source_path) || !nzchar(bench_root)) {
  stop("set COMPRESSOR_BENCH_SOURCE and COMPRESSOR_BENCH_ROOT to external-storage paths")
}
stopifnot(file.exists(source_path), runs >= 5L, nzchar(python), file.exists(python))
dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
bridge_root <- file.path(bench_root, "tmp")
dir.create(bridge_root, recursive = TRUE, showWarnings = FALSE)
options(CompreSSoR.python = python, CompreSSoR.tempdir = bridge_root)

store_path <- file.path(bench_root, "finngen-full-pcodec.cpr")
csv_path <- file.path(result_root, "pcodec-full-api-runs.csv")
json_path <- file.path(result_root, "pcodec-full-api-benchmark.json")
roundtrip_path <- file.path(result_root, "pcodec-full-api-roundtrip.json")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "scripts", "benchmark_package_api_final.R"),
                mustWork = TRUE)
}
backend_path <- system.file("python", "compressor_pcodec.py", package = "CompreSSoR")
native_path <- file.path(getwd(), "src", "codec_native.cpp")
git_sha <- tryCatch(trimws(system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
                    error = function(e) NA_character_)
git_dirty <- tryCatch(length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0L,
                      error = function(e) NA)
source_sha256 <- digest::digest(source_path, algo = "sha256", file = TRUE)
implementation <- list(
  package = as.character(utils::packageVersion("CompreSSoR")),
  git_sha = git_sha,
  git_worktree_dirty = git_dirty,
  benchmark_script_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE),
  pcodec_backend_sha256 = digest::digest(backend_path, algo = "sha256", file = TRUE),
  native_reader_sha256 = digest::digest(native_path, algo = "sha256", file = TRUE),
  source_sha256 = source_sha256
)
benchmark_signature <- digest::digest(implementation, algo = "sha256", serialize = TRUE)
records <- if (file.exists(csv_path)) {
  existing <- fread(csv_path)
  if (!"benchmark_signature" %in% names(existing) ||
      any(existing$benchmark_signature != benchmark_signature)) {
    stop("existing benchmark records belong to a different implementation; use a fresh COMPRESSOR_BENCH_RESULTS directory")
  }
  lapply(seq_len(nrow(existing)), function(i) existing[i])
} else {
  list()
}

is_complete <- function(format_name, workload, run) {
  if (!length(records)) return(FALSE)
  table <- rbindlist(records, fill = TRUE)
  any(table$format == format_name & table$workload == workload & table$run == run)
}

directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                      include.dirs = FALSE, no.. = TRUE)
  as.numeric(sum(file.info(files)$size, na.rm = TRUE))
}

save_records <- function() {
  table <- rbindlist(records, fill = TRUE)
  fwrite(table, csv_path)
  summaries <- table[, .(
    repeats = .N,
    median_seconds = median(seconds),
    min_seconds = min(seconds),
    max_seconds = max(seconds),
    rows = max(rows),
    bytes = max(bytes)
  ), by = .(format, workload)]
  payload <- list(
    schema_version = "1.0.0",
    measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    machine = Sys.info()[c("nodename", "sysname", "release", "machine")],
    R = R.version.string,
    data_table_threads = getDTthreads(),
    python = trimws(system2(python, c("-c", shQuote(
      "import importlib.metadata as m, platform, numpy; print(f'Python {platform.python_version()}; numpy {numpy.__version__}; pcodec {m.version(\"pcodec\")}; zstandard {m.version(\"zstandard\")}')"
    )), stdout = TRUE)),
    runs_required = runs,
    benchmark_signature = benchmark_signature,
    implementation = implementation,
    source = list(
      file = if (public_labels) "FinnGen SNP core.tsv.gz on external SSD" else source_path,
      bytes = unname(file.info(source_path)$size),
      sha256 = source_sha256,
      columns = c("chrom", "pos", "a1", "a2", "beta", "se", "eaf", "p")
    ),
    store = list(
      path = if (public_labels) "FinnGen full Pcodec store on external SSD" else store_path,
      bytes = if (dir.exists(store_path)) directory_bytes(store_path) else NA_real_,
      manifest_sha256 = if (file.exists(file.path(store_path, "manifest.json"))) {
        digest::digest(file.path(store_path, "manifest.json"), algo = "sha256", file = TRUE)
      } else NA_character_,
      self_contained_identity = TRUE,
      representation = "lossless GRCh38 position+REF+ALT key; Z9/EAF8/SE6 Pcodec streams"
    ),
    summaries = summaries,
    runs = table
  )
  write_json(payload, json_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
}

record_run <- function(format_name, workload, run, expression, bytes = NA_real_) {
  gc()
  value <- NULL
  elapsed <- system.time(value <- force(expression))[["elapsed"]]
  rows <- as.numeric(value$rows)
  checksum <- as.numeric(value$checksum)
  records[[length(records) + 1L]] <<- data.table(
    format = format_name,
    workload = workload,
    run = as.integer(run),
    seconds = as.numeric(elapsed),
    rows = rows,
    checksum = checksum,
    bytes = as.numeric(bytes),
    benchmark_signature = benchmark_signature
  )
  save_records()
  message(sprintf("[%s] %s run %d/%d: %.3f s (%s rows)",
                  format_name, workload, run, runs, elapsed,
                  format(rows, scientific = FALSE, big.mark = ",")))
  invisible(value)
}

touch_pcodec <- function(columns, region = NULL, variants = NULL) {
  x <- read_sumstats(store_path, columns = columns, region = region, variants = variants)
  numeric_names <- intersect(c("base_pair_location", "z", "standard_error",
                               "effect_allele_frequency", "beta", "p_value"), names(x))
  checksum <- sum(vapply(x[numeric_names], function(column) {
    sum(as.numeric(column), na.rm = TRUE)
  }, numeric(1)))
  list(rows = nrow(x), checksum = checksum)
}

read_source <- function(include_p = FALSE, region = NULL, sparse_rows = NULL) {
  selected <- c("chrom", "pos", "a1", "a2", "beta", "se", "eaf")
  if (include_p) selected <- c(selected, "p")
  x <- fread(source_path, select = selected, showProgress = FALSE)
  if (!is.null(region)) {
    keep <- x$chrom == region$chrom & x$pos >= region$start & x$pos <= region$end
    x <- x[keep]
  }
  if (!is.null(sparse_rows)) x <- x[sparse_rows + 1L]
  z <- x$beta / x$se
  checksum <- sum(x$pos, z, x$se, x$eaf, na.rm = TRUE)
  if (include_p) checksum <- checksum + sum(x$p, na.rm = TRUE)
  list(rows = nrow(x), checksum = checksum)
}

source_bytes <- unname(file.info(source_path)$size)

# Compression is intentionally end-to-end: import, validation, R-to-Python
# bridge, encoding, checksums, manifest finalisation, and atomic replacement.
for (i in seq_len(runs)) {
  if (is_complete("CompreSSoR Pcodec", "compress_end_to_end", i)) next
  result <- record_run("CompreSSoR Pcodec", "compress_end_to_end", i, {
    store <- compress_sumstats(
      source_path,
      store_path,
      reference = NULL,
      mode = "convert",
      backend = "pcodec",
      profile = "standard",
      assume_grch38_ref_alt = TRUE,
      overwrite = dir.exists(store_path)
    )
    list(rows = store$manifest$n_rows, checksum = store$manifest$n_rows)
  })
  records[[length(records)]]$bytes <- directory_bytes(store_path)
  save_records()
}

manifest <- open_compressor(store_path)$manifest
n_rows <- as.integer(manifest$n_rows)
sparse_rows <- unique(as.integer(floor(seq(0, n_rows - 1, length.out = 1000L))))
region <- list(chrom = "1", start = 1e6L, end = 2e6L)
core_columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                  "z", "standard_error", "effect_allele_frequency")
derived_columns <- c(core_columns, "beta", "p_value")
store_bytes <- directory_bytes(store_path)

for (i in seq_len(runs)) {
  if (!is_complete("CompreSSoR Pcodec", "full_core", i))
    record_run("CompreSSoR Pcodec", "full_core", i,
               touch_pcodec(core_columns), bytes = store_bytes)
  if (!is_complete("CompreSSoR Pcodec", "full_with_beta_p", i))
    record_run("CompreSSoR Pcodec", "full_with_beta_p", i,
               touch_pcodec(derived_columns), bytes = store_bytes)
  if (!is_complete("CompreSSoR Pcodec", "region_chr1_1mb", i))
    record_run("CompreSSoR Pcodec", "region_chr1_1mb", i,
               touch_pcodec(core_columns, region = "chr1:1000000-2000000"), bytes = store_bytes)
  if (!is_complete("CompreSSoR Pcodec", "sparse_1000_rows", i))
    record_run("CompreSSoR Pcodec", "sparse_1000_rows", i,
               touch_pcodec(core_columns, variants = sparse_rows), bytes = store_bytes)
}

for (i in seq_len(runs)) {
  if (!is_complete("TSV gzip", "full_core", i))
    record_run("TSV gzip", "full_core", i,
               read_source(include_p = FALSE), bytes = source_bytes)
  if (!is_complete("TSV gzip", "full_with_beta_p", i))
    record_run("TSV gzip", "full_with_beta_p", i,
               read_source(include_p = TRUE), bytes = source_bytes)
  if (!is_complete("TSV gzip", "region_chr1_1mb", i))
    record_run("TSV gzip", "region_chr1_1mb", i,
               read_source(include_p = FALSE, region = region), bytes = source_bytes)
  if (!is_complete("TSV gzip", "sparse_1000_rows", i))
    record_run("TSV gzip", "sparse_1000_rows", i,
               read_source(include_p = FALSE, sparse_rows = sparse_rows), bytes = source_bytes)
}

for (i in seq_len(runs)) {
  if (!is_complete("CompreSSoR Pcodec", "validate_checksums", i))
    record_run("CompreSSoR Pcodec", "validate_checksums", i, {
      result <- validate_compressor(store_path, full = FALSE)
      stopifnot(isTRUE(result$valid))
      list(rows = result$rows, checksum = result$files_checked)
    }, bytes = store_bytes)
  if (!is_complete("CompreSSoR Pcodec", "validate_all_frames", i))
    record_run("CompreSSoR Pcodec", "validate_all_frames", i, {
      result <- validate_compressor(store_path, full = TRUE)
      stopifnot(isTRUE(result$valid))
      list(rows = result$rows, checksum = result$frames_checked)
    }, bytes = store_bytes)
}

# A real-data round-trip audit compares 10,000 evenly spaced rows after sorting
# by the format's exact self-contained identity key.
source <- fread(source_path, select = c("chrom", "pos", "a1", "a2", "beta", "se", "eaf"),
                showProgress = FALSE)
chrom_lengths <- c(
  `1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
  `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
  `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
  `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
  `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
  `21` = 46709983, `22` = 50818468, X = 156040895, Y = 57227415
)
offsets <- c(0, head(cumsum(chrom_lengths), -1L))
names(offsets) <- names(chrom_lengths)
base_code <- c(A = 0, C = 1, G = 2, T = 3)
key <- (unname(offsets[as.character(source$chrom)]) + source$pos - 1) * 16 +
  unname(base_code[source$a2]) * 4 + unname(base_code[source$a1])
if (is.unsorted(key, strictly = TRUE)) {
  source_order <- order(key, method = "radix")
} else {
  source_order <- seq_len(nrow(source))
}
audit_rows <- unique(as.integer(floor(seq(0, n_rows - 1, length.out = 10000L))))
expected <- source[source_order[audit_rows + 1L]]
observed <- read_sumstats(store_path, variants = audit_rows, columns = derived_columns)
expected_z <- expected$beta / expected$se
identity_ok <- identical(as.character(observed$chromosome), as.character(expected$chrom)) &&
  identical(as.integer(observed$base_pair_location), as.integer(expected$pos)) &&
  identical(as.character(observed$effect_allele), as.character(expected$a1)) &&
  identical(as.character(observed$other_allele), as.character(expected$a2))
roundtrip <- list(
  benchmark_signature = benchmark_signature,
  implementation = implementation,
  rows_checked = length(audit_rows),
  identity_exact = identity_ok,
  max_abs_z_error = max(abs(observed$z - expected_z), na.rm = TRUE),
  max_abs_eaf_error = max(abs(observed$effect_allele_frequency - expected$eaf), na.rm = TRUE),
  max_relative_se_error = max(abs(observed$standard_error / expected$se - 1), na.rm = TRUE),
  max_abs_derived_beta_error = max(abs(observed$beta - expected$beta), na.rm = TRUE),
  p_is_exactly_derived_from_decoded_z = isTRUE(all.equal(
    observed$p_value, 2 * pnorm(-abs(observed$z)), tolerance = 1e-14
  ))
)
stopifnot(
  isTRUE(roundtrip$identity_exact),
  roundtrip$max_abs_z_error <= 0.007,
  roundtrip$max_abs_eaf_error <= 0.0041,
  roundtrip$max_relative_se_error <= 0.012,
  isTRUE(roundtrip$p_is_exactly_derived_from_decoded_z)
)
write_json(roundtrip, roundtrip_path, auto_unbox = TRUE, pretty = TRUE)
save_records()
message("Final benchmark and real-data round-trip audit completed successfully")
