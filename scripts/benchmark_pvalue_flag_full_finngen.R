#!/usr/bin/env Rscript

# Full-file comparison of reconstructed-p filtering versus an aligned binary
# significance flag. Raw input, the full store, flag streams, and logs stay
# outside the synced project; only compact summaries are retained.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

source_root <- normalizePath(
  Sys.getenv("COMPRESSOR_SOURCE_ROOT", unset = "."), mustWork = TRUE
)
source_env <- new.env(parent = asNamespace("CompreSSoR"))
source_files <- list.files(file.path(source_root, "R"), pattern = "[.]R$",
                           full.names = TRUE)
for (source_file in sort(source_files)) sys.source(source_file, source_env)

compress_sumstats <- source_env$compress_sumstats
open_compressor <- source_env$open_compressor
read_sumstats <- source_env$read_sumstats
pcodec_native_append_stream <- source_env$pcodec_native_append_stream
pcodec_native_decompress <- source_env$pcodec_native_decompress

input_path <- Sys.getenv(
  "COMPRESSOR_PVALUE_FULL_INPUT",
  unset = "/Volumes/crucial_x9/CompreSSoR-benchmarks/raw/finngen_r2_ANTIDEPRESSANTS.gz"
)
bench_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_FULL_BENCH_ROOT",
  unset = file.path(tempdir(), "compressor-pvalue-full-finngen")
)
result_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_FULL_RESULTS",
  unset = file.path(getwd(), "inst/benchmarks/pvalue-sidechannel-macmini-20260806")
)
thresholds <- as.numeric(strsplit(
  Sys.getenv("COMPRESSOR_PVALUE_FULL_THRESHOLDS", unset = "1e-4,1e-5,1e-6,1e-7,1e-8"),
  ",", fixed = TRUE
)[[1L]])
repeats <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_FULL_RUNS", unset = "5"))
threads <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_FULL_THREADS", unset = "4"))
rebuild <- identical(Sys.getenv("COMPRESSOR_PVALUE_FULL_REBUILD", unset = "0"), "1")
build <- Sys.getenv("COMPRESSOR_PVALUE_FULL_BUILD", unset = "GRCh38")
dataset_label <- Sys.getenv(
  "COMPRESSOR_PVALUE_FULL_LABEL", unset = basename(input_path)
)

# This is the downstream MR projection: identity, beta, SE, and EAF. The
# p-filter and flag paths both fetch exactly this same projection after
# identifying the selected row IDs.
mr_columns <- c(
  "chromosome", "base_pair_location", "effect_allele", "other_allele",
  "beta", "standard_error", "effect_allele_frequency"
)

if (!file.exists(input_path)) stop("FinnGen input not found: ", input_path)
if (!length(thresholds) || any(!is.finite(thresholds) | thresholds <= 0 |
                                thresholds > 1) || repeats < 1L || threads < 1L) {
  stop("invalid thresholds, repeats, or threads")
}

dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
store_path <- file.path(bench_root, "finngen-full-current.cpr")
flag_root <- file.path(bench_root, "flags")
dir.create(flag_root, recursive = TRUE, showWarnings = FALSE)

column <- function(data, primary, alternatives = character()) {
  candidates <- c(primary, alternatives)
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) stop("FinnGen input is missing column: ", primary)
  data[[hit[[1L]]]]
}

read_command <- if (grepl("[.]gz$|[.]bgz$", input_path, ignore.case = TRUE)) {
  paste("gzip -dc", shQuote(normalizePath(input_path)))
} else {
  paste("cat", shQuote(normalizePath(input_path)))
}
raw <- fread(cmd = read_command, data.table = TRUE, showProgress = FALSE)
prepared <- data.table(
  chromosome = as.character(column(raw, "#chrom", c("chromosome", "chr"))),
  position = suppressWarnings(as.integer(column(raw, "pos", c("base_pair_location", "position")))),
  ref = toupper(as.character(column(raw, "ref", c("other_allele", "other allele", "REF")))),
  alt = toupper(as.character(column(raw, "alt", c("effect_allele", "effect allele", "ALT")))),
  beta = suppressWarnings(as.numeric(column(raw, "beta"))),
  se = suppressWarnings(as.numeric(column(raw, "sebeta", c("standard_error", "se")))),
  eaf = suppressWarnings(as.numeric(column(raw, "maf", c("eaf", "effect_allele_frequency")))),
  p_value = suppressWarnings(as.numeric(column(raw, "pval", c("p_value", "p-value"))))
)
rm(raw)
gc(FALSE)

valid <- !is.na(prepared$chromosome) & is.finite(prepared$position) &
  !is.na(prepared$ref) & !is.na(prepared$alt) &
  nchar(prepared$ref) == 1L & nchar(prepared$alt) == 1L &
  prepared$ref %in% c("A", "C", "G", "T") &
  prepared$alt %in% c("A", "C", "G", "T") & prepared$ref != prepared$alt &
  is.finite(prepared$beta) & is.finite(prepared$se) & prepared$se > 0 &
  is.finite(prepared$eaf) & prepared$eaf >= 0 & prepared$eaf <= 1
prepared <- prepared[valid]
rm(valid)
gc(FALSE)

duplicate <- duplicated(prepared, by = c("chromosome", "position", "ref", "alt"))
if (any(duplicate)) prepared <- prepared[!duplicate]
rm(duplicate)
prepared[, chromosome_order := match(
  chromosome, c(as.character(seq_len(22L)), "X", "Y")
)]
setorder(prepared, chromosome_order, position, ref, alt)

source_p <- prepared$p_value
sumstats <- data.frame(
  chromosome = prepared$chromosome,
  base_pair_location = prepared$position,
  reference_allele = prepared$ref,
  alternate_allele = prepared$alt,
  effect_allele = prepared$alt,
  other_allele = prepared$ref,
  beta = prepared$beta,
  standard_error = prepared$se,
  effect_allele_frequency = prepared$eaf,
  z = prepared$beta / prepared$se,
  stringsAsFactors = FALSE
)
rm(prepared)
gc(FALSE)
n <- nrow(sumstats)

directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  if (!length(files)) return(0)
  as.numeric(sum(file.info(files)$size, na.rm = TRUE))
}

write_full_store <- function(path) {
  started <- proc.time()[["elapsed"]]
  store <- compress_sumstats(
    sumstats, path, input_build = build, store_build = build,
    selection = "full", backend = "pcodec", qc = "none", threads = threads,
    overwrite = TRUE, pvalue_flag = FALSE
  )
  list(store = store, elapsed_seconds = unname(proc.time()[["elapsed"]] - started))
}

existing_store <- FALSE
if (dir.exists(store_path) && !rebuild) {
  manifest <- tryCatch(
    read_json(file.path(store_path, "manifest.json"), simplifyVector = FALSE),
    error = function(e) NULL
  )
  manifest_rows <- if (!is.null(manifest$n_rows)) manifest$n_rows else
    if (!is.null(manifest$rows)) manifest$rows else NA
  existing_store <- !is.null(manifest) &&
    identical(as.character(manifest$format_version), "0.4.5-pcodec-native") &&
    identical(as.integer(manifest_rows), as.integer(n))
}
if (!existing_store) unlink(store_path, recursive = TRUE, force = TRUE)
full_result <- if (existing_store) {
  list(store = open_compressor(store_path), elapsed_seconds = NA_real_)
} else {
  write_full_store(store_path)
}
full_store <- full_result$store
if (!identical(as.integer(full_store$manifest$n_rows), as.integer(n))) {
  stop("full store row count does not match prepared input")
}
store_bytes <- directory_bytes(store_path)

read_stats <- function() {
  timer <- system.time({
    value <- read_sumstats(
      full_store, columns = c("z", "standard_error"), threads = threads
    )
    observed <- nrow(value)
  })[["elapsed"]]
  if (!identical(as.integer(observed), as.integer(n))) {
    stop("ordinary read row count failed")
  }
  as.numeric(timer)
}

read_p_and_filter <- function(threshold, fetch = FALSE) {
  timer <- system.time({
    value <- read_sumstats(full_store, columns = "p_value", threads = threads)
    p <- as.numeric(value$p_value)
    row_ids <- which(is.finite(p) & p <= threshold) - 1L
    if (fetch) {
      fetched <- read_sumstats(
        full_store, variants = row_ids,
        columns = mr_columns, threads = 1L
      )
      fetched_rows <- nrow(fetched)
      if (!identical(sort(names(fetched)), sort(mr_columns))) {
        stop("MR projection did not return the requested columns")
      }
    } else {
      fetched_rows <- 0L
    }
  })[["elapsed"]]
  list(seconds = as.numeric(timer), rows = length(row_ids),
       fetched_rows = as.integer(fetched_rows))
}

read_flag_row_ids <- function(path, info) {
  blob <- readBin(path, raw(), n = file.info(path)$size)
  pieces <- lapply(info$blocks, function(block) {
    first <- as.integer(block$offset) + 1L
    last <- first + as.integer(block$length) - 1L
    pcodec_native_decompress(blob[first:last], as.integer(block$values), "u8")
  })
  flags <- unlist(pieces, use.names = FALSE)
  which(flags != 0L) - 1L
}

ordinary_before <- median(vapply(seq_len(repeats), function(i) read_stats(), numeric(1)))

results <- lapply(thresholds, function(threshold) {
  source_rows <- which(is.finite(source_p) & source_p <= threshold)
  label <- paste0("p-", format(threshold, scientific = TRUE))
  flag_dir <- file.path(flag_root, label)
  dir.create(flag_dir, recursive = TRUE, showWarnings = FALSE)
  flag_path <- file.path(flag_dir, "significant.pco")

  flag_info <- NULL
  flag_write_runs <- numeric(repeats)
  for (i in seq_len(repeats)) {
    flag <- integer(n)
    flag[source_rows] <- 1L
    started <- proc.time()[["elapsed"]]
    flag_info <- pcodec_native_append_stream(
      flag, flag_path, dtype = "u8", workers = threads
    )
    flag_write_runs[[i]] <- unname(proc.time()[["elapsed"]] - started)
  }
  write_json(list(
    codec = "Pcodec uint8",
    semantic = "aligned_binary_significance_flag",
    threshold = threshold,
    rows = n,
    source_hit_rows = length(source_rows),
    stream = flag_info
  ), file.path(flag_dir, "flag.json"), auto_unbox = TRUE,
  pretty = TRUE, digits = 17)
  flag_bytes <- directory_bytes(flag_dir)

  p_filter_runs <- lapply(seq_len(repeats), function(i) {
    read_p_and_filter(threshold, fetch = FALSE)
  })
  p_filter_fetch_runs <- lapply(seq_len(repeats), function(i) {
    read_p_and_filter(threshold, fetch = TRUE)
  })
  flag_only_runs <- lapply(seq_len(repeats), function(i) {
    started <- proc.time()[["elapsed"]]
    row_ids <- read_flag_row_ids(flag_path, flag_info)
    elapsed <- unname(proc.time()[["elapsed"]] - started)
    list(seconds = elapsed, rows = length(row_ids))
  })
  flag_fetch_runs <- lapply(seq_len(repeats), function(i) {
    started <- proc.time()[["elapsed"]]
    row_ids <- read_flag_row_ids(flag_path, flag_info)
    fetched <- read_sumstats(
      full_store, variants = row_ids,
      columns = mr_columns, threads = 1L
    )
    if (!identical(sort(names(fetched)), sort(mr_columns))) {
      stop("MR projection did not return the requested columns")
    }
    elapsed <- unname(proc.time()[["elapsed"]] - started)
    list(seconds = elapsed, rows = length(row_ids), fetched_rows = nrow(fetched))
  })

  p_filter_rows <- unique(vapply(p_filter_runs, function(x) x$rows, integer(1)))
  p_filter_fetch_rows <- unique(vapply(
    p_filter_fetch_runs, function(x) x$rows, integer(1)
  ))
  flag_rows <- unique(vapply(flag_only_runs, function(x) x$rows, integer(1)))
  flag_fetch_rows <- unique(vapply(flag_fetch_runs, function(x) x$rows, integer(1)))
  if (length(p_filter_rows) != 1L || length(p_filter_fetch_rows) != 1L ||
      length(flag_rows) != 1L || length(flag_fetch_rows) != 1L ||
      any(flag_rows != length(source_rows)) ||
      any(flag_fetch_rows != length(source_rows))) {
    stop("full-file row selection varied or did not match source flags")
  }

  data.table(
    threshold = threshold,
    source_hit_rows = length(source_rows),
    p_filter_hit_rows = p_filter_rows,
    flag_hit_rows = flag_rows,
    store_rows = n,
    store_bytes = store_bytes,
    flag_bytes = flag_bytes,
    flag_write_seconds = median(flag_write_runs),
    p_filter_seconds = median(vapply(
      p_filter_runs, function(x) x$seconds, numeric(1)
    )),
    p_filter_fetch_seconds = median(vapply(
      p_filter_fetch_runs, function(x) x$seconds, numeric(1)
    )),
    flag_only_seconds = median(vapply(
      flag_only_runs, function(x) x$seconds, numeric(1)
    )),
    flag_fetch_seconds = median(vapply(
      flag_fetch_runs, function(x) x$seconds, numeric(1)
    )),
    p_filter_fetch_rows = p_filter_fetch_rows,
    flag_fetch_rows = flag_fetch_rows,
    flag_additional_fraction = flag_bytes / store_bytes,
    ordinary_read_before_seconds = ordinary_before
  )
})

ordinary_after <- median(vapply(seq_len(repeats), function(i) read_stats(), numeric(1)))
records <- rbindlist(results, fill = TRUE)
records[, ordinary_read_after_seconds := ordinary_after]
records[, ordinary_read_delta_seconds := ordinary_after - ordinary_read_before_seconds]

csv_path <- file.path(result_root, "pvalue-sidechannel-full-finngen-summary.csv")
md_path <- file.path(result_root, "pvalue-sidechannel-full-finngen-summary.md")
json_path <- file.path(result_root, "pvalue-sidechannel-full-finngen-provenance.json")
fwrite(records, csv_path)

lines <- c(
  "# Full-file p-value versus flag-read benchmark",
  "",
  sprintf("Dataset: %s.", dataset_label),
  sprintf("Full valid-SNP store: %s rows; store bytes: %s; native threads: %d.",
          format(n, big.mark = ","), format(round(store_bytes), big.mark = ","), threads),
  sprintf("Ordinary z/SE read before flag streams: %.4f s; after: %.4f s; delta: %.4f s.",
          ordinary_before, ordinary_after, ordinary_after - ordinary_before),
  "P-filter reads reconstructed p from the full store and filters in R.",
  paste0("The two MR paths fetch the same columns: ", paste(mr_columns, collapse = ", "), "."),
  "Flag inspection decodes the full aligned uint8 flag stream; flag + MR fetch then reads selected rows by native row ID.",
  "P filtering + MR fetch reads reconstructed p for the whole store, filters in R, then reads the same selected MR columns.",
  "",
  "| Threshold | Source hits | P-filter hits | Flag bytes | P-filter + MR fetch s | Flag inspection + MR fetch s |",
  "|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(records)), function(i) {
    row <- records[i]
    sprintf("| p <= %s | %s | %s | %s | %.4f | %.4f |",
            format(row$threshold, scientific = TRUE),
            format(round(row$source_hit_rows), big.mark = ","),
            format(round(row$p_filter_hit_rows), big.mark = ","),
            format(round(row$flag_bytes), big.mark = ","),
            row$p_filter_fetch_seconds, row$flag_fetch_seconds)
  }, character(1)),
  "",
  "The flag bytes include the Pcodec uint8 payload plus its small metadata/index.",
  "The p-filter path is the fastest current whole-file API path for this question; p is reconstructed because it is not stored in the native payload.",
  ""
)
writeLines(lines, md_path, useBytes = TRUE)

source_commit <- tryCatch(
  system2("git", c("-C", source_root, "rev-parse", "HEAD"),
          stdout = TRUE, stderr = FALSE)[[1L]],
  error = function(e) NA_character_
)
provenance <- list(
  benchmark = "pvalue_flag_full_finngen",
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source = normalizePath(input_path),
  store = normalizePath(store_path),
  rows = n,
  thresholds = thresholds,
  dataset_label = dataset_label,
  build = build,
  mr_columns = mr_columns,
  repeats = repeats,
  threads = threads,
  rebuild = rebuild,
  store_bytes = store_bytes,
  build_seconds = full_result$elapsed_seconds,
  source_commit = source_commit,
  r = R.version.string,
  platform = Sys.info()[["sysname"]],
  architecture = R.version$arch
)
write_json(provenance, json_path, auto_unbox = TRUE, pretty = TRUE, digits = 17)

print(records)
cat("Wrote", csv_path, "and", md_path, "\n")
