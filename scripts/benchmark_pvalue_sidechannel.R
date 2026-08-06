#!/usr/bin/env Rscript

# Mac mini design benchmark for p-value association side-channels.
# Raw input, stores, and logs stay outside the synced project; only compact
# summaries are retained under inst/benchmarks/.

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
  "COMPRESSOR_PVALUE_SIDECHANNEL_INPUT",
  unset = "/Volumes/crucial_x9/CompreSSoR-benchmarks/raw/finngen_r2_ANTIDEPRESSANTS.gz"
)
bench_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_SIDECHANNEL_BENCH_ROOT",
  unset = file.path(tempdir(), "compressor-pvalue-sidechannel")
)
result_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_SIDECHANNEL_RESULTS",
  unset = file.path(getwd(), "inst/benchmarks/pvalue-sidechannel-macmini-20260806")
)
rows <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_SIDECHANNEL_ROWS", unset = "100000"))
thresholds <- as.numeric(strsplit(
  Sys.getenv("COMPRESSOR_PVALUE_SIDECHANNEL_THRESHOLDS",
             unset = "1e-7,1e-6,1e-5"),
  ",", fixed = TRUE
)[[1L]])
filter_threshold <- 1e-5
synthetic_hit_count <- as.integer(Sys.getenv(
  "COMPRESSOR_PVALUE_SIDECHANNEL_SYNTHETIC_HITS", unset = "10000"
))
repeats <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_SIDECHANNEL_RUNS", unset = "5"))
threads <- as.integer(Sys.getenv("COMPRESSOR_PVALUE_SIDECHANNEL_THREADS", unset = "4"))
max_threshold <- max(thresholds)

if (!file.exists(input_path)) stop("FinnGen input not found: ", input_path)
if (rows < 1000L || repeats < 1L || threads < 1L ||
    synthetic_hit_count < 1L || synthetic_hit_count > rows ||
    !length(thresholds) || any(!is.finite(thresholds) | thresholds <= 0 |
                                thresholds > 1)) {
  stop("invalid rows, repeats, threads, or p-value thresholds")
}

dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
run_root <- file.path(bench_root, "runs")
unlink(run_root, recursive = TRUE, force = TRUE)
dir.create(run_root, recursive = TRUE, showWarnings = FALSE)

# Emit the first requested number of valid SNPs and every valid SNP at the
# largest threshold. This keeps the R fixture small while retaining the real
# low-p rows needed for the threshold sweep.
awk_program <- paste(
  "BEGIN { first = 0 }",
  "NR == 1 { print; next }",
  "{ ref = $3; alt = $4; p = $7 + 0;",
  "  valid = length(ref) == 1 && length(alt) == 1 &&",
  "    ref ~ /^[ACGT]$/ && alt ~ /^[ACGT]$/ && ref != alt;",
  "  if (valid && first < n) { print; first++; next }",
  "  if (valid && p <= maxp) print",
  "}"
)
fixture_cmd <- paste(
  "gzip -dc", shQuote(normalizePath(input_path)), "| awk -F '\\t'",
  "-v", paste0("n=", rows), "-v", paste0("maxp=", format(max_threshold, scientific = TRUE)),
  shQuote(awk_program)
)
raw <- fread(cmd = fixture_cmd, data.table = FALSE, showProgress = FALSE)

column <- function(data, primary, alternatives = character()) {
  candidates <- c(primary, alternatives)
  hit <- candidates[candidates %in% names(data)]
  if (!length(hit)) stop("FinnGen input is missing column: ", primary)
  data[[hit[[1L]]]]
}

ref <- toupper(as.character(column(raw, "ref", "other_allele")))
alt <- toupper(as.character(column(raw, "alt", "effect_allele")))
beta <- suppressWarnings(as.numeric(column(raw, "beta")))
se <- suppressWarnings(as.numeric(column(raw, "sebeta", "standard_error")))
eaf <- suppressWarnings(as.numeric(column(raw, "maf", "eaf")))
p_value <- suppressWarnings(as.numeric(column(raw, "pval", "p_value")))
chromosome <- as.character(column(raw, "#chrom", "chromosome"))
position <- suppressWarnings(as.integer(column(raw, "pos", "base_pair_location")))
valid <- !is.na(chromosome) & is.finite(position) &
  !is.na(ref) & !is.na(alt) & nchar(ref) == 1L & nchar(alt) == 1L &
  ref %in% c("A", "C", "G", "T") & alt %in% c("A", "C", "G", "T") &
  ref != alt & is.finite(beta) & is.finite(se) & se > 0 &
  is.finite(eaf) & eaf >= 0 & eaf <= 1
raw <- raw[valid, , drop = FALSE]
ref <- ref[valid]
alt <- alt[valid]
beta <- beta[valid]
se <- se[valid]
eaf <- eaf[valid]
p_value <- p_value[valid]
chromosome <- chromosome[valid]
position <- position[valid]
key <- paste(chromosome, position, ref, alt, sep = ":")
keep_unique <- !duplicated(key)
raw <- raw[keep_unique, , drop = FALSE]
ref <- ref[keep_unique]
alt <- alt[keep_unique]
beta <- beta[keep_unique]
se <- se[keep_unique]
eaf <- eaf[keep_unique]
p_value <- p_value[keep_unique]
chromosome <- chromosome[keep_unique]
position <- position[keep_unique]

base_rows <- seq_len(min(rows, nrow(raw)))
low_rows <- which(is.finite(p_value) & p_value <= max_threshold)
non_low_base <- setdiff(base_rows, low_rows)
needed_base <- rows - length(low_rows)
if (needed_base < 0L || length(non_low_base) < needed_base) {
  stop("the fixture could not retain the requested rows and low-p rows")
}
keep <- sort(unique(c(low_rows, head(non_low_base, needed_base))))
if (length(keep) != rows) stop("fixture row count is not exactly requested")

sumstats <- data.frame(
  chromosome = chromosome[keep],
  base_pair_location = position[keep],
  reference_allele = ref[keep],
  alternate_allele = alt[keep],
  effect_allele = alt[keep],
  other_allele = ref[keep],
  beta = beta[keep],
  standard_error = se[keep],
  effect_allele_frequency = eaf[keep],
  z = beta[keep] / se[keep],
  stringsAsFactors = FALSE
)
fixture_p <- p_value[keep]
synthetic_hit_rows <- unique(round(seq.int(1L, rows,
                                             length.out = synthetic_hit_count)))

directory_bytes <- function(path) {
  files <- list.files(path, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  if (!length(files)) return(0)
  as.numeric(sum(file.info(files)$size, na.rm = TRUE))
}

write_native <- function(data, path) {
  started <- proc.time()[["elapsed"]]
  store <- compress_sumstats(
    data, path, input_build = "GRCh38", store_build = "GRCh38",
    backend = "pcodec", qc = "none", threads = threads,
    overwrite = TRUE, pvalue_flag = FALSE
  )
  list(store = store, elapsed_seconds = unname(proc.time()[["elapsed"]] - started))
}

read_flag <- function(path, info) {
  blob <- readBin(path, raw(), n = file.info(path)$size)
  pieces <- lapply(info$blocks, function(block) {
    first <- as.integer(block$offset) + 1L
    last <- first + as.integer(block$length) - 1L
    pcodec_native_decompress(blob[first:last], as.integer(block$values), "u8")
  })
  unlist(pieces, use.names = FALSE)
}

run_read <- function(store, columns) {
  timer <- system.time({
    value <- read_sumstats(store, columns = columns)
  })[["elapsed"]]
  as.numeric(timer)
}

run_pvalue_filter <- function(store, threshold) {
  timer <- system.time({
    value <- read_sumstats(store, columns = "p_value")
    selected <- is.finite(value$p_value) & value$p_value <= threshold
    selected_rows <- sum(selected)
  })[["elapsed"]]
  list(seconds = as.numeric(timer), rows = as.integer(selected_rows))
}

one_run <- function(run_id) {
  run_dir <- file.path(run_root, sprintf("run-%02d", run_id))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  baseline_path <- file.path(run_dir, "baseline.cpr")
  baseline_result <- write_native(sumstats, baseline_path)
  baseline_store <- baseline_result$store
  baseline_bytes <- directory_bytes(baseline_path)
  baseline_read <- median(vapply(seq_len(repeats), function(i) {
    run_read(baseline_store, c("z", "standard_error"))
  }, numeric(1)))
  pvalue_filter_runs <- lapply(seq_len(repeats), function(i) {
    run_pvalue_filter(baseline_store, filter_threshold)
  })
  pvalue_filter_seconds <- median(vapply(
    pvalue_filter_runs, function(value) value$seconds, numeric(1)
  ))
  pvalue_filter_rows <- unique(vapply(
    pvalue_filter_runs, function(value) value$rows, integer(1)
  ))
  if (length(pvalue_filter_rows) != 1L) stop("p-value filter row count varied")

  scenarios <- c(
    lapply(thresholds, function(threshold) list(
      scenario = paste0("p-", format(threshold, scientific = TRUE)),
      selection = "source_p_value",
      threshold = threshold,
      hit_rows = which(is.finite(fixture_p) & fixture_p <= threshold)
    )),
    list(list(
      scenario = paste0("designated-", synthetic_hit_count),
      selection = "synthetic_cardinality",
      threshold = NA_real_,
      hit_rows = synthetic_hit_rows
    ))
  )

  rows_out <- lapply(scenarios, function(scenario) {
    threshold <- scenario$threshold
    hit_rows <- scenario$hit_rows
    hit_count <- length(hit_rows)
    hit_data <- sumstats[hit_rows, , drop = FALSE]
    label <- scenario$scenario
    duplicate_path <- file.path(run_dir, paste0(label, ".cpr"))
    flag_dir <- file.path(run_dir, paste0(label, "-flag"))
    dir.create(flag_dir, recursive = TRUE, showWarnings = FALSE)

    duplicate_result <- write_native(hit_data, duplicate_path)
    duplicate_bytes <- directory_bytes(duplicate_path)
    duplicate_read <- median(vapply(seq_len(repeats), function(i) {
      run_read(duplicate_result$store, c("z", "standard_error"))
    }, numeric(1)))

    flag <- integer(rows)
    flag[hit_rows] <- 1L
    flag_path <- file.path(flag_dir, "significant.pco")
    flag_started <- proc.time()[["elapsed"]]
    flag_info <- pcodec_native_append_stream(
      flag, flag_path, dtype = "u8", workers = threads
    )
    flag_write <- unname(proc.time()[["elapsed"]] - flag_started)
    flag_meta <- list(
      codec = "Pcodec uint8",
      semantic = "aligned_binary_significance_flag",
      threshold = threshold,
      selection = scenario$selection,
      rows = rows,
      hit_rows = hit_count,
      stream = flag_info
    )
    write_json(flag_meta, file.path(flag_dir, "flag.json"),
               auto_unbox = TRUE, pretty = TRUE, digits = 17)
    flag_bytes <- directory_bytes(flag_dir)
    flag_read <- median(vapply(seq_len(repeats), function(i) {
      started <- proc.time()[["elapsed"]]
      decoded <- read_flag(flag_path, flag_info)
      elapsed <- unname(proc.time()[["elapsed"]] - started)
      if (length(decoded) != rows || sum(decoded) != hit_count) {
        stop("flag round-trip failed")
      }
      elapsed
    }, numeric(1)))

    data.table(
      run = run_id,
      scenario = scenario$scenario,
      selection = scenario$selection,
      threshold = threshold,
      hit_rows = hit_count,
      hit_fraction = hit_count / rows,
      pvalue_filter_threshold = filter_threshold,
      pvalue_filter_rows = pvalue_filter_rows,
      pvalue_filter_seconds = pvalue_filter_seconds,
      baseline_bytes = baseline_bytes,
      duplicate_sidecar_bytes = duplicate_bytes,
      flag_stream_bytes = flag_bytes,
      duplicate_additional_bytes = duplicate_bytes,
      flag_additional_bytes = flag_bytes,
      duplicate_total_bytes = baseline_bytes + duplicate_bytes,
      flag_total_bytes = baseline_bytes + flag_bytes,
      baseline_write_seconds = baseline_result$elapsed_seconds,
      duplicate_write_seconds = duplicate_result$elapsed_seconds,
      flag_write_seconds = flag_write,
      duplicate_total_write_seconds = baseline_result$elapsed_seconds +
        duplicate_result$elapsed_seconds,
      flag_total_write_seconds = baseline_result$elapsed_seconds + flag_write,
      baseline_read_seconds = baseline_read,
      duplicate_read_seconds = duplicate_read,
      flag_decode_seconds = flag_read,
      duplicate_additional_fraction = if (baseline_bytes) duplicate_bytes / baseline_bytes else NA_real_,
      flag_additional_fraction = if (baseline_bytes) flag_bytes / baseline_bytes else NA_real_
    )
  })
  rbindlist(rows_out, fill = TRUE)
}

records <- rbindlist(lapply(seq_len(repeats), one_run), fill = TRUE)
summary <- records[, lapply(.SD, median), by = .(scenario, selection, threshold),
                   .SDcols = setdiff(names(records),
                                     c("run", "scenario", "selection", "threshold"))]
setorder(summary, selection, threshold, na.last = TRUE)

run_csv_path <- file.path(result_root, "pvalue-sidechannel-macmini-runs.csv")
csv_path <- file.path(result_root, "pvalue-sidechannel-macmini-summary.csv")
md_path <- file.path(result_root, "pvalue-sidechannel-macmini-summary.md")
fwrite(records, run_csv_path)
fwrite(summary, csv_path)

lines <- c(
  "# FinnGen p-value side-channel benchmark",
  "",
  sprintf("Mac mini; %s valid prepared GRCh38 FinnGen rows.",
          format(rows, big.mark = ",")),
  "The fixture includes all source valid-SNP rows through p <= 1e-5.",
  sprintf("The p-value filter path reads reconstructed p and filters at p <= %.0e.",
          filter_threshold),
  "",
  "| Scenario | Selected rows | Fraction | Baseline B | Duplicate extra B | Flag extra B | P-filter read+filter s | Duplicate extra write s | Flag extra write s |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(summary)), function(i) {
    row <- summary[i]
    sprintf("| %s | %s | %.4f%% | %s | %s | %s | %.4f | %.4f | %.4f |",
            row$scenario,
            format(round(row$hit_rows), big.mark = ","),
            100 * row$hit_fraction,
            format(round(row$baseline_bytes), big.mark = ","),
            format(round(row$duplicate_additional_bytes), big.mark = ","),
            format(round(row$flag_additional_bytes), big.mark = ","),
            row$pvalue_filter_seconds,
            row$duplicate_write_seconds,
            row$flag_write_seconds)
  }, character(1)),
  "",
  "Duplicate extra bytes are a second current native Pcodec store containing the selected rows.",
  "Flag extra bytes are a Pcodec uint8 stream plus its small metadata/index, modelled as an additional main-store stream.",
  sprintf("The designated-%d row scenario is a cardinality stress test, not source-level significance.",
          synthetic_hit_count),
  "Write times are incremental side-channel encoding times; all rows are already prepared and the benchmark uses qc=none.",
  ""
)
writeLines(lines, md_path, useBytes = TRUE)

source_commit <- tryCatch(
  system2("git", c("-C", source_root, "rev-parse", "HEAD"),
          stdout = TRUE, stderr = FALSE)[[1L]],
  error = function(e) NA_character_
)
provenance <- list(
  benchmark = "pvalue_sidechannel_macmini",
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  source = normalizePath(input_path),
  source_rows_in_fixture = rows,
  fixture_source_threshold = max_threshold,
  pvalue_filter_threshold = filter_threshold,
  synthetic_hit_count = synthetic_hit_count,
  fixture_hits = lapply(thresholds, function(threshold) {
    list(threshold = threshold,
         rows = sum(is.finite(fixture_p) & fixture_p <= threshold))
  }),
  thresholds = thresholds,
  repeats = repeats,
  threads = threads,
  source_commit = source_commit,
  r = R.version.string,
  platform = Sys.info()[["sysname"]],
  architecture = R.version$arch
)
write_json(provenance,
           file.path(result_root, "pvalue-sidechannel-macmini-provenance.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 17)

print(summary)
cat("Wrote", csv_path, "and", md_path, "\n")
