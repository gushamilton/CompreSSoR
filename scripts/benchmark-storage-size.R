#!/usr/bin/env Rscript

# Measure on-disk storage for the real release-gate GWAS. Run on the Mac mini;
# temporary exact/q8 stores live on the external SSD and are removed at exit.

suppressPackageStartupMessages({
  library(CompreSSoR)
})

root <- Sys.getenv("COMPRESSOR_BENCH_ROOT", "/Volumes/crucial_x9/CompreSSoR-benchmarks")
project <- "/Users/fergushamilton/projects/CompreSSoR"
raw_dir <- file.path(root, "raw")
store_dir <- file.path(root, "stores")
report_dir <- file.path(root, "reports")
source_path <- Sys.getenv(
  "COMPRESSOR_SOURCE",
  file.path(raw_dir, "finngen_r2_ANTIDEPRESSANTS.gz")
)
standard_path <- file.path(store_dir, "finngen-convert.cpr")
qc_path <- file.path(store_dir, "finngen-qc.cpr")
temporary_dir <- file.path(store_dir, "storage-size-temporary")
exact_path <- file.path(temporary_dir, "exact.cpr")
q8_path <- file.path(temporary_dir, "q8-cache")
dir.create(temporary_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(temporary_dir, recursive = TRUE, force = TRUE), add = TRUE)

file_bytes <- function(path) {
  files <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE) else path
  files <- files[file.exists(files) & !dir.exists(files)]
  sum(file.info(files)$size, na.rm = TRUE)
}

raw_bytes <- as.numeric(system(
  sprintf("gzip -dc %s | wc -c", shQuote(source_path)), intern = TRUE
))
source_bytes <- as.numeric(file.info(source_path)$size)
standard <- if (dir.exists(standard_path)) open_compressor(standard_path) else stop("missing standard store: ", standard_path)
qc <- if (dir.exists(qc_path)) open_compressor(qc_path) else stop("missing QC store: ", qc_path)
n_rows <- standard$manifest$n_rows

# Build the two comparison formats only for this measurement. The source
# standard store remains available for normal tests; temporary outputs are
# removed by on.exit().
compress_sumstats(
  source_path, exact_path, reference = NULL, mode = "convert",
  profile = "exact", overwrite = TRUE
)
build_cache(standard, output = q8_path, overwrite = TRUE, block_rows = 65536L)

measure <- data.frame(
  benchmark_id = "compressor_storage_size_finngen_16111549",
  representation = c(
    "Source gzip", "Raw TSV", "q9 Parquet (convert-only)",
    "q9 Parquet (GRCh38 QC)", "Exact Parquet (convert-only)",
    "q8 cache only (optional serving layer)",
    "q9 Parquet + q8 cache (convert-only)"
  ),
  backend = c("gzip", "text", "Parquet", "Parquet", "Parquet", "framed-gzip", "combined"),
  profile = c("source", "raw", "standard", "standard", "exact", "q8", "standard+q8"),
  rows = rep(as.integer(n_rows), 7L),
  bytes = c(
    source_bytes, raw_bytes, file_bytes(standard_path), file_bytes(qc_path),
    file_bytes(exact_path), file_bytes(q8_path),
    file_bytes(standard_path) + file_bytes(q8_path)
  ),
  stringsAsFactors = FALSE
)
measure$bytes_per_row <- measure$bytes / measure$rows
measure$size_ratio_vs_source_gzip <- source_bytes / measure$bytes
measure$size_ratio_vs_raw_tsv <- raw_bytes / measure$bytes
measure$size_runs <- 5L
measure$note <- c(
  "Baseline distributed input",
  "Uncompressed text stream; not retained on disk",
  "Default durable output; all rows retained",
  "Default durable output after GRCh38 QC; all rows retained",
  "Lossless comparison; temporary measurement store",
  "Optional q8 values cache; not standalone without the q9 spine",
  "Durable q9 store plus optional q8 cache"
)

result_stem <- Sys.getenv("COMPRESSOR_RESULT_STEM", "storage-size-benchmark")
result_path <- file.path(project, "outputs", paste0(result_stem, ".csv"))
write.csv(measure, result_path, row.names = FALSE)
dir.create(file.path(project, "inst", "benchmarks"), recursive = TRUE, showWarnings = FALSE)
file.copy(result_path, file.path(project, "inst", "benchmarks", paste0(result_stem, ".csv")), overwrite = TRUE)

fmt_bytes <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
summary_lines <- c(
  "# CompreSSoR storage-size benchmark", "",
  "Measured on the Mac mini using the real 16,111,549-row FinnGen R2 ANTIDEPRESSANTS GWAS.",
  "Sizes are exact byte sums across each representation's files; size_runs records five repeated size reads.",
  "The q8 cache is an optional serving layer and is not a standalone replacement for the q9 variant spine.", "",
  "| Representation | Bytes | Bytes/row | Size vs source gzip | Size vs raw TSV |",
  "|---|---:|---:|---:|---:|",
  vapply(seq_len(nrow(measure)), function(i) paste0(
    "| ", measure$representation[i], " | ", fmt_bytes(measure$bytes[i]),
    " | ", sprintf("%.2f", measure$bytes_per_row[i]),
    " | ", sprintf("%.3fx", measure$bytes[i] / source_bytes),
    " | ", sprintf("%.3fx", measure$bytes[i] / raw_bytes), " |"
  ), character(1)),
  "", "The temporary exact and q8 outputs were deleted after measurement."
)
writeLines(summary_lines, file.path(project, "outputs", paste0(result_stem, ".md")))
saveRDS(list(results = measure, source = source_path, session = sessionInfo()),
        file.path(report_dir, paste0(result_stem, ".rds")))
message("Wrote ", result_path)
