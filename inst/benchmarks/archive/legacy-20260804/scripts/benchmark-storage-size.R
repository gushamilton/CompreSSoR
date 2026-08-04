#!/usr/bin/env Rscript

# Measure on-disk storage for the real release-gate GWAS. Run on the Mac mini;
# temporary exact/q8 stores live on the external SSD and are removed at exit.

suppressPackageStartupMessages({
  library(CompreSSoR)
})

root <- Sys.getenv("COMPRESSOR_BENCH_ROOT", unset = "")
project <- normalizePath(getwd(), mustWork = TRUE)
raw_dir <- file.path(root, "raw")
store_dir <- file.path(root, "stores")
report_dir <- file.path(root, "reports")
source_path <- Sys.getenv("COMPRESSOR_SOURCE", unset = "")
reference_path <- Sys.getenv("COMPRESSOR_REFERENCE", unset = "")
if (!nzchar(root) || !nzchar(source_path) || !nzchar(reference_path)) {
  stop("set COMPRESSOR_BENCH_ROOT, COMPRESSOR_SOURCE, and COMPRESSOR_REFERENCE")
}
standard_path <- file.path(store_dir, "finngen-shared.cpr")
temporary_dir <- file.path(store_dir, "storage-size-temporary")
inline_path <- file.path(temporary_dir, "convert-inline.cpr")
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
reference <- list(
  id = "bp_shared_hm3", build = "GRCh38",
  source = "MR-atlas shared GRCh38 dictionary", variants = reference_path,
  cache_dir = file.path(root, "cache", "reference")
)

# Rebuild the measured stores with the current writer, so the benchmark cannot
# silently compare a previous layout with a new implementation. The standard
# store is reference-anchored; inline storage is retained only as a comparison.
invisible(compress_sumstats(source_path, standard_path, reference = reference, mode = "qc",
                            profile = "standard", keep_extras = FALSE, overwrite = TRUE))
invisible(compress_sumstats(source_path, inline_path, reference = NULL, mode = "convert",
                            profile = "standard", keep_extras = FALSE, overwrite = TRUE))
standard <- open_compressor(standard_path)
n_rows <- standard$manifest$n_rows

# Build the two comparison formats only for this measurement. The source
# standard store remains available for normal tests; temporary outputs are
# removed by on.exit().
invisible(compress_sumstats(
  source_path, exact_path, reference = reference, mode = "qc",
  profile = "exact", overwrite = TRUE
))
build_cache(standard, output = q8_path, overwrite = TRUE, block_rows = 65536L)

canonical <- CompreSSoR:::reference_table(reference)
canonical_path <- attr(canonical, "reference_metadata")$normalized_cache_path
canonical_bytes <- file_bytes(canonical_path)
canonical_rows <- nrow(canonical)

standard_payload_bytes <- file_bytes(file.path(standard_path, "values.parquet")) +
  file_bytes(file.path(standard_path, "exceptions.parquet"))
standard_index_bytes <- file_bytes(file.path(standard_path, "variants.parquet"))
standard_unmatched_bytes <- file_bytes(file.path(standard_path, "unmatched.parquet"))

measure <- data.frame(
  benchmark_id = "compressor_storage_size_finngen_16111549",
  representation = c(
    "Source gzip", "Raw TSV", "Shared BP canonical spine",
    "q9 numeric payload", "q9 shared reference index",
    "q9 unresolved-variant table", "q9 shared-reference store",
    "q9 inline store (comparison)", "Exact shared-reference store",
    "q8 cache only (optional serving layer)",
    "q9 shared-reference + q8 cache"
  ),
  backend = c("gzip", "text", "Parquet", "Parquet", "Parquet", "Parquet", "Parquet", "Parquet", "Parquet", "framed-gzip", "combined"),
  profile = c("source", "raw", "shared-reference", "q9-payload", "q9-index", "q9-unmatched", "standard", "inline", "exact", "q8", "standard+q8"),
  rows = c(as.integer(n_rows), as.integer(n_rows), as.integer(canonical_rows),
           rep(as.integer(n_rows), 8L)),
  bytes = c(
    source_bytes, raw_bytes, canonical_bytes, standard_payload_bytes,
    standard_index_bytes, standard_unmatched_bytes, file_bytes(standard_path),
    file_bytes(inline_path),
    file_bytes(exact_path), file_bytes(q8_path), file_bytes(standard_path) + file_bytes(q8_path)
  ),
  stringsAsFactors = FALSE
)
measure$bytes_per_row <- measure$bytes / measure$rows
measure$size_vs_source_gzip <- measure$bytes / source_bytes
measure$size_vs_raw_tsv <- measure$bytes / raw_bytes
measure$size_runs <- 5L
measure$note <- c(
  "Baseline distributed input",
  "Uncompressed text stream; not retained on disk",
  "One shared BP canonical GRCh38 spine; not duplicated per study",
  "q9 values, p-values and exceptions only",
  "Per-study zero-based indices into the shared canonical spine",
  "Local identity rows retained for variants absent/ambiguous in the spine",
  "Default durable reference-anchored output",
  "No reference; full inline identity comparison",
  "Lossless reference-anchored comparison",
  "Optional q8 values cache; not standalone without the q9 spine",
  "Durable q9 store plus optional q8 cache and shared spine"
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
  "The BP canonical spine is shared and hash-pinned; it is not copied into the per-study store.",
  "The q8 cache is an optional serving layer and is not a standalone replacement for the q9 payload/index.", "",
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

study_counts <- c(1L, 5L, 100L)
amortization <- data.frame(
  benchmark_id = "compressor_storage_amortization_finngen_16111549",
  studies = study_counts,
  shared_canonical_bytes_per_study = canonical_bytes / study_counts,
  q9_shared_store_bytes_per_study = file_bytes(standard_path) + canonical_bytes / study_counts,
  q9_shared_plus_q8_bytes_per_study = file_bytes(standard_path) + file_bytes(q8_path) + canonical_bytes / study_counts,
  source_gzip_bytes_per_study = source_bytes,
  stringsAsFactors = FALSE
)
write.csv(amortization, file.path(project, "outputs", "storage-amortization.csv"), row.names = FALSE)
file.copy(file.path(project, "outputs", "storage-amortization.csv"),
          file.path(project, "inst", "benchmarks", "storage-amortization.csv"), overwrite = TRUE)
writeLines(c(
  "# Shared-spine storage amortization", "",
  "The BP canonical spine is paid once and amortized across studies.", "",
  "| Studies | Shared canonical bytes/study | q9 shared store bytes/study | q9 + q8 bytes/study |",
  "|---:|---:|---:|---:|",
  vapply(seq_len(nrow(amortization)), function(i) paste0(
    "| ", amortization$studies[i], " | ", fmt_bytes(amortization$shared_canonical_bytes_per_study[i]),
    " | ", fmt_bytes(amortization$q9_shared_store_bytes_per_study[i]),
    " | ", fmt_bytes(amortization$q9_shared_plus_q8_bytes_per_study[i]), " |"
  ), character(1))
), file.path(project, "outputs", "storage-amortization.md"))
saveRDS(list(results = measure, source = source_path, session = sessionInfo()),
        file.path(report_dir, paste0(result_stem, ".rds")))
unlink(temporary_dir, recursive = TRUE, force = TRUE)
message("Wrote ", result_path)
