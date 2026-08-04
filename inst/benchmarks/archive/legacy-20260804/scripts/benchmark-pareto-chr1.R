#!/usr/bin/env Rscript

# Same-data Pareto benchmark for the current package and the main competing
# representations. It intentionally keeps one comparison contract: every
# plotted format stores its variant identity in the file. Reference-anchored
# numeric projections are not part of this final benchmark.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(arrow)
})

source_path <- Sys.getenv("COMPRESSOR_PARETO_SOURCE", unset = "")
output_root <- Sys.getenv("COMPRESSOR_PARETO_ROOT", unset = "pareto-chr1")
runs <- as.integer(Sys.getenv("COMPRESSOR_PARETO_RUNS", unset = "5"))
if (!nzchar(source_path) || !file.exists(source_path)) {
  stop("set COMPRESSOR_PARETO_SOURCE to the FinnGen chr1 TSV.gz")
}
if (!is.finite(runs) || runs < 5L) stop("runs must be at least five")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
format_root <- file.path(output_root, "formats")
dir.create(format_root, recursive = TRUE, showWarnings = FALSE)

file_bytes <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE,
                        all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE & file.exists(files)]
    return(sum(file.info(files)$size, na.rm = TRUE))
  }
  unname(file.info(path)$size)
}

read_input <- function() {
  x <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
             showProgress = FALSE)
  setnames(x, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p"))
  x[, z := beta / se]
  x
}

checksum <- function(x) {
  if (is.null(x)) stop("benchmark reader returned NULL")
  numeric_columns <- vapply(x, is.numeric, logical(1))
  values <- if (is.data.frame(x)) {
    if (is.data.table(x)) x[, numeric_columns, with = FALSE] else
      x[, numeric_columns, drop = FALSE]
  } else {
    Filter(is.numeric, x)
  }
  sum(vapply(values, function(value) sum(as.numeric(value), na.rm = TRUE), numeric(1)))
}

timed <- function(expr) {
  gc(FALSE)
  value <- NULL
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(seconds = as.numeric(elapsed), checksum = checksum(value))
}

source <- read_input()
n <- nrow(source)
source_bytes <- unname(file.info(source_path)$size)
self <- source[, .(chrom, pos, ref, alt, z, se, eaf)]

write_timed <- function(label, fun) {
  started <- proc.time()[["elapsed"]]
  fun()
  elapsed <- proc.time()[["elapsed"]] - started
  data.table(format = label, write_seconds = as.numeric(elapsed))
}

write_records <- list()
add_write <- function(label, fun) {
  write_records[[length(write_records) + 1L]] <<- write_timed(label, fun)
}

paths <- list(
  raw_tsv = file.path(format_root, "source.tsv"),
  parquet_self = file.path(format_root, "parquet-self"),
  vcf = file.path(format_root, "sumstats.vcf.gz"),
  pcodec = file.path(format_root, "compressor.cpr")
)

add_write("TSV uncompressed", function() {
  fwrite(source[, .(chrom, pos, ref, alt, beta, se, eaf, p)], paths$raw_tsv)
})
add_write("Parquet self-contained exact", function() {
  dir.create(paths$parquet_self, showWarnings = FALSE)
  arrow::write_parquet(self, file.path(paths$parquet_self, "values.parquet"),
                       compression = "zstd", compression_level = 9L,
                       write_statistics = TRUE, chunk_size = 65536L)
})
add_write("VCF bgzip + Tabix", function() {
  plain <- file.path(format_root, "sumstats.vcf")
  header <- c(
    "##fileformat=VCFv4.3", "##contig=<ID=1>",
    "##INFO=<ID=BETA,Number=1,Type=Float,Description=Beta>",
    "##INFO=<ID=SE,Number=1,Type=Float,Description=SE>",
    "##INFO=<ID=EAF,Number=1,Type=Float,Description=EAF>",
    "##INFO=<ID=P,Number=1,Type=Float,Description=P>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )
  writeLines(header, plain)
  vcf <- source[, .(chrom, pos, id = ".", ref, alt, qual = ".", filter = "PASS",
                    info = sprintf("BETA=%s;SE=%s;EAF=%s;P=%s", beta, se, eaf, p))]
  fwrite(vcf, plain, sep = "\t", col.names = FALSE, append = TRUE)
  status <- system2("bcftools", c("sort", "-m", "4G", "-Oz", "-o",
                                   paths$vcf, plain))
  if (!identical(status, 0L)) stop("bcftools sort failed")
  status <- system2("bcftools", c("index", "-t", "-f", paths$vcf))
  if (!identical(status, 0L)) stop("bcftools index failed")
  unlink(plain)
})
add_write("CompreSSoR native Pcodec", function() {
  compress_sumstats(source[, .(chromosome = chrom, base_pair_location = pos,
                               effect_allele = alt, other_allele = ref,
                               beta, standard_error = se,
                               effect_allele_frequency = eaf, z)],
                    paths$pcodec, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
})

# Fail early if the actual package store is not a faithful keyed round trip.
# The numerical profile is bounded-lossy by design, so validate against the
# documented tolerances rather than demanding bitwise equality.
decoded <- read_sumstats(
  paths$pcodec,
  columns = c("chromosome", "base_pair_location", "effect_allele",
              "other_allele", "z", "standard_error",
              "effect_allele_frequency", "beta")
)
source_key <- paste(source$chrom, source$pos, source$ref, source$alt, sep = ":")
decoded_key <- paste(decoded$chromosome, decoded$base_pair_location,
                     decoded$other_allele, decoded$effect_allele, sep = ":")
decoded_index <- match(source_key, decoded_key)
if (anyNA(decoded_index) || anyDuplicated(decoded_key)) {
  stop("CompreSSoR native round trip changed or duplicated variant identity")
}
validation <- data.table(
  metric = c("rows", "max_abs_z_error", "max_abs_se_error", "max_abs_eaf_error",
             "max_se_relative_error", "max_abs_beta_error",
             "max_beta_error_over_quantisation_bound"),
  value = c(
    nrow(decoded),
    max(abs(decoded$z[decoded_index] - source$z), na.rm = TRUE),
    max(abs(decoded$standard_error[decoded_index] - source$se), na.rm = TRUE),
    max(abs(decoded$effect_allele_frequency[decoded_index] - source$eaf), na.rm = TRUE),
    max(abs(decoded$standard_error[decoded_index] - source$se) /
        pmax(abs(source$se), 1e-8), na.rm = TRUE),
    max(abs(decoded$beta[decoded_index] - source$beta), na.rm = TRUE),
    {
      z_tol <- 7 / (2 * (2^9 - 2))
      se_rel_tol <- 2^(4 / 254) - 1
      beta_bound <- 1.02 * (abs(source$se) * z_tol +
                            abs(source$z) * abs(source$se) * se_rel_tol) + 1e-6
      max(abs(decoded$beta[decoded_index] - source$beta) / beta_bound, na.rm = TRUE)
    }
  )
)
validation_value <- function(label) validation$value[match(label, validation$metric)]
if (validation_value("rows") != nrow(source) ||
    validation_value("max_abs_z_error") > 0.02 ||
    validation_value("max_abs_eaf_error") > 0.004 ||
    validation_value("max_se_relative_error") > (2^(4 / 254) - 1) * 1.01 ||
    validation_value("max_beta_error_over_quantisation_bound") > 1) {
  stop("CompreSSoR native round-trip tolerance failed")
}

read_vcf <- function(path) {
  x <- fread(cmd = paste("bcftools query -f '%INFO/BETA\\t%INFO/SE\\t%INFO/EAF\\n'",
                         shQuote(path)), header = FALSE, showProgress = FALSE)
  setnames(x, c("beta", "se", "eaf"))
  x[, z := beta / se]
  x[, .(z, se, eaf)]
}

readers <- list(
  "TSV gzip" = function() {
    x <- read_input()
    x[, .(z, se, eaf)]
  },
  "TSV uncompressed" = function() {
    x <- fread(paths$raw_tsv, showProgress = FALSE)
    x[, .(z = beta / se, se, eaf)]
  },
  "Parquet self-contained exact" = function() {
    x <- as.data.table(arrow::read_parquet(file.path(paths$parquet_self, "values.parquet")))
    x[, .(z, se, eaf)]
  },
  "VCF bgzip + Tabix" = function() read_vcf(paths$vcf),
  "CompreSSoR Pcodec self-contained" = function() read_sumstats(
    paths$pcodec, columns = c("chromosome", "base_pair_location", "effect_allele",
                              "other_allele", "z", "standard_error",
                              "effect_allele_frequency"))
)

key_formats <- c(
  "TSV gzip", "TSV uncompressed", "Parquet self-contained exact",
  "VCF bgzip + Tabix", "CompreSSoR Pcodec self-contained"
)
readers <- readers[key_formats]

contracts <- c(
  "TSV gzip" = "self-contained text",
  "TSV uncompressed" = "self-contained text",
  "Parquet self-contained exact" = "self-contained",
  "VCF bgzip + Tabix" = "self-contained",
  "CompreSSoR Pcodec self-contained" = "self-contained"
)

storage <- data.table(
  format = names(readers),
  contract = unname(contracts[names(readers)]),
  storage_bytes = c(
    source_bytes,
    file_bytes(paths$raw_tsv),
    file_bytes(paths$parquet_self),
    file_bytes(paths$vcf) + file_bytes(paste0(paths$vcf, ".tbi")),
    file_bytes(paths$pcodec)
  ), source_bytes = source_bytes
)
storage <- storage[format %in% key_formats]
storage[, compression_ratio := source_bytes / storage_bytes]

access_records <- list()
for (label in names(readers)) {
  for (run in seq_len(runs)) {
    measured <- timed(readers[[label]]())
    access_records[[length(access_records) + 1L]] <- data.table(
      format = label, run = run, seconds = measured$seconds,
      checksum = measured$checksum, rows = n
    )
  }
}
access_records <- rbindlist(access_records)
access_summary <- access_records[, .(
  runs = .N, median_seconds = median(seconds), min_seconds = min(seconds),
  max_seconds = max(seconds), checksum_median = median(checksum), rows = max(rows)
), by = format]

summary <- merge(storage, access_summary, by = "format", sort = FALSE)
summary[, benchmark_id := "compressor_chr1_pareto_same_data_v1"]
setcolorder(summary, c("benchmark_id", "format", "contract", "storage_bytes",
                       "source_bytes", "compression_ratio", "median_seconds",
                       "min_seconds", "max_seconds", "runs", "rows",
                       "checksum_median"))

fwrite(summary, file.path(output_root, "pareto-summary.csv"))
fwrite(access_records, file.path(output_root, "pareto-access-runs.csv"))
write_summary <- rbindlist(write_records)[format %in% c(key_formats,
                                                          "CompreSSoR native Pcodec")]
fwrite(write_summary, file.path(output_root, "pareto-write.csv"))
fwrite(validation, file.path(output_root, "pareto-validation.csv"))
writeLines(c(
  "# Same-data chr1 Pareto benchmark",
  "",
  paste0("Source: ", normalizePath(source_path)),
  paste0("Rows: ", format(n, big.mark = ",")),
  paste0("Source bytes: ", source_bytes),
  paste0("Runs per access format: ", runs),
  "Only formats that store a variant identity key are included. All access
checks materialise numeric values and consume a checksum."
), file.path(output_root, "README.txt"))
print(summary)
