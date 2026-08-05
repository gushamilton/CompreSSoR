#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
output_root <- if (length(args)) args[[1L]] else
  "/Volumes/crucial_x9/CompreSSoR-benchmarks/issue18-validation-20260805"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
store_root <- file.path(output_root, "stores")
dir.create(store_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(CompreSSoR))

make_fixture <- function(n, build = "GRCh38") {
  data.frame(
    chromosome = rep("1", n),
    base_pair_location = seq_len(n),
    effect_allele = rep("A", n),
    other_allele = rep("C", n),
    beta = sin(seq_len(n) / 1000),
    standard_error = 0.1 + (seq_len(n) %% 31) / 1000,
    effect_allele_frequency = 0.1 + (seq_len(n) %% 800) / 1000,
    stringsAsFactors = FALSE
  )
}

write_store <- function(data, path, threads) {
  started <- proc.time()[["elapsed"]]
  store <- compress_sumstats(
    data, path, reference = NULL, backend = "pcodec",
    assume_grch38_ref_alt = TRUE, overwrite = TRUE,
    input_build = "GRCh38", store_build = "GRCh38",
    threads = threads, row_policy = "error"
  )
  elapsed <- proc.time()[["elapsed"]] - started
  list(store = store, elapsed_seconds = unname(elapsed))
}

payload_files <- c("position.pco", "substitution.pco", "z.pco", "eaf.pco",
                   "se.pco", "exceptions.bin")
sha256_files <- function(path) {
  stats::setNames(
    vapply(payload_files, function(name) {
      digest::digest(file.path(path, name), algo = "sha256", file = TRUE)
    }, character(1)),
    payload_files
  )
}

n <- 200000L
fixture <- make_fixture(n)
serial_result <- write_store(fixture, file.path(store_root, "threads-1.cpr"), 1L)
parallel_result <- write_store(fixture, file.path(store_root, "threads-4.cpr"), 4L)

serial <- serial_result$store
parallel <- parallel_result$store
serial_validation <- validate_compressor(serial, full = TRUE)
parallel_validation <- validate_compressor(parallel, full = TRUE)
serial_hashes <- sha256_files(serial$path)
parallel_hashes <- sha256_files(parallel$path)

region <- read_sumstats(parallel,
                        region = "chr1:100-200",
                        columns = c("chromosome", "base_pair_location", "z"))
key_read <- read_sumstats(
  parallel, variants = c("1:100:C:A", "1:200:C:A"),
  columns = c("chromosome", "base_pair_location", "effect_allele", "other_allele")
)

grch37_data <- make_fixture(20000L, build = "GRCh37")
grch37_path <- file.path(store_root, "grch37.cpr")
grch37 <- compress_sumstats(
  grch37_data, grch37_path, reference = NULL, backend = "pcodec",
  assume_grch38_ref_alt = TRUE, overwrite = TRUE,
  input_build = "GRCh37", store_build = "GRCh37",
  threads = 2L, row_policy = "error"
)
grch37_validation <- validate_compressor(grch37, full = TRUE)
grch37_read <- read_sumstats(
  grch37, variants = "1:100:C:A",
  columns = c("chromosome", "base_pair_location", "effect_allele", "other_allele")
)

summary <- list(
  gate = "issue-18-macmini",
  timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  platform = list(r = R.version.string, arch = R.version$arch,
                  os = R.version$os),
  rows = n,
  serial_seconds = serial_result$elapsed_seconds,
  parallel_seconds = parallel_result$elapsed_seconds,
  serial_valid = isTRUE(serial_validation$valid),
  parallel_valid = isTRUE(parallel_validation$valid),
  deterministic_payload = identical(serial_hashes, parallel_hashes),
  serial_hashes = serial_hashes,
  parallel_hashes = parallel_hashes,
  region_rows = nrow(region),
  key_rows = nrow(key_read),
  writer = parallel$manifest$writer %||% NULL,
  grch37 = list(
    valid = isTRUE(grch37_validation$valid), rows = nrow(grch37_read),
    genome_build = grch37$manifest$genome_build,
    identity_table = as.character(
      grch37$manifest$identity$chromosome_table_id %||% NA_character_
    )
  )
)
summary$pass <- isTRUE(summary$serial_valid) && isTRUE(summary$parallel_valid) &&
  isTRUE(summary$deterministic_payload) && summary$region_rows == 101L &&
summary$key_rows == 2L && isTRUE(summary$grch37$valid) &&
  summary$grch37$rows == 1L && nzchar(summary$grch37$identity_table)
jsonlite::write_json(summary, file.path(output_root, "issue18-macmini-gate.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = 17)
if (!isTRUE(summary$pass)) stop("issue-18 Mac mini gate failed", call. = FALSE)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = 17), "\n")
