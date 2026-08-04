#!/usr/bin/env Rscript

# Compare the native 0.4 implementation before/after the semantic SE8 change
# on one deterministic fixture. The fixture is generated inside each clean R
# process so the two package builds see the same rows and source bytes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("usage: benchmark-native-se-codec.R LIBRARY OUTPUT_CSV [STORE_DIR]",
       call. = FALSE)
}
library(CompreSSoR, lib.loc = args[[1L]])
library(data.table)

n <- 1000000L
set.seed(20260804)
position <- seq.int(100001L, length.out = n)
reference <- rep(c("A", "C", "G", "T"), length.out = n)
alternate <- rep(c("C", "G", "T", "A"), length.out = n)
eaf <- pmin(0.999, pmax(0.001, rbeta(n, 1.4, 3.2)))
standard_error <- exp(rnorm(n, log(0.05), 0.45))
z <- rnorm(n)
data <- data.frame(
  chromosome = rep("1", n), base_pair_location = position,
  effect_allele = alternate, other_allele = reference,
  beta = z * standard_error, standard_error = standard_error,
  effect_allele_frequency = eaf, z = z,
  p_value = 2 * pnorm(-abs(z)),
  variant_id = paste0("1:", position, ":", reference, ":", alternate),
  rsid = paste0("rs", position), stringsAsFactors = FALSE
)

root <- if (length(args) >= 3L) args[[3L]] else tempfile("native-se-codec-")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
source_plain <- file.path(root, "source.tsv")
source_gz <- paste0(source_plain, ".gz")
store_path <- file.path(root, "store.cpr")
fwrite(data, source_plain, sep = "\t", quote = FALSE, na = "NA")
system2("gzip", c("-f", source_plain))

block_rows <- as.integer(Sys.getenv("COMPRESSOR_NATIVE_BLOCK_ROWS", "32768"))
write_elapsed <- system.time({
  compress_sumstats(data, store_path, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE,
                    block_rows = block_rows)
})[["elapsed"]]
store_files <- list.files(store_path, recursive = TRUE, full.names = TRUE)
store_bytes <- sum(file.info(store_files)$size)

timed <- function(fun) {
  unname(as.numeric(system.time(invisible(fun()))[["elapsed"]]))
}
source_times <- vapply(seq_len(5L), function(i) timed(function() {
  data.table::fread(cmd = paste("gzip -dc", shQuote(source_gz)),
                    data.table = FALSE, showProgress = FALSE,
                    check.names = FALSE)
}), numeric(1))
native_times <- vapply(seq_len(5L), function(i) timed(function() {
  read_sumstats(store_path)
}), numeric(1))
region_times <- vapply(seq_len(5L), function(i) timed(function() {
  read_sumstats(store_path, region = "1:100001-110000",
                columns = c("z", "standard_error", "effect_allele_frequency"))
}), numeric(1))
keys <- compressor_variant_key(data$chromosome, data$base_pair_location,
                               data$other_allele, data$effect_allele)
sparse_keys <- keys[seq.int(1001L, n, by = floor(n / 25L))][seq_len(25L)]
sparse_times <- vapply(seq_len(5L), function(i) timed(function() {
  read_sumstats(store_path, variants = sparse_keys,
                columns = c("z", "standard_error", "effect_allele_frequency"))
}), numeric(1))

manifest <- jsonlite::read_json(file.path(store_path, "manifest.json"), simplifyVector = FALSE)
result <- data.frame(
  package_version = as.character(packageVersion("CompreSSoR")),
  host = Sys.info()[["nodename"]], rows = n, runs = 5L,
  source_gzip_bytes = file.info(source_gz)$size,
  native_store_bytes = store_bytes,
  compression_ratio_vs_source_gzip = file.info(source_gz)$size / store_bytes,
  source_median_seconds = median(source_times),
  native_median_seconds = median(native_times),
  native_region_median_seconds = median(region_times),
  native_sparse25_median_seconds = median(sparse_times),
  write_seconds = write_elapsed,
  native_format = manifest$format_version,
  block_rows = block_rows,
  codec = manifest$codec$name,
  exception_rows = manifest$semantic_codec$exception_rows,
  se_bits = manifest$semantic_codec$se_bits,
  se_residual_range = paste(unlist(manifest$semantic_codec$se_residual_range), collapse = ":"),
  stringsAsFactors = FALSE
)
write.csv(result, args[[2L]], row.names = FALSE, quote = TRUE)
write.csv(data.frame(access = c(rep("source_gzip_full_read", 5L),
                                rep("native_full_read", 5L)),
                     run = rep(seq_len(5L), 2L),
                     seconds = c(source_times, native_times)),
          sub("\\.csv$", "-runs.csv", args[[2L]]), row.names = FALSE)
print(result)
print(data.frame(source_times = source_times, native_times = native_times,
                 region_times = region_times, sparse25_times = sparse_times))
