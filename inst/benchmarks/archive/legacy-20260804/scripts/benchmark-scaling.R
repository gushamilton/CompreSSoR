#!/usr/bin/env Rscript

# Real-GWAS size scaling benchmark: TSV.gz, Parquet q9, native Pcodec.
# Each size is written once and read three times. Formats are removed after
# each size so BP scratch is not filled with redundant benchmark copies.
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(CompreSSoR)
})

source_path <- Sys.getenv("COMPRESSOR_SCALING_SOURCE", unset = "")
result_path <- Sys.getenv("COMPRESSOR_SCALING_RESULT", unset = "")
scratch_root <- Sys.getenv("COMPRESSOR_SCALING_SCRATCH", unset = "")
runs <- as.integer(Sys.getenv("COMPRESSOR_SCALING_RUNS", unset = "3"))
threads <- as.integer(Sys.getenv("COMPRESSOR_SCALING_PCODEC_THREADS", unset = "4"))
if (!file.exists(source_path) || !nzchar(result_path) || !nzchar(scratch_root)) {
  stop("source, result and scratch paths are required")
}
if (!is.finite(runs) || runs < 3L || !is.finite(threads) || threads < 1L) {
  stop("invalid runs or threads")
}

data <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
              showProgress = FALSE)
setnames(data, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p"))
data[, chrom := toupper(sub("^CHR", "", as.character(chrom)))]
data[, z := beta / se]

grch38_lengths <- c(
  `1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
  `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
  `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
  `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
  `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
  `21` = 46709983, `22` = 50818468)
grch38_offsets <- c(0, cumsum(as.numeric(grch38_lengths)))[seq_along(grch38_lengths)]
names(grch38_offsets) <- names(grch38_lengths)
data[, global_pos := unname(grch38_offsets[chrom]) + as.numeric(pos) - 1]
base_codes <- c(A = 0L, C = 1L, G = 2L, T = 3L)
data[, substitution := as.integer(base_codes[ref]) * 4L + as.integer(base_codes[alt])]

z_min <- -8
z_max <- 8
se_min <- min(log2(data$se), na.rm = TRUE)
se_max <- max(log2(data$se), na.rm = TRUE)
quantise_linear <- function(x, lo, hi, bits) {
  as.integer(round(pmin(1, pmax(0, (x - lo) / (hi - lo))) * (2^bits - 1)))
}
quantise_eaf <- function(x) {
  as.integer(round(asin(sqrt(pmin(1, pmax(0, x)))) / (pi / 2) * 255))
}
decode_linear <- function(x, lo, hi, bits) {
  lo + as.numeric(x) / (2^bits - 1) * (hi - lo)
}
decode_eaf <- function(x) sin(as.numeric(x) / 255 * pi / 2)^2
file_bytes <- function(path) {
  if (dir.exists(path)) {
    f <- list.files(path, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
    f <- f[file.exists(f) & !file.info(f)$isdir]
    return(sum(file.info(f)$size, na.rm = TRUE))
  }
  as.numeric(file.info(path)$size)
}
checksum <- function(x) sum(vapply(x, function(v) sum(as.numeric(v), na.rm = TRUE), 0.0))
read_tsv <- function(path) {
  x <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(path))), showProgress = FALSE)
  checksum(x[, .(pos, beta, se, eaf)])
}
read_parquet_q9 <- function(path) {
  x <- as.data.table(arrow::read_parquet(file.path(path, "values.parquet")))
  checksum(data.table(pos = x$pos, sub = x$sub,
                      z = decode_linear(x$z_code, z_min, z_max, 9L),
                      se = 2^decode_linear(x$se_code, se_min, se_max, 8L),
                      eaf = decode_eaf(x$eaf_code)))
}
read_pcodec <- function(path) {
  x <- read_sumstats(path, columns = c("global_position", "substitution", "z",
                                       "standard_error", "effect_allele_frequency"),
                     threads = threads)
  checksum(x)
}
write_parquet_q9 <- function(path, x) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  table <- data.table(
    pos = as.numeric(x$global_pos), sub = as.integer(x$substitution),
    z_code = quantise_linear(x$z, z_min, z_max, 9L),
    eaf_code = quantise_eaf(x$eaf),
    se_code = quantise_linear(log2(x$se), se_min, se_max, 8L))
  arrow::write_parquet(table, file.path(path, "values.parquet"),
                       compression = "zstd", compression_level = 9L,
                       write_statistics = TRUE, chunk_size = 65536L)
}
write_pcodec <- function(path, x) {
  compress_sumstats(x[, .(chromosome = chrom, base_pair_location = pos,
                          effect_allele = alt, other_allele = ref, beta,
                          standard_error = se,
                          effect_allele_frequency = eaf, z)],
                    path, reference = NULL, mode = "convert", backend = "pcodec",
                    profile = "standard", assume_grch38_ref_alt = TRUE,
                    overwrite = TRUE)
}

size_text <- Sys.getenv("COMPRESSOR_SCALING_SIZES", unset = "")
if (nzchar(size_text)) {
  # Keep requested sizes numeric: e.g. 10,000 * 1,000,000 overflows
  # 32-bit R integers even though the source row count is perfectly valid.
  sizes <- as.numeric(strsplit(size_text, ",", fixed = TRUE)[[1L]]) * 1000000
} else {
  sizes <- seq.int(1000000L, floor(nrow(data) / 1000000L) * 1000000L,
                   by = 1000000L)
  if (tail(sizes, 1L) < nrow(data)) sizes <- c(sizes, nrow(data))
}
sizes <- sort(unique(sizes[sizes <= nrow(data) & sizes > 0L]))
if (!length(sizes)) stop("no requested sizes are available")
dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
dir.create(scratch_root, recursive = TRUE, showWarnings = FALSE)
records <- list()
for (n in sizes) {
  x <- data[seq_len(n)]
  size_root <- file.path(scratch_root, paste0(n, "-rows"))
  dir.create(size_root, recursive = TRUE, showWarnings = FALSE)
  candidates <- list(
    list(id = "tsv_gz", label = "TSV.gz", path = file.path(size_root, "data.tsv.gz"),
         writer = function(path) fwrite(x[, .(chrom, pos, ref, alt, beta, se, eaf, p)], path, compress = "gzip"),
         reader = read_tsv),
    list(id = "parquet_q9", label = "Parquet q9", path = file.path(size_root, "parquet-q9"),
         writer = function(path) write_parquet_q9(path, x), reader = read_parquet_q9),
    list(id = "pcodec", label = "CompreSSoR Pcodec", path = file.path(size_root, "pcodec.cpr"),
         writer = function(path) write_pcodec(path, x), reader = read_pcodec)
  )
  for (candidate in candidates) {
    gc(FALSE)
    write_start <- proc.time()[["elapsed"]]
    tryCatch(candidate$writer(candidate$path), error = function(e) stop(candidate$label, " write failed: ", conditionMessage(e)))
    write_seconds <- proc.time()[["elapsed"]] - write_start
    storage_bytes <- file_bytes(candidate$path)
    values <- numeric(runs)
    for (run in seq_len(runs)) {
      gc(FALSE)
      value <- NULL
      elapsed <- system.time(value <- candidate$reader(candidate$path))[["elapsed"]]
      values[[run]] <- as.numeric(elapsed)
      if (!is.finite(value)) stop(candidate$label, " checksum failed")
    }
    records[[length(records) + 1L]] <- data.table(
      rows = n, format_id = candidate$id, format = candidate$label,
      runs = runs, storage_bytes = storage_bytes,
      write_seconds = write_seconds,
      read_median_seconds = median(values),
      read_min_seconds = min(values), read_max_seconds = max(values),
      checksum = value)
    unlink(candidate$path, recursive = TRUE, force = TRUE)
  }
  unlink(size_root, recursive = TRUE, force = TRUE)
  cat("completed", n, "rows\n")
}
result <- rbindlist(records)
tsv_sizes <- result[format_id == "tsv_gz", .(rows, source_tsvgz_bytes = storage_bytes)]
result <- merge(result, tsv_sizes, by = "rows", all.x = TRUE, sort = FALSE)
result[, relative_to_tsvgz := source_tsvgz_bytes / storage_bytes]
fwrite(result, result_path)
cat("wrote", result_path, "\n")
print(result[order(rows, format_id)])
