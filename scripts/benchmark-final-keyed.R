#!/usr/bin/env Rscript

# Final self-contained format screen for real FinnGen chr1 data.
# Every candidate contains the same exact identity: position plus directed
# REF->ALT substitution. No shared spine or external reference is counted.
# Each Slurm array task is one independent benchmark repetition.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(jsonlite)
  library(CompreSSoR)
})

source_path <- Sys.getenv("COMPRESSOR_FINAL_SOURCE", unset = "")
result_root <- Sys.getenv("COMPRESSOR_FINAL_ROOT", unset = "final-keyed")
run_id <- as.integer(Sys.getenv("COMPRESSOR_FINAL_RUN", unset = "1"))
if (!nzchar(source_path) || !file.exists(source_path)) {
  stop("set COMPRESSOR_FINAL_SOURCE to the FinnGen chr1 TSV.gz")
}
if (!is.finite(run_id) || run_id < 1L) stop("invalid COMPRESSOR_FINAL_RUN")

run_root <- file.path(result_root, paste0("run-", run_id))
format_root <- file.path(run_root, "formats")
dir.create(format_root, recursive = TRUE, showWarnings = FALSE)

file_bytes <- function(path) {
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE,
                        all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
    files <- files[file.exists(files) & !file.info(files)$isdir]
    return(sum(file.info(files)$size, na.rm = TRUE))
  }
  if (!file.exists(path)) return(NA_real_)
  as.numeric(file.info(path)$size)
}

source_data <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
                     showProgress = FALSE)
setnames(source_data, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf", "p"))
source_data[, z := beta / se]
n <- nrow(source_data)
source_bytes <- as.numeric(file.info(source_path)$size)

base_codes <- c(A = 0L, C = 1L, G = 2L, T = 3L)
substitution <- function(ref, alt) {
  as.integer(base_codes[ref]) * 4L + as.integer(base_codes[alt])
}
source_sub <- substitution(source_data$ref, source_data$alt)

read_keyed_table <- function(x, z = NULL) {
  if (is.null(z)) z <- x$beta / x$se
  data.table(pos = as.integer(x$pos), sub = as.integer(x$sub),
             z = as.numeric(z), se = as.numeric(x$se),
             eaf = as.numeric(x$eaf))
}

checksum <- function(x) {
  sum(as.numeric(x$pos), as.numeric(x$sub), as.numeric(x$z),
      as.numeric(x$se), as.numeric(x$eaf), na.rm = TRUE)
}

quant_spec <- list(
  z_min = -8, z_max = 8,
  eaf_min = 0, eaf_max = 1,
  se_min = min(log2(source_data$se), na.rm = TRUE),
  se_max = max(log2(source_data$se), na.rm = TRUE)
)

quantise_linear <- function(x, lo, hi, bits) {
  levels <- 2^bits - 1
  as.integer(round(pmin(1, pmax(0, (x - lo) / (hi - lo))) * levels))
}
decode_linear <- function(code, lo, hi, bits) {
  levels <- 2^bits - 1
  lo + as.numeric(code) / levels * (hi - lo)
}
quantise_eaf <- function(x, bits) {
  as.integer(round(asin(sqrt(pmin(1, pmax(0, x)))) / (pi / 2) * (2^bits - 1)))
}
decode_eaf <- function(code, bits) {
  sin(as.numeric(code) / (2^bits - 1) * pi / 2)^2
}

header_raw <- function(mode, n, bits_z = 0L, bits_eaf = 0L, bits_se = 0L) {
  description <- paste0(
    "CompreSSoR-FINAL-KEYED-v1;mode=", mode, ";n=", n,
    ";z_bits=", bits_z, ";eaf_bits=", bits_eaf, ";se_bits=", bits_se,
    ";z_min=", quant_spec$z_min, ";z_max=", quant_spec$z_max,
    ";eaf_min=", quant_spec$eaf_min, ";eaf_max=", quant_spec$eaf_max,
    ";se_min=", quant_spec$se_min, ";se_max=", quant_spec$se_max)
  raw <- charToRaw(substr(description, 1L, 255L))
  c(raw, raw(256L - length(raw)))
}

write_stream <- function(con, x, type) {
  if (identical(type, "u8")) writeBin(as.raw(x), con, size = 1L)
  else if (identical(type, "u16")) writeBin(as.integer(x), con, size = 2L, endian = "little")
  else if (identical(type, "u32")) writeBin(as.integer(x), con, size = 4L, endian = "little")
  else if (identical(type, "f32")) writeBin(as.numeric(x), con, size = 4L, endian = "little")
  else if (identical(type, "f64")) writeBin(as.numeric(x), con, size = 8L, endian = "little")
  else stop("unknown stream type: ", type)
}

read_stream <- function(con, n, type) {
  if (identical(type, "u8")) as.integer(readBin(con, raw(), n = n, size = 1L))
  else if (identical(type, "u16")) readBin(con, integer(), n = n, size = 2L, endian = "little")
  else if (identical(type, "u32")) readBin(con, integer(), n = n, size = 4L, endian = "little")
  else if (identical(type, "f32")) readBin(con, numeric(), n = n, size = 4L, endian = "little")
  else if (identical(type, "f64")) readBin(con, numeric(), n = n, size = 8L, endian = "little")
  else stop("unknown stream type: ", type)
}

binary_vectors <- function(mode) {
  if (mode == "f64") {
    list(bits_z = 0L, bits_eaf = 0L, bits_se = 0L,
         z = list(source_data$z, "f64"), eaf = list(source_data$eaf, "f64"),
         se = list(source_data$se, "f64"))
  } else if (mode == "f32") {
    list(bits_z = 0L, bits_eaf = 0L, bits_se = 0L,
         z = list(source_data$z, "f32"), eaf = list(source_data$eaf, "f32"),
         se = list(source_data$se, "f32"))
  } else if (mode %in% c("q8", "q9", "q12", "q16")) {
    bz <- if (mode == "q8") 8L else if (mode == "q9") 9L else if (mode == "q12") 12L else 9L
    bs <- if (mode == "q16") 16L else 8L
    z_code <- quantise_linear(source_data$z, quant_spec$z_min, quant_spec$z_max, bz)
    eaf_code <- quantise_eaf(source_data$eaf, 8L)
    se_code <- quantise_linear(log2(source_data$se), quant_spec$se_min,
                               quant_spec$se_max, bs)
    list(bits_z = bz, bits_eaf = 8L, bits_se = bs,
         z = list(z_code, if (bz <= 8L) "u8" else "u16"),
         eaf = list(eaf_code, "u8"),
         se = list(se_code, if (bs <= 8L) "u8" else "u16"))
  } else stop("unknown binary mode: ", mode)
}

write_binary <- function(path, mode, framed = FALSE, zstd = FALSE) {
  vectors <- binary_vectors(mode)
  raw_path <- if (zstd) paste0(path, ".raw") else path
  con <- file(raw_path, open = "wb")
  writeBin(header_raw(mode, n, vectors$bits_z, vectors$bits_eaf, vectors$bits_se), con)
  write_stream(con, source_data$pos, "u32")
  write_stream(con, source_sub, "u8")
  if (!framed) {
    write_stream(con, vectors$z[[1L]], vectors$z[[2L]])
    write_stream(con, vectors$eaf[[1L]], vectors$eaf[[2L]])
    write_stream(con, vectors$se[[1L]], vectors$se[[2L]])
  } else {
    frame_rows <- 65536L
    for (start in seq.int(1L, n, by = frame_rows)) {
      inside <- start:min(n, start + frame_rows - 1L)
      write_stream(con, length(inside), "u32")
      write_stream(con, vectors$z[[1L]][inside], vectors$z[[2L]])
      write_stream(con, vectors$eaf[[1L]][inside], vectors$eaf[[2L]])
      write_stream(con, vectors$se[[1L]][inside], vectors$se[[2L]])
    }
  }
  close(con)
  if (zstd) {
    status <- system2("zstd", c("-q", "-f", "-19", raw_path, "-o", path))
    if (!identical(status, 0L)) stop("zstd failed")
    unlink(raw_path)
  }
}

read_binary_connection <- function(con, mode, framed = FALSE) {
  readBin(con, raw(), n = 256L)
  vectors <- binary_vectors(mode)
  pos <- read_stream(con, n, "u32")
  sub <- read_stream(con, n, "u8")
  if (!framed) {
    z_code <- read_stream(con, n, vectors$z[[2L]])
    eaf_code <- read_stream(con, n, vectors$eaf[[2L]])
    se_code <- read_stream(con, n, vectors$se[[2L]])
  } else {
    z_code <- integer(n); eaf_code <- integer(n); se_code <- integer(n)
    offset <- 1L
    while (offset <= n) {
      frame_n <- read_stream(con, 1L, "u32")
      inside <- offset:(offset + frame_n - 1L)
      z_code[inside] <- read_stream(con, frame_n, vectors$z[[2L]])
      eaf_code[inside] <- read_stream(con, frame_n, vectors$eaf[[2L]])
      se_code[inside] <- read_stream(con, frame_n, vectors$se[[2L]])
      offset <- offset + frame_n
    }
  }
  if (mode %in% c("f64", "f32")) {
    z <- z_code; eaf <- eaf_code; se <- se_code
  } else {
    z <- decode_linear(z_code, quant_spec$z_min, quant_spec$z_max, vectors$bits_z)
    eaf <- decode_eaf(eaf_code, vectors$bits_eaf)
    se <- 2^decode_linear(se_code, quant_spec$se_min, quant_spec$se_max, vectors$bits_se)
  }
  data.table(pos = as.integer(pos), sub = as.integer(sub), z = as.numeric(z),
             se = as.numeric(se), eaf = as.numeric(eaf))
}

read_binary <- function(path, mode, framed = FALSE, zstd = FALSE) {
  if (zstd) con <- pipe(paste("zstd -q -d -c", shQuote(path)), open = "rb")
  else con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  read_binary_connection(con, mode, framed)
}

write_tsv <- function(path) {
  fwrite(source_data[, .(chrom, pos, ref, alt, beta, se, eaf, p)], path)
}
read_tsv <- function(path) {
  x <- fread(path, showProgress = FALSE)
  x[, sub := substitution(ref, alt)]
  read_keyed_table(x, x$beta / x$se)
}
read_tsv_gz <- function(path) {
  x <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(path))), showProgress = FALSE)
  x[, sub := substitution(ref, alt)]
  read_keyed_table(x, x$beta / x$se)
}

write_vcf <- function(path) {
  plain <- paste0(path, ".plain")
  header <- c("##fileformat=VCFv4.3", "##contig=<ID=1>",
              "##INFO=<ID=BETA,Number=1,Type=Float,Description=Beta>",
              "##INFO=<ID=SE,Number=1,Type=Float,Description=SE>",
              "##INFO=<ID=EAF,Number=1,Type=Float,Description=EAF>",
              "##INFO=<ID=P,Number=1,Type=Float,Description=P>",
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  writeLines(header, plain)
  x <- source_data[, .(chrom, pos, id = paste(chrom, pos, ref, alt, sep = ":"),
                       ref, alt, qual = ".", filter = "PASS",
                       info = sprintf("BETA=%s;SE=%s;EAF=%s;P=%s", beta, se, eaf, p))]
  fwrite(x, plain, sep = "\t", col.names = FALSE, append = TRUE)
  status <- system2("bcftools", c("sort", "-m", "4G", "-Oz", "-o", path, plain))
  if (!identical(status, 0L)) stop("bcftools sort failed")
  status <- system2("bcftools", c("index", "-t", "-f", path))
  if (!identical(status, 0L)) stop("bcftools index failed")
  unlink(plain)
}
read_vcf <- function(path) {
  query <- paste0("bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/BETA\\t%INFO/SE\\t%INFO/EAF\\n' ", shQuote(path))
  x <- fread(cmd = query, header = FALSE, showProgress = FALSE)
  setnames(x, c("chrom", "pos", "ref", "alt", "beta", "se", "eaf"))
  x[, sub := substitution(ref, alt)]
  read_keyed_table(x, x$beta / x$se)
}

self_table <- source_data[, .(pos = as.integer(pos), sub = source_sub,
                              z, se, eaf)]
q_table <- function(bits_z) {
  data.table(pos = as.integer(source_data$pos), sub = source_sub,
             z_code = quantise_linear(source_data$z, quant_spec$z_min,
                                      quant_spec$z_max, bits_z),
             eaf_code = quantise_eaf(source_data$eaf, 8L),
             se_code = quantise_linear(log2(source_data$se), quant_spec$se_min,
                                      quant_spec$se_max, 8L))
}
write_parquet <- function(path, table, metadata = NULL,
                          compression = "zstd", compression_level = 9L) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  args <- list(x = table, sink = file.path(path, "values.parquet"),
               compression = compression, write_statistics = TRUE,
               chunk_size = 65536L)
  if (compression %in% c("zstd", "gzip")) args$compression_level <- compression_level
  do.call(arrow::write_parquet, args)
  if (!is.null(metadata)) writeLines(toJSON(metadata, auto_unbox = TRUE),
                                     file.path(path, "metadata.json"))
}
read_parquet_exact <- function(path) as.data.table(
  arrow::read_parquet(file.path(path, "values.parquet")))
read_parquet_q <- function(path, bits_z, bits_se = 8L) {
  x <- as.data.table(arrow::read_parquet(file.path(path, "values.parquet")))
  data.table(pos = x$pos, sub = x$sub,
             z = decode_linear(x$z_code, quant_spec$z_min, quant_spec$z_max, bits_z),
             se = 2^decode_linear(x$se_code, quant_spec$se_min, quant_spec$se_max, bits_se),
             eaf = decode_eaf(x$eaf_code, 8L))
}

parquet_f32_table <- function() {
  arrow::Table$create(
    pos = as.integer(source_data$pos), sub = as.integer(source_sub),
    z = arrow::Array$create(as.numeric(source_data$z), type = arrow::float32()),
    se = arrow::Array$create(as.numeric(source_data$se), type = arrow::float32()),
    eaf = arrow::Array$create(as.numeric(source_data$eaf), type = arrow::float32())
  )
}

q_table_profile <- function(bits_z, bits_se = 8L) {
  data.table(pos = as.integer(source_data$pos), sub = source_sub,
             z_code = quantise_linear(source_data$z, quant_spec$z_min,
                                      quant_spec$z_max, bits_z),
             eaf_code = quantise_eaf(source_data$eaf, 8L),
             se_code = quantise_linear(log2(source_data$se), quant_spec$se_min,
                                       quant_spec$se_max, bits_se))
}

write_pcodec <- function(path) {
  compress_sumstats(source_data[, .(chromosome = chrom, base_pair_location = pos,
                                   effect_allele = alt, other_allele = ref,
                                   beta, standard_error = se,
                                   effect_allele_frequency = eaf, z)],
                    path, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
}
read_pcodec <- function(path) {
  x <- read_sumstats(path, columns = c("base_pair_location", "effect_allele",
                                       "other_allele", "z", "standard_error",
                                       "effect_allele_frequency"))
  data.table(pos = as.integer(x$base_pair_location),
             sub = substitution(x$other_allele, x$effect_allele),
             z = as.numeric(x$z), se = as.numeric(x$standard_error),
             eaf = as.numeric(x$effect_allele_frequency))
}

candidate <- function(id, label, family, path, writer, reader) {
  list(id = id, label = label, family = family, path = path,
       writer = writer, reader = reader)
}

candidates <- list(
  candidate("tsv_gzip", "TSV gzip", "Text", source_path, NULL, read_tsv_gz),
  candidate("tsv_plain", "TSV uncompressed", "Text",
            file.path(format_root, "source.tsv"), write_tsv, read_tsv),
  candidate("vcf_tabix", "VCF bgzip + Tabix", "VCF + Tabix",
            file.path(format_root, "source.vcf.gz"), write_vcf, read_vcf),
  candidate("parquet_exact", "Parquet exact self-contained", "Parquet",
            file.path(format_root, "parquet-exact"),
            function(path) write_parquet(path, self_table), read_parquet_exact),
  candidate("parquet_f32", "Parquet float32 self-contained", "Parquet",
            file.path(format_root, "parquet-f32"),
            function(path) write_parquet(path, parquet_f32_table()), read_parquet_exact),
  candidate("parquet_q8", "Parquet q8 self-contained", "Parquet",
            file.path(format_root, "parquet-q8"),
            function(path) write_parquet(path, q_table(8L), quant_spec),
            function(path) read_parquet_q(path, 8L)),
  candidate("parquet_q9", "Parquet q9 self-contained", "Parquet",
            file.path(format_root, "parquet-q9"),
            function(path) write_parquet(path, q_table(9L), quant_spec),
            function(path) read_parquet_q(path, 9L)),
  candidate("parquet_q12", "Parquet q12 self-contained", "Parquet",
            file.path(format_root, "parquet-q12"),
            function(path) write_parquet(path, q_table_profile(12L), quant_spec),
            function(path) read_parquet_q(path, 12L)),
  candidate("parquet_q16_u8", "Parquet q9/eaf8/se16 self-contained", "Parquet",
            file.path(format_root, "parquet-q16-u8"),
            function(path) write_parquet(path, q_table_profile(9L, 16L), quant_spec),
            function(path) read_parquet_q(path, 9L, 16L)),
  candidate("parquet_q8_snappy", "Parquet q8 + Snappy", "Parquet",
            file.path(format_root, "parquet-q8-snappy"),
            function(path) write_parquet(path, q_table(8L), compression = "snappy"),
            function(path) read_parquet_q(path, 8L)),
  candidate("parquet_q8_gzip", "Parquet q8 + Gzip", "Parquet",
            file.path(format_root, "parquet-q8-gzip"),
            function(path) write_parquet(path, q_table(8L), compression = "gzip", compression_level = 9L),
            function(path) read_parquet_q(path, 8L)),
  candidate("raw_f64", "Raw float64 keyed binary", "Binary",
            file.path(format_root, "raw-f64.bin"),
            function(path) write_binary(path, "f64"),
            function(path) read_binary(path, "f64")),
  candidate("raw_f64_zstd", "Raw float64 keyed binary + Zstd", "Binary",
            file.path(format_root, "raw-f64.zst"),
            function(path) write_binary(path, "f64", zstd = TRUE),
            function(path) read_binary(path, "f64", zstd = TRUE)),
  candidate("raw_f32", "Raw float32 keyed binary", "Binary",
            file.path(format_root, "raw-f32.bin"),
            function(path) write_binary(path, "f32"),
            function(path) read_binary(path, "f32")),
  candidate("raw_f32_zstd", "Raw float32 keyed binary + Zstd", "Binary",
            file.path(format_root, "raw-f32.zst"),
            function(path) write_binary(path, "f32", zstd = TRUE),
            function(path) read_binary(path, "f32", zstd = TRUE)),
  candidate("raw_f32_framed", "Raw float32 framed keyed binary", "Framed binary",
            file.path(format_root, "raw-f32-framed.bin"),
            function(path) write_binary(path, "f32", framed = TRUE),
            function(path) read_binary(path, "f32", framed = TRUE)),
  candidate("raw_f32_framed_zstd", "Raw float32 framed keyed binary + Zstd", "Framed binary",
            file.path(format_root, "raw-f32-framed.zst"),
            function(path) write_binary(path, "f32", framed = TRUE, zstd = TRUE),
            function(path) read_binary(path, "f32", framed = TRUE, zstd = TRUE)),
  candidate("q8_binary", "q8 keyed binary", "Binary",
            file.path(format_root, "q8.bin"),
            function(path) write_binary(path, "q8"),
            function(path) read_binary(path, "q8")),
  candidate("q8_binary_zstd", "q8 keyed binary + Zstd", "Binary",
            file.path(format_root, "q8.zst"),
            function(path) write_binary(path, "q8", zstd = TRUE),
            function(path) read_binary(path, "q8", zstd = TRUE)),
  candidate("q9_binary", "q9 keyed binary", "Binary",
            file.path(format_root, "q9.bin"),
            function(path) write_binary(path, "q9"),
            function(path) read_binary(path, "q9")),
  candidate("q9_binary_zstd", "q9 keyed binary + Zstd", "Binary",
            file.path(format_root, "q9.zst"),
            function(path) write_binary(path, "q9", zstd = TRUE),
            function(path) read_binary(path, "q9", zstd = TRUE)),
  candidate("q12_binary", "q12 keyed binary", "Binary",
            file.path(format_root, "q12.bin"),
            function(path) write_binary(path, "q12"),
            function(path) read_binary(path, "q12")),
  candidate("q12_binary_zstd", "q12 keyed binary + Zstd", "Binary",
            file.path(format_root, "q12.zst"),
            function(path) write_binary(path, "q12", zstd = TRUE),
            function(path) read_binary(path, "q12", zstd = TRUE)),
  candidate("q8_framed", "q8 framed keyed binary", "Framed binary",
            file.path(format_root, "q8-framed.bin"),
            function(path) write_binary(path, "q8", framed = TRUE),
            function(path) read_binary(path, "q8", framed = TRUE)),
  candidate("q8_framed_zstd", "q8 framed keyed binary + Zstd", "Framed binary",
            file.path(format_root, "q8-framed.zst"),
            function(path) write_binary(path, "q8", framed = TRUE, zstd = TRUE),
            function(path) read_binary(path, "q8", framed = TRUE, zstd = TRUE)),
  candidate("q9_framed", "q9 framed keyed binary", "Framed binary",
            file.path(format_root, "q9-framed.bin"),
            function(path) write_binary(path, "q9", framed = TRUE),
            function(path) read_binary(path, "q9", framed = TRUE)),
  candidate("q9_framed_zstd", "q9 framed keyed binary + Zstd", "Framed binary",
            file.path(format_root, "q9-framed.zst"),
            function(path) write_binary(path, "q9", framed = TRUE, zstd = TRUE),
            function(path) read_binary(path, "q9", framed = TRUE, zstd = TRUE)),
  candidate("q12_framed", "q12 framed keyed binary", "Framed binary",
            file.path(format_root, "q12-framed.bin"),
            function(path) write_binary(path, "q12", framed = TRUE),
            function(path) read_binary(path, "q12", framed = TRUE)),
  candidate("q12_framed_zstd", "q12 framed keyed binary + Zstd", "Framed binary",
            file.path(format_root, "q12-framed.zst"),
            function(path) write_binary(path, "q12", framed = TRUE, zstd = TRUE),
            function(path) read_binary(path, "q12", framed = TRUE, zstd = TRUE)),
  candidate("q16_u8", "q16/u8 keyed binary", "Binary",
            file.path(format_root, "q16-u8.bin"),
            function(path) write_binary(path, "q16"),
            function(path) read_binary(path, "q16")),
  candidate("q16_u8_framed", "q9/eaf8/se16 framed keyed binary", "Framed binary",
            file.path(format_root, "q16-u8-framed.bin"),
            function(path) write_binary(path, "q16", framed = TRUE),
            function(path) read_binary(path, "q16", framed = TRUE)),
  candidate("q16_u8_zstd", "q16/u8 keyed binary + Zstd", "Binary",
            file.path(format_root, "q16-u8.zst"),
            function(path) write_binary(path, "q16", zstd = TRUE),
            function(path) read_binary(path, "q16", zstd = TRUE)),
  candidate("q16_u8_framed_zstd", "q9/eaf8/se16 framed keyed binary + Zstd", "Framed binary",
            file.path(format_root, "q16-u8-framed.zst"),
            function(path) write_binary(path, "q16", framed = TRUE, zstd = TRUE),
            function(path) read_binary(path, "q16", framed = TRUE, zstd = TRUE)),
  candidate("pcodec_native", "CompreSSoR native Pcodec", "CompreSSoR",
            file.path(format_root, "compressor.cpr"), write_pcodec, read_pcodec)
)

# A final-candidate rerun can select exactly one already-defined format.  This
# keeps the measurement contract identical to the full matrix while avoiding
# another expensive rebuild of every historical comparator.
selected_id <- Sys.getenv("COMPRESSOR_FINAL_ONLY", unset = "")
if (nzchar(selected_id)) {
  candidates <- Filter(function(item) identical(item$id, selected_id), candidates)
  if (!length(candidates)) {
    stop("COMPRESSOR_FINAL_ONLY did not match a candidate: ", selected_id)
  }
}

validate <- function(value) {
  if (nrow(value) != n) stop("row count mismatch")
  source_key <- paste(source_data$pos, source_sub, sep = ":")
  value_key <- paste(value$pos, value$sub, sep = ":")
  if (anyNA(value$pos) || anyNA(value$sub) || anyDuplicated(source_key) ||
      anyDuplicated(value_key) || !setequal(value_key, source_key)) {
    stop("identity mismatch")
  }
  # VCF/Tabix and the native Pcodec store are allowed to reorder records by
  # their canonical key.  Compare values after aligning on the exact key.
  source_order <- order(source_data$pos, source_sub, method = "radix")
  value_order <- order(value$pos, value$sub, method = "radix")
  value <- value[value_order, , drop = FALSE]
  source_data_ordered <- source_data[source_order, , drop = FALSE]
  data.table(
    max_abs_z_error = max(abs(value$z - source_data_ordered$z), na.rm = TRUE),
    max_abs_se_error = max(abs(value$se - source_data_ordered$se), na.rm = TRUE),
    max_abs_eaf_error = max(abs(value$eaf - source_data_ordered$eaf), na.rm = TRUE),
    checksum = checksum(value)
  )
}

records <- list()
for (item in candidates) {
  write_seconds <- NA_real_; read_seconds <- NA_real_; storage_bytes <- NA_real_
  status <- "OK"; error <- ""
  validation <- data.table(max_abs_z_error = NA_real_, max_abs_se_error = NA_real_,
                           max_abs_eaf_error = NA_real_, checksum = NA_real_)
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    if (!is.null(item$writer)) item$writer(item$path)
    NULL
  }, error = function(e) e)
  write_seconds <- proc.time()[["elapsed"]] - started
  if (inherits(result, "error")) {
    status <- "WRITE_ERROR"; error <- conditionMessage(result)
  } else {
    storage_bytes <- file_bytes(item$path)
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(item$reader(item$path), error = function(e) e)
    read_seconds <- proc.time()[["elapsed"]] - started
    if (inherits(result, "error")) {
      status <- "READ_ERROR"; error <- conditionMessage(result)
    } else {
      check <- tryCatch(validate(result), error = function(e) e)
      if (inherits(check, "error")) {
        status <- "VALIDATION_ERROR"; error <- conditionMessage(check)
      } else validation <- check
    }
  }
  records[[length(records) + 1L]] <- data.table(
    run_id = run_id, format_id = item$id, format = item$label,
    family = item$family, identity_included = TRUE,
    source_bytes = source_bytes, storage_bytes = storage_bytes,
    compression_ratio = source_bytes / storage_bytes,
    write_seconds = as.numeric(write_seconds),
    full_read_seconds = as.numeric(read_seconds), status = status,
    error = error, rows = n,
    max_abs_z_error = validation$max_abs_z_error,
    max_abs_se_error = validation$max_abs_se_error,
    max_abs_eaf_error = validation$max_abs_eaf_error,
    checksum = validation$checksum
  )
}

records <- rbindlist(records, fill = TRUE)
fwrite(records, file.path(run_root, "final-keyed-run.csv"))
fwrite(data.table(
  format_id = vapply(candidates, `[[`, character(1), "id"),
  format = vapply(candidates, `[[`, character(1), "label"),
  family = vapply(candidates, `[[`, character(1), "family"),
  identity_included = TRUE,
  benchmark_contract = "exact position + directed REF->ALT key included"
), file.path(run_root, "candidate-registry.csv"))

if (!identical(Sys.getenv("COMPRESSOR_FINAL_KEEP_FORMATS", "0"), "1")) {
  unlink(format_root, recursive = TRUE)
}
print(records[, .(format, storage_bytes, compression_ratio, full_read_seconds, status)])
