#!/usr/bin/env Rscript

# Prepare one larger, immutable real-GWAS input for the scaling benchmark.
# The benchmark intentionally stops at the number of valid SNP rows present
# in this public FinnGen trait rather than padding or duplicating observations.
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

source_path <- Sys.getenv("COMPRESSOR_SCALING_RAW_SOURCE", unset = "")
output_path <- Sys.getenv("COMPRESSOR_SCALING_SOURCE", unset = "")
metadata_path <- Sys.getenv("COMPRESSOR_SCALING_METADATA", unset = "")
target_rows <- as.integer(Sys.getenv("COMPRESSOR_SCALING_ROWS", unset = "15000000"))
if (!nzchar(source_path) || !file.exists(source_path)) stop("raw source is missing")
if (!nzchar(output_path) || !nzchar(metadata_path)) stop("output paths are required")
if (!is.finite(target_rows) || target_rows < 1L) stop("invalid target row count")

raw <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
             select = c("#chrom", "pos", "ref", "alt", "pval", "beta", "sebeta", "maf"),
             showProgress = FALSE)
setnames(raw, c("chrom", "pos", "ref", "alt", "p", "beta", "se", "eaf"))
raw[, chrom := toupper(sub("^CHR", "", as.character(chrom)))]
raw[, ref := toupper(as.character(ref))]
raw[, alt := toupper(as.character(alt))]
base <- c("A", "C", "G", "T")
keep <- raw$chrom %in% as.character(1:22) &
  nchar(raw$ref) == 1L & nchar(raw$alt) == 1L &
  raw$ref %in% base & raw$alt %in% base & raw$ref != raw$alt &
  is.finite(raw$pos) & is.finite(raw$beta) & is.finite(raw$se) & raw$se > 0 &
  is.finite(raw$eaf) & raw$eaf >= 0 & raw$eaf <= 1 & is.finite(raw$p)
raw <- raw[keep]
raw[, key := paste(chrom, pos, ref, alt, sep = ":")]
raw <- raw[!duplicated(key)]
setorder(raw, chrom, pos, ref, alt)
if (nrow(raw) < target_rows) {
  stop("only ", nrow(raw), " valid SNPs; requested ", target_rows)
}
out <- raw[seq_len(target_rows), .(chrom, pos, ref, alt, beta, se, eaf, p)]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(out, output_path, sep = "\t", compress = "gzip")
metadata <- list(
  source = normalizePath(source_path),
  source_sha256 = digest(source_path, algo = "sha256", file = TRUE),
  output = normalizePath(output_path),
  output_sha256 = digest(output_path, algo = "sha256", file = TRUE),
  rows = nrow(out),
  columns = names(out),
  source_mapping = list(eaf = "FinnGen maf, used as the benchmark EAF field"),
  filters = c("autosomes 1-22", "biallelic A/C/G/T SNPs", "finite beta/se/EAF/P", "se > 0", "unique CHR:POS:REF:ALT"),
  ordering = "lexicographic chromosome then position/ref/alt; prefix after filtering"
)
write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE)
cat("prepared", nrow(out), "rows; bytes", file.info(output_path)$size,
    "; sha256", metadata$output_sha256, "\n")
