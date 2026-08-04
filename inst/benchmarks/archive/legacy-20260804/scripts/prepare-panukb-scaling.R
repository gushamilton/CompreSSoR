#!/usr/bin/env Rscript

# Prepare a large public Pan-UKB quantitative GWAS for size-scaling tests.
# The source is GRCh37. This benchmark retains the source coordinates because
# it measures representation scaling, not the CompreSSoR liftover pathway.
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

source_path <- Sys.getenv("COMPRESSOR_PANUKB_SOURCE", unset = "")
output_path <- Sys.getenv("COMPRESSOR_SCALING_SOURCE", unset = "")
metadata_path <- Sys.getenv("COMPRESSOR_SCALING_METADATA", unset = "")
target_rows <- as.integer(Sys.getenv("COMPRESSOR_SCALING_ROWS", unset = "21000000"))
if (!file.exists(source_path) || !nzchar(output_path) || !nzchar(metadata_path)) {
  stop("source, output and metadata paths are required")
}
if (!is.finite(target_rows) || target_rows < 1L) stop("invalid target row count")

raw <- fread(cmd = paste("gzip -dc", shQuote(normalizePath(source_path))),
             select = c("chr", "pos", "ref", "alt", "af_EUR", "beta_EUR",
                         "se_EUR", "neglog10_pval_EUR"),
             showProgress = FALSE)
setnames(raw, c("chrom", "pos", "ref", "alt", "eaf", "beta", "se", "p"))
raw[, chrom := toupper(sub("^CHR", "", as.character(chrom)))]
raw[, ref := toupper(as.character(ref))]
raw[, alt := toupper(as.character(alt))]
base <- c("A", "C", "G", "T")
grch38_lengths <- c(
  `1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555,
  `5` = 181538259, `6` = 170805979, `7` = 159345973, `8` = 145138636,
  `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
  `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345,
  `17` = 83257441, `18` = 80373285, `19` = 58617616, `20` = 64444167,
  `21` = 46709983, `22` = 50818468)
keep <- raw$chrom %in% as.character(1:22) &
  nchar(raw$ref) == 1L & nchar(raw$alt) == 1L &
  raw$ref %in% base & raw$alt %in% base & raw$ref != raw$alt &
  is.finite(raw$pos) & is.finite(raw$beta) & is.finite(raw$se) & raw$se > 0 &
  is.finite(raw$eaf) & raw$eaf >= 0 & raw$eaf <= 1 & is.finite(raw$p) &
  raw$pos <= unname(grch38_lengths[raw$chrom])
raw <- raw[keep]
raw[, key := paste(chrom, pos, ref, alt, sep = ":")]
raw <- raw[!duplicated(key)]
if (nrow(raw) < target_rows) {
  stop("only ", nrow(raw), " valid autosomal biallelic SNPs; requested ", target_rows)
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
  source_assembly = "GRCh37",
  source_mapping = list(
    eaf = "Pan-UKB af_EUR",
    beta = "Pan-UKB beta_EUR",
    se = "Pan-UKB se_EUR",
    p = "Pan-UKB neglog10_pval_EUR (retained as the benchmark p field)"
  ),
  filters = c("autosomes 1-22", "biallelic A/C/G/T SNPs", "finite beta/se/EAF/logP", "se > 0", "unique CHR:POS:REF:ALT", "positions within GRCh38 primary chromosome lengths for collision-free benchmark key"),
  ordering = "source order after filtering"
)
write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE)
cat("prepared", nrow(out), "rows; bytes", file.info(output_path)$size,
    "; sha256", metadata$output_sha256, "\n")
