#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: build_bundled_core_hm3_panel.R COMBINED.tsv.gz STAGING-MANIFEST.json OUTPUT.cpr",
       call. = FALSE)
}
combined_path <- args[[1L]]
staging_manifest_path <- args[[2L]]
output_path <- args[[3L]]
combined <- fread(combined_path, showProgress = FALSE)
if (!all(c("variant_id", "hm3") %in% names(combined))) {
  stop("combined panel must contain variant_id and hm3", call. = FALSE)
}
staging_manifest <- fromJSON(staging_manifest_path, simplifyVector = FALSE)
core_source <- staging_manifest$panels$core
hm3_source <- staging_manifest$panels$hm3
source_provenance <- list(
  core = list(
    source = core_source$source,
    source_sha256 = core_source$source_sha256,
    source_rows = as.integer(core_source$source_rows),
    valid_core_rows = as.integer(core_source$prepared_rows),
    excluded_invalid_rows = as.integer(core_source$excluded_invalid_rows),
    prepared_sha256 = core_source$prepared_sha256
  ),
  hm3 = list(
    source = hm3_source$source,
    source_sha256 = hm3_source$source_sha256,
    source_rows = as.integer(hm3_source$source_rows),
    valid_hm3_rows = as.integer(hm3_source$prepared_rows),
    prepared_sha256 = hm3_source$prepared_sha256
  ),
  hm3_rows_in_core_universe = as.integer(sum(combined$hm3))
)
panel <- write_variant_panel(
  combined, output_path, build = "GRCh38",
  source_provenance = source_provenance, overwrite = TRUE
)
files <- file.info(list.files(output_path, full.names = TRUE))
cat("rows=", panel$manifest$n_rows,
    " hm3_rows=", panel$manifest$hm3_rows,
    " bytes=", sum(files$size),
    " output=", output_path, "\n", sep = "")
