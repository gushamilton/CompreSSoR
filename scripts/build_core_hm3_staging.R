#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: build_core_hm3_staging.R CORE.tsv.gz HM3.tsv.gz OUTPUT.tsv.gz", call. = FALSE)
}
core_path <- args[[1L]]
hm3_path <- args[[2L]]
output_path <- args[[3L]]
if (!file.exists(core_path) || !file.exists(hm3_path)) {
  stop("staged core and HM3 inputs must exist", call. = FALSE)
}
core <- fread(core_path, select = "variant_id", showProgress = FALSE)
hm3 <- fread(hm3_path, select = "variant_id", showProgress = FALSE)
core[, hm3 := as.integer(variant_id %chin% hm3[["variant_id"]])]
temporary <- sub("[.]gz$", "", output_path, ignore.case = TRUE)
fwrite(core, temporary, sep = "\t", quote = FALSE)
status <- system2("gzip", c("-n", "-f", temporary))
if (!identical(status, 0L) || !file.exists(output_path)) {
  stop("could not write combined staged panel", call. = FALSE)
}
cat("rows=", nrow(core), " hm3_rows=", sum(core[["hm3"]]),
    " output=", output_path, "\n", sep = "")
