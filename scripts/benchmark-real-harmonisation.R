#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(CompreSSoR))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: benchmark-real-harmonisation.R INPUT REFERENCE OUTPUT")
input <- args[[1L]]
reference <- list(
  id = "bp-real-benchmark",
  build = "GRCh38",
  cache_dir = "/tmp/CompreSSoR-speed-reference-cache-final",
  variants = list(path = args[[2L]])
)
results <- list()
for (threads in c(1L, 4L)) {
  gc()
  tm <- system.time({
    got <- harmonise_sumstats(input, reference = reference, mode = "qc", chrom_threads = threads)
  })
  results[[length(results) + 1L]] <- data.frame(
    scenario = if (threads == 1L) "qc_serial" else "qc_chrom_parallel",
    chrom_threads = threads,
    input_rows = nrow(got),
    aligned_rows = sum(got$harmonisation_status == "aligned"),
    unresolved_rows = sum(got$harmonisation_status != "aligned"),
    elapsed_seconds = unname(tm[["elapsed"]]),
    stringsAsFactors = FALSE
  )
  rm(got)
}
results <- do.call(rbind, results)
dir.create(dirname(args[[3L]]), recursive = TRUE, showWarnings = FALSE)
write.csv(results, args[[3L]], row.names = FALSE)
print(results, row.names = FALSE)
