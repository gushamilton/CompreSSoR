#!/usr/bin/env Rscript

## Compact input-projection evidence for issue #27.
##
## Usage:
##   /usr/bin/time -v Rscript scripts/benchmark-input-projection.R input.tsv.gz [summary.tsv]
##
## The input is read through the native-core projection path. Keep the source
## and any /usr/bin/time log outside the repository; only the small summary
## table belongs in a benchmark results directory.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("usage: benchmark-input-projection.R INPUT [SUMMARY.tsv]", call. = FALSE)
}
input <- args[[1L]]
if (!file.exists(input)) stop("input does not exist: ", input, call. = FALSE)

library(CompreSSoR)
started <- unname(proc.time()[["elapsed"]])
projected <- CompreSSoR:::read_sumstats_input(
  input, project_columns = TRUE, core_only = TRUE
)
elapsed <- unname(proc.time()[["elapsed"]]) - started
metadata <- attr(projected, "input_read_metadata")
summary <- data.frame(
  rows = nrow(projected),
  columns_before = metadata$columns_before,
  columns_read = metadata$columns_read,
  elapsed_seconds = elapsed,
  stringsAsFactors = FALSE
)

if (length(args) == 2L) {
  write.table(summary, args[[2L]], sep = "\t", quote = FALSE,
              row.names = FALSE)
} else {
  write.table(summary, stdout(), sep = "\t", quote = FALSE,
              row.names = FALSE)
}
