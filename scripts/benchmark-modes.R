#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(CompreSSoR))

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "outputs/modes-edge-benchmark.csv"
set.seed(42)
n_reference_per_chr <- 1000L
repeats <- 25L

make_reference_chr <- function(chr) {
  pos <- seq.int(100001L, length.out = n_reference_per_chr)
  data.frame(
    chromosome = as.character(chr),
    base_pair_location = pos,
    effect_allele = rep(c("A", "C", "G", "T"), length.out = n_reference_per_chr),
    other_allele = rep(c("C", "G", "T", "A"), length.out = n_reference_per_chr),
    effect_allele_frequency = 0.1 + (seq_len(n_reference_per_chr) %% 80) / 100,
    variant_id = paste(chr, pos, sep = ":"),
    stringsAsFactors = FALSE
  )
}

reference <- do.call(rbind, lapply(1:22, make_reference_chr))
input <- reference[rep(seq_len(nrow(reference)), each = repeats), , drop = FALSE]
input$beta <- sin(seq_len(nrow(input)) / 31) / 5
input$standard_error <- 0.02 + (seq_len(nrow(input)) %% 17) / 1000
input$p_value <- 2 * pnorm(-abs(input$beta / input$standard_error))
input$annotation <- rep(c("a", "b", "c"), length.out = nrow(input))
row.names(input) <- NULL
panel <- reference[seq(1L, nrow(reference), by = 2L), c("variant_id"), drop = FALSE]

time_one <- function(expr) {
  unname(system.time(force(expr))[["elapsed"]])
}
results <- list()
for (replicate in seq_len(3L)) {
  serial <- NULL
  parallel <- NULL
  serial_seconds <- time_one(serial <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 1L))
  parallel_seconds <- time_one(parallel <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 4L))
  core_seconds <- time_one(core <- harmonise_sumstats(input, reference, mode = "core", variant_set = panel, chrom_threads = 4L))
  results[[length(results) + 1L]] <- data.frame(
    scenario = c("qc_serial", "qc_chrom_parallel", "core_chrom_parallel"),
    replicate = replicate,
    input_rows = nrow(input),
    output_rows = c(nrow(serial), nrow(parallel), nrow(core)),
    elapsed_seconds = c(serial_seconds, parallel_seconds, core_seconds),
    stringsAsFactors = FALSE
  )
  stopifnot(identical(serial$harmonisation_status, parallel$harmonisation_status))
  stopifnot(nrow(serial) == nrow(parallel), nrow(core) == nrow(input) / 2)
}
results <- do.call(rbind, results)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output, row.names = FALSE)
cat("wrote", output, "rows", nrow(results), "input_rows", nrow(input), "\n")
