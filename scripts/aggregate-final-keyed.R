#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv(
  "COMPRESSOR_FINAL_RESULT_ROOT",
  "/user/work/fh6520/CompreSSoR-bp-thread-test/results/final-keyed-10m-compact-20260804"
)
pcodec_summary <- Sys.getenv(
  "COMPRESSOR_PCODEC_SUMMARY",
  "/user/work/fh6520/CompreSSoR-bp-thread-test/results/pcodec-10m-thread-compare-18268241-summary.csv"
)

files <- list.files(root, pattern = "final-keyed-run[.]csv$", recursive = TRUE,
                    full.names = TRUE)
stopifnot(length(files) > 0L)
x <- rbindlist(lapply(files, fread), fill = TRUE)
main <- x[status == "OK",
  .(runs = .N,
    source_bytes = median(source_bytes),
    storage_bytes = median(storage_bytes),
    compression_ratio = median(compression_ratio),
    write_median_seconds = median(write_seconds),
    read_median_seconds = median(full_read_seconds),
    rows = median(rows),
    status = "OK"),
  by = .(format_id, format, family)]

p <- fread(pcodec_summary)
p4 <- p[threads == 4]
stopifnot(nrow(p4) == 1L)
pcodec <- data.table(
  format_id = "pcodec_native",
  format = "CompreSSoR native Pcodec (4 threads)",
  family = "CompreSSoR",
  runs = p4$runs,
  source_bytes = unique(x$source_bytes),
  storage_bytes = p4$storage_bytes,
  compression_ratio = unique(x$source_bytes) / p4$storage_bytes,
  write_median_seconds = NA_real_,
  read_median_seconds = p4$median_seconds,
  rows = p4$rows,
  status = "OK"
)

summary <- rbind(main, pcodec, fill = TRUE)
summary[, pareto_frontier := vapply(seq_len(.N), function(i) {
  !any(read_median_seconds <= read_median_seconds[i] &
       compression_ratio >= compression_ratio[i] &
       (read_median_seconds < read_median_seconds[i] |
        compression_ratio > compression_ratio[i]))
}, logical(1))]

fwrite(summary, file.path(root, "whole-file-10m-summary.csv"))
fwrite(summary[pareto_frontier == TRUE], file.path(root, "whole-file-10m-frontier.csv"))

print(summary[order(read_median_seconds),
  .(format, family, runs, rows,
    storage_MB = round(storage_bytes / 1e6, 3),
    compression_ratio = round(compression_ratio, 3),
    read_s = round(read_median_seconds, 3), pareto_frontier)])
cat("\nVALID_MAIN_ROWS=", nrow(main),
    " TOTAL=", nrow(summary),
    " FRONTIER=", sum(summary$pareto_frontier), "\n", sep = "")
