#!/usr/bin/env Rscript
# Usage: Rscript aggregate_csr_perf.R <results_dir_or_csvs...> [--baseline=baseline]
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(TRUE)
base <- sub("^--baseline=", "", grep("^--baseline=", args, value = TRUE)[1]); if (is.na(base)) base <- "baseline"
paths <- grep("^--", args, value = TRUE, invert = TRUE)
files <- unlist(lapply(paths, function(p) if (dir.exists(p)) list.files(p, "\\.csv$", full.names = TRUE) else p))
d <- rbindlist(lapply(files, fread), fill = TRUE)
s <- d[, .(n = .N, nodes = uniqueN(hostname), median_s = median(wall_s), min_s = min(wall_s),
           max_s = max(wall_s), peak_rss_mb = median(peak_rss_kb) / 1024,
           rows_out = rows_out[1], checksums = uniqueN(signif(checksum, 12))),
       by = .(op, label, git_sha = substr(git_sha, 1, 7))]
b <- s[label == base, .(op, base_s = median_s, base_rss = peak_rss_mb, base_rows = rows_out)]
s <- merge(s, b, by = "op", all.x = TRUE)
s[, `:=`(speedup = round(base_s / median_s, 2), rss_ratio = round(peak_rss_mb / base_rss, 2),
         rows_match = rows_out == base_rows)]
setorder(s, op, label)
print(s[, .(op, label, git_sha, n, nodes, median_s, min_s, max_s, peak_rss_mb = round(peak_rss_mb),
            speedup, rss_ratio, rows_match, checksums)], nrows = 500)
# Correctness gate: identical checksum across labels for every read op.
chk <- d[!grepl("^compress", op), .(k = uniqueN(signif(checksum, 12))), by = op]
if (any(chk$k > 1)) { cat("\nCHECKSUM MISMATCH:\n"); print(chk[k > 1]); quit(status = 1) }
