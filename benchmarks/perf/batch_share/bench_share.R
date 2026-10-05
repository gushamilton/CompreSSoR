#!/usr/bin/env Rscript
# read_sumstats_batch() over N stores with one shared key list: one pass per store vs panel sharing.
# Rscript bench_share.R LIB SET N K T MODE   (SET = showcase | ukbppp; MODE = one | share | auto)
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressMessages({library(CompreSSoR); library(data.table)})
SET <- a[2]; N <- as.integer(a[3]); K <- as.integer(a[4]); T <- as.integer(a[5]); MODE <- a[6]
if (MODE == "one") options(CompreSSoR.batch_share_panels = FALSE)
if (MODE == "share") options(CompreSSoR.batch_share_panels = TRUE)
if (SET == "showcase") {   # 20 traits (traits_v7) or 50 (traits25_v7), all on the same 1000G panel
  ids <- if (N <= 20) c(sprintf("exp%02d", 1:10), sprintf("out%02d", 1:10)) else c(sprintf("exp%02d", 1:25), sprintf("out%02d", 1:25))
  P <- file.path("/user/work/fh6520/showcase/prep", if (N <= 20) "traits_v7" else "traits25_v7", "cpr", paste0(ids, ".cpr"))[seq_len(N)]
  kf <- sprintf("/user/work/fh6520/showcase/storage/multi/lists/keys_%d.tsv", K)
} else if (SET == "mix") {   # the 20 showcase stores (one panel) plus N - 20 UKB-PPP stores (own panels)
  ids <- c(sprintf("exp%02d", 1:10), sprintf("out%02d", 1:10))
  all <- sort(list.files("/user/work/fh6520/pulling-proteins/results/production-pilot-20260806/stores", pattern = "\\.cpr$", full.names = TRUE))
  set.seed(1); P <- c(file.path("/user/work/fh6520/showcase/prep/traits_v7/cpr", paste0(ids, ".cpr")), sample(all, N - 20L))
  kf <- sprintf("/user/work/fh6520/showcase/storage/multi/lists/keys_%d.tsv", K)
} else {
  all <- sort(list.files("/user/work/fh6520/pulling-proteins/results/production-pilot-20260806/stores", pattern = "\\.cpr$", full.names = TRUE))
  set.seed(1); P <- if (N >= length(all)) all else sort(sample(all, N))
  kf <- sprintf("/user/work/fh6520/showcase/sharebench/ukb_keys_%d.tsv", K)
}
keys <- fread(kf, colClasses = "character"); v <- paste(keys$chrom, keys$pos, keys$ref, keys$alt, sep = ":")
cols <- c("chromosome", "base_pair_location", "other_allele", "effect_allele", "beta", "standard_error", "effect_allele_frequency", "p_value")
share_used <- tryCatch(CompreSSoR:::pcodec_batch_share_panels(length(P), T, TRUE, as.list(P)), error = function(e) NA)
t <- system.time(r <- read_sumstats_batch(P, variants = v, columns = cols, threads = T))[["elapsed"]]
dig <- digest::digest(lapply(r, function(x) { x <- as.data.frame(x); attributes(x) <- attributes(x)[c("names", "row.names", "class")]; x }), algo = "sha1")
cat(sprintf("RESULT,%s,%s,%d,%d,%d,%s,%s,%.3f,%d,%s,%s\n", packageVersion("CompreSSoR"), SET, length(P), K, T, MODE, share_used, t,
            sum(vapply(r, nrow, 0L)), dig, Sys.info()[["nodename"]]))
