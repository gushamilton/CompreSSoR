# Random keys from the first UKB-PPP store's variant panel (K = 1000, 6000, 100000).
.libPaths(c(commandArgs(TRUE)[1], .libPaths())); suppressMessages({library(CompreSSoR); library(data.table)})
s <- sort(list.files("/user/work/fh6520/pulling-proteins/results/production-pilot-20260806/stores", pattern = "\\.cpr$", full.names = TRUE))[1]
x <- as.data.table(read_sumstats(s, columns = c("chromosome", "base_pair_location", "reference_allele", "alternate_allele"), threads = 4))
for (k in c(1000L, 6000L, 100000L)) { set.seed(20261005L + k); y <- x[sort(sample.int(.N, k))]
  fwrite(y[, .(chrom = chromosome, pos = base_pair_location, ref = reference_allele, alt = alternate_allele)], sprintf("/user/work/fh6520/showcase/sharebench/ukb_keys_%d.tsv", k), sep = "\t") }
