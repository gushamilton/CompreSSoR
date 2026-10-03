# Regenerates fixtures/legacy-z9se6-a27b32d.cpr and its reference reads.
# Run from tests/testthat with CompreSSoR a27b32d (the last Z9/SE6-only
# writer) first on .libPaths():
#   Rscript fixtures/make-legacy-z9se6.R
suppressMessages(library(CompreSSoR))
source("helper-legacy-store.R")
out <- file.path("fixtures", "legacy-z9se6-a27b32d.cpr")
unlink(out, recursive = TRUE)
compress_sumstats(legacy_store_fixture_data(), out, qc = "none",
                  pvalue_flag = TRUE, pvalue_order = TRUE,
                  pvalue_order_threshold = 0.05, overwrite = TRUE)
reads <- legacy_store_reads(out)
saveRDS(reads, file.path("fixtures", "legacy-z9se6-a27b32d-reads.rds"),
        version = 3L, compress = "xz")
cat("CompreSSoR", as.character(packageVersion("CompreSSoR")), "\n")
