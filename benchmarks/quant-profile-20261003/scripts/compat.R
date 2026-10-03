# Rscript compat.R LIB OUT.rds STORE...   reads every store with the CompreSSoR in LIB
a <- commandArgs(TRUE)
.libPaths(c(a[1], "/user/work/fh6520/r_packages", .libPaths()))
suppressMessages(library(CompreSSoR)); source("/user/work/fh6520/cs-quant/scripts/helper-legacy-store.R")
out <- list(version = as.character(packageVersion("CompreSSoR")), lib = find.package("CompreSSoR"))
for (s in a[-(1:2)]) {
  r <- legacy_store_reads(s)
  out[[s]] <- vapply(r, function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE), "")
}
saveRDS(out, a[2])
