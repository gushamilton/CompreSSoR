#!/usr/bin/env Rscript
# One grid cell: Rscript grid.R DATASET PROFILE
#   DATASET: finngen | exp01 | out01 ...;  PROFILE: z9/eaf8/se6 ...
# 3 timed writes (from an in-memory prepared table, 8 threads), 3 timed full
# reads, then fidelity of the last store against the source.
a <- commandArgs(TRUE); DS <- a[1]; PROF <- a[2]
Q <- "/user/work/fh6520/cs-quant"
.libPaths(c(file.path(Q, Sys.getenv("CSQ_LIB", "lib")), "/user/work/fh6520/r_packages",
            "/software/local/languages/miniforge3/envs/r-4.5.1/lib/R/library"))
suppressPackageStartupMessages({library(data.table); library(CompreSSoR)})
TH <- 8L; setDTthreads(TH)
options(CompreSSoR.native_profile = PROF)
tag <- paste0(DS, "__", gsub("/", "_", PROF))
work <- file.path(Q, "work", "grid", tag); dir.create(work, recursive = TRUE, showWarnings = FALSE)
res <- file.path(Q, "results", "grid"); dir.create(res, recursive = TRUE, showWarnings = FALSE)
if (DS == "finngen") {
  d <- fread(cmd = paste("pigz -dc", shQuote("/user/work/fh6520/CompreSSoR-bp-thread-test/external/benchmark-10m/finngen_10m_snps.tsv.gz")),
             nThread = TH, colClasses = list(character = c("chrom", "ref", "alt")))
  src <- data.frame(chromosome = d$chrom, base_pair_location = d$pos,
                    reference_allele = d$ref, alternate_allele = d$alt,
                    effect_allele = d$alt, other_allele = d$ref, beta = d$beta,
                    standard_error = d$se, effect_allele_frequency = d$eaf,
                    p_value = d$p, stringsAsFactors = FALSE)
} else {
  d <- fread(file.path("/user/work/fh6520/showcase/prep/traits/tsv", paste0(DS, ".tsv.gz")),
             nThread = TH, colClasses = list(character = c("chromosome", "effect_allele", "other_allele")))
  d[, `:=`(reference_allele = other_allele, alternate_allele = effect_allele)]   # as convert_tsv_to_cpr()
  src <- as.data.frame(d)
}
rm(d); invisible(gc())
n <- nrow(src)
st <- file.path(work, "store.cpr")
wt <- numeric(3)
for (r in 1:3) {
  unlink(st, recursive = TRUE); invisible(gc())
  wt[r] <- system.time(compress_sumstats(src, st, input_build = "GRCh38", store_build = "GRCh38",
                                         threads = TH, overwrite = TRUE))[["elapsed"]]
}
files <- list.files(st, full.names = TRUE)
fsz <- setNames(file.info(files)$size, basename(files))
rt <- numeric(3); x <- NULL
for (r in 1:3) { rm(x); invisible(gc()); rt[r] <- system.time(x <- read_sumstats(st, threads = TH))[["elapsed"]] }
m <- open_compressor(st)$manifest
stopifnot(identical(m$semantic_codec$name, PROF), nrow(x) == n)
setDT(x)
kx <- paste(x$chromosome, x$base_pair_location, x$other_allele, x$effect_allele, sep = ":")
ks <- paste(src$chromosome, src$base_pair_location, src$other_allele, src$effect_allele, sep = ":")
mi <- match(ks, kx); stopifnot(!anyNA(mi)); x <- x[mi]; rm(kx, ks, mi); invisible(gc())
src_z <- src$beta / src$standard_error
src_p_exact <- 2 * pnorm(-abs(src_z))
err1 <- function(field, ref, new) {
  ad <- abs(new - ref); ok <- is.finite(ad)
  rel <- ad / abs(ref); rel[!is.finite(rel) | ref == 0] <- NA
  data.table(field = field, n = sum(ok), n_na_new = sum(!is.finite(new) & is.finite(ref)),
             max_abs = max(ad[ok]), p999_abs = quantile(ad[ok], 0.999, names = FALSE),
             max_rel = max(rel, na.rm = TRUE), p999_rel = quantile(rel, 0.999, na.rm = TRUE, names = FALSE))
}
errs <- rbindlist(list(
  err1("beta", src$beta, x$beta), err1("se", src$standard_error, x$standard_error),
  err1("z", src_z, x$z), err1("eaf", src$effect_allele_frequency, x$effect_allele_frequency),
  err1("nlp_vs_source_p", -log10(src$p_value), -log10(x$p_value)),
  err1("nlp_vs_exact_z", -log10(src_p_exact), -log10(x$p_value)),
  err1("beta_in_se_units", rep(0, n), (x$beta - src$beta) / src$standard_error)))
flips <- rbindlist(lapply(c("source_p", "exact_z_p"), function(ref) {
  pr <- if (ref == "source_p") src$p_value else src_p_exact
  rbindlist(lapply(c(5e-8, 1e-5, 0.01), function(th) {
    s <- !is.na(pr) & pr <= th; r <- !is.na(x$p_value) & x$p_value <= th
    data.table(reference = ref, threshold = th, n_ref = sum(s), n_recon = sum(r),
               lost = sum(s & !r), gained = sum(!s & r), flips = sum(s != r),
               flip_pct = 100 * sum(s != r) / max(1, sum(s)))
  }))
}))
meta <- data.table(dataset = DS, profile = PROF, n = n, bytes = sum(fsz),
                   mb = sum(fsz) / 1e6, bytes_per_variant = sum(fsz) / n,
                   z_bytes = fsz[["z.pco"]], se_bytes = fsz[["se.pco"]], eaf_bytes = fsz[["eaf.pco"]],
                   exc_bytes = fsz[["exceptions.bin"]],
                   ident_bytes = fsz[["position.pco"]] + fsz[["substitution.pco"]],
                   domain_bytes = sum(fsz[grepl("^pvalue_", names(fsz))]),
                   exception_rows = m$semantic_codec$exception_rows,
                   write_s_1 = wt[1], write_s_2 = wt[2], write_s_3 = wt[3], write_s_med = median(wt),
                   read_s_1 = rt[1], read_s_2 = rt[2], read_s_3 = rt[3], read_s_med = median(rt),
                   host = Sys.info()[["nodename"]], codec = m$codec$name)
fwrite(meta, file.path(res, paste0(tag, "_meta.csv")))
fwrite(cbind(dataset = DS, profile = PROF, errs), file.path(res, paste0(tag, "_errors.csv")))
fwrite(cbind(dataset = DS, profile = PROF, flips), file.path(res, paste0(tag, "_flips.csv")))
print(meta); print(errs); print(flips)
unlink(work, recursive = TRUE)
cat("CELL_OK\n")
