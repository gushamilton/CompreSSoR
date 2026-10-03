#!/usr/bin/env Rscript
# MR impact of the semantic profile: Rscript mr.R ARM   (ARM = exact | z9/eaf8/se6 | ...)
# Instruments are fixed (10x10 C8 arm of the showcase e2e, rep1) so only the
# stored values differ between arms.  The exact arm runs the same
# fast_mr_compressed() code with its store reader replaced by the exact TSV
# values at the same canonical keys.
ARM <- commandArgs(TRUE)[1]
Q <- "/user/work/fh6520/cs-quant"; SHOW <- "/user/work/fh6520/showcase"
.libPaths(c(file.path(Q, Sys.getenv("CSQ_LIB", "lib")), "/user/work/fh6520/r_packages",
            "/software/local/languages/miniforge3/envs/r-4.5.1/lib/R/library"))
suppressPackageStartupMessages({library(data.table); library(CompreSSoR); library(fastMR)})
source(file.path(SHOW, "prep/common.R"))     # convert_tsv_to_cpr(); setup_libs() is not called
TH <- 8L; setDTthreads(TH)
expo <- sprintf("exp%02d", 1:10); outc <- sprintf("out%02d", 1:10); traits <- c(expo, outc)
tsv <- setNames(file.path(SHOW, "prep/traits/tsv", paste0(traits, ".tsv.gz")), traits)
inst <- readRDS(file.path(SHOW, "e2e/results/rep1/10x10_C8/result.rds"))$instruments[expo]
methods <- c("ivw", "egger", "weighted_median", "simple_mode", "weighted_mode")
tag <- if (ARM == "exact") "exact" else gsub("/", "_", ARM)
res <- file.path(Q, "results", "mr"); dir.create(res, recursive = TRUE, showWarnings = FALSE)
info <- data.table(arm = ARM, host = Sys.info()[["nodename"]],
                   compressor = as.character(packageVersion("CompreSSoR")),
                   fastmr = as.character(packageVersion("fastMR")))
if (ARM == "exact") {
  keys_all <- unique(unlist(inst, use.names = FALSE))
  exact <- lapply(traits, function(id) {
    d <- fread(tsv[[id]], nThread = TH, select = c("chromosome", "base_pair_location", "effect_allele",
                                                   "other_allele", "beta", "standard_error"),
               colClasses = list(character = c("chromosome", "effect_allele", "other_allele")))
    d[, variant_key := CompreSSoR::compressor_variant_key(chromosome, base_pair_location, other_allele, effect_allele)]
    d <- d[variant_key %in% keys_all]
    as.data.frame(d)
  })
  names(exact) <- traits
  # Store paths only pass fastMR validation; values come from `exact`.
  cpr <- setNames(file.path(SHOW, "prep/traits/cpr", paste0(traits, ".cpr")), traits)
  path_to_id <- setNames(traits, unname(cpr))
  fake_io <- function(paths, keys, columns, io_threads) {
    out <- lapply(seq_along(paths), function(i) {
      d <- exact[[path_to_id[[paths[[i]]]]]]
      d <- d[d$variant_key %in% keys[[i]], , drop = FALSE]
      d[unique(c(columns, "variant_key"))]
    })
    attr(out, "source_bytes_read") <- NA_real_
    out
  }
  utils::assignInNamespace("fastmr_io_map", fake_io, ns = "fastMR")
} else {
  options(CompreSSoR.native_profile = ARM)
  work <- file.path(Q, "work", "mr", tag); dir.create(work, recursive = TRUE, showWarnings = FALSE)
  cpr <- setNames(file.path(work, paste0(traits, ".cpr")), traits)
  conv <- parallel::mclapply(traits, function(id) {
    tm <- convert_tsv_to_cpr(tsv[[id]], cpr[[id]], threads = 4L, tmpdir = work)
    stopifnot(identical(open_compressor(cpr[[id]])$manifest$semantic_codec$name, ARM))
    data.table(trait = id, prepare_s = tm[["prepare"]], encode_s = tm[["encode"]],
               bytes = sum(file.info(list.files(cpr[[id]], full.names = TRUE))$size))
  }, mc.cores = 4L, mc.preschedule = FALSE)
  bad <- !vapply(conv, is.data.frame, logical(1)); if (any(bad)) stop(conv[[which(bad)[1]]])
  fwrite(rbindlist(conv), file.path(res, paste0(tag, "_conversion.csv")))
}
t0 <- proc.time()[["elapsed"]]
mr <- fast_mr_compressed(cpr[expo], cpr[outc], inst, methods = methods, nboot = 1000, seed = 1,
                         threads = TH, io_threads = TH)
info$mr_s <- proc.time()[["elapsed"]] - t0
mr <- as.data.table(as.data.frame(mr)); mr[, arm := ARM]
fwrite(mr, file.path(res, paste0(tag, "_mr.csv")))
fwrite(info, file.path(res, paste0(tag, "_info.csv")))
print(info); print(mr[, .N, by = method])
cat("MR_OK\n")
