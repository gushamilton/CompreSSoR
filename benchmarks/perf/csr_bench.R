#!/usr/bin/env Rscript
# CompreSSoR hot-path benchmark (one operation per process so the wrapper's
# /usr/bin/time peak RSS is attributable to that operation).
#
# Usage:
#   Rscript csr_bench.R build  <store_dir> <n_rows> <traits> <out.csv> [threads]
#   Rscript csr_bench.R op     <store_dir> <op_name> <reps> <out.csv>
#   Rscript csr_bench.R list-ops
#
# Env: CSR_LABEL (baseline/optimised), CSR_GIT_SHA (else git -C $CSR_REPO),
#      CSR_REPO. Emits one CSV row per rep: env columns + op, n_rows, traits,
#      rep, wall_s, rows_out, checksum, r_max_mb (R heap high-water mark).
#      The shell wrapper appends peak_rss_kb.
#
# Ops replicate how fastMR (R/compressed.R, R/clumping.R) calls CompreSSoR, so
# baseline and optimised SHAs can be compared on identical public calls.

here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
source(file.path(here, "csr_bench_common.R"))

MR <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error")
ID <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")

store_paths <- function(dir) {
  p <- sort(list.files(dir, pattern = "^trait[0-9]+\\.cpr$", full.names = TRUE))
  if (!length(p)) stop("no trait*.cpr stores in ", dir)
  p
}

# Deterministic query sets derived from the store itself (seeded), cached.
query_sets <- function(dir) {
  cache <- file.path(dir, "queries.rds")
  if (file.exists(cache)) return(readRDS(cache))
  st <- open_compressor(store_paths(dir)[1])
  n <- as.integer(st$manifest$n_rows)
  set.seed(42)
  rows1k <- sort(sample.int(n, 1000L) - 1L)
  rows100k <- sort(sample.int(n, 100000L) - 1L)
  rows10k <- sort(sample.int(n, 10000L) - 1L)
  idk <- read_sumstats(st, variants = rows10k, columns = ID)
  keys10k <- compressor_variant_key(idk$chromosome, idk$base_pair_location,
                                    idk$other_allele, idk$effect_allele)
  # cis windows: 50 random gene-like loci, +/-300 kb around a 30 kb body
  chr <- sample(1:22, 50, replace = TRUE, prob = csr_grch38_lengths)
  start <- vapply(chr, function(c) sample.int(csr_grch38_lengths[c] - 2e6, 1) + 5e5, 0)
  cis <- sprintf("chr%d:%d-%d", chr, as.integer(start - 3e5), as.integer(start + 3e4 + 3e5))
  q <- list(n = n, rows1k = rows1k, rows100k = rows100k,
            keys1k = keys10k[sort(sample.int(10000L, 1000L))], keys10k = keys10k, cis = cis)
  saveRDS(q, cache)
  q
}

ops <- list(
  # ---- whole-store decode ---------------------------------------------------
  full_mr_t1 = function(d, q) read_sumstats(store_paths(d)[1], columns = MR, threads = 1),
  full_mr_t4 = function(d, q) read_sumstats(store_paths(d)[1], columns = MR, threads = 4),
  full_p_t4  = function(d, q) read_sumstats(store_paths(d)[1], columns = "p_value", threads = 4),
  # ---- row / key subsets (fast_mr_compressed path uses canonical keys) ------
  rows_1k    = function(d, q) read_sumstats(store_paths(d)[1], variants = q$rows1k, columns = MR),
  rows_100k  = function(d, q) read_sumstats(store_paths(d)[1], variants = q$rows100k, columns = MR),
  keys_1k    = function(d, q) read_sumstats(store_paths(d)[1], variants = q$keys1k, columns = MR),
  keys_10k   = function(d, q) read_sumstats(store_paths(d)[1], variants = q$keys10k, columns = MR),
  # ---- regional (fastMR#4/#8 cis: region read then p<=0.01 filter) ---------
  cis50_p01  = function(d, q) {
    p <- store_paths(d)[1]
    out <- lapply(q$cis, function(r) {
      x <- read_sumstats(p, region = r, columns = c(ID, "p_value"))
      x[is.finite(x$p_value) & x$p_value <= 0.01, , drop = FALSE]
    })
    do.call(rbind, out)
  },
  # ---- genome-wide candidate extraction (fast_clump_compressed) -------------
  # fastMR candidate_source="pvalue_flag": flag rows -> ID + p
  cand_flag  = function(d, q) {
    p <- store_paths(d)[1]
    rows <- read_pvalue_flag(p)
    read_sumstats(p, variants = rows, columns = c(ID, "p_value"))
  },
  # fastMR candidate_source="full": all row IDs passed as variants (current code)
  cand_full_rowids = function(d, q) {
    p <- store_paths(d)[1]
    x <- read_sumstats(p, variants = seq.int(0L, q$n - 1L), columns = c(ID, "p_value"),
                       threads = 1L)
    x[is.finite(x$p_value) & x$p_value <= 5e-8, , drop = FALSE]
  },
  # best current public path for an arbitrary threshold: p-only decode, filter, fetch
  cand_p01_decode_filter = function(d, q) {
    p <- store_paths(d)[1]
    pv <- read_sumstats(p, columns = "p_value", threads = 4)$p_value
    rows <- which(is.finite(pv) & pv <= 0.01) - 1L
    read_sumstats(p, variants = rows, columns = c(ID, "p_value"))
  },
  order_exact = function(d, q) read_pvalue_order(store_paths(d)[1]),
  # ---- multi-trait (fast_mr_compressed / fast_clump_compressed over K stores)
  batch_keys1k_t1 = function(d, q) read_sumstats_batch(store_paths(d), q$keys1k, columns = MR, threads = 1),
  batch_keys1k_t4 = function(d, q) read_sumstats_batch(store_paths(d), q$keys1k, columns = MR, threads = 4),
  multi_cand_flag = function(d, q) lapply(store_paths(d), function(p) {
    read_sumstats(p, variants = read_pvalue_flag(p), columns = c(ID, "p_value"))
  }),
  multi_cis50_p01 = function(d, q) lapply(store_paths(d), function(p) {
    do.call(rbind, lapply(q$cis, function(r) {
      x <- read_sumstats(p, region = r, columns = c(ID, "p_value"))
      x[is.finite(x$p_value) & x$p_value <= 0.01, , drop = FALSE]
    }))
  })
)

emit <- function(rows, out) {
  rows <- rbindlist(rows, fill = TRUE)
  fwrite(rows, out, append = file.exists(out))
}

args <- commandArgs(TRUE)
mode <- args[1]
if (identical(mode, "list-ops")) { cat(names(ops), sep = "\n"); quit(save = "no") }
if (identical(mode, "build")) {
  dir <- args[2]; n <- as.numeric(args[3]); traits <- as.integer(args[4])
  out <- args[5]; threads <- as.integer(if (length(args) >= 6) args[6] else 4L)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (k in seq_len(traits)) {
    d <- as.data.frame(csr_synthetic_gwas(n, seed = 1L, trait = k))
    gc(reset = TRUE)
    t <- system.time(compress_sumstats(d, file.path(dir, sprintf("trait%02d.cpr", k)),
                                       qc = "none", overwrite = TRUE,
                                       pvalue_order = TRUE, threads = threads))[["elapsed"]]
    m <- open_compressor(file.path(dir, sprintf("trait%02d.cpr", k)))$manifest
    rows[[k]] <- cbind(csr_env_row(), data.table(
      op = "compress_none", n_rows = n, traits = 1L, rep = k, wall_s = t,
      rows_out = m$n_rows, checksum = csr_dir_bytes(file.path(dir, sprintf("trait%02d.cpr", k))),
      r_max_mb = { g <- gc(); sum(g[, ncol(g)]) }, threads = threads))
    rm(d); gc()
  }
  emit(rows, out)
  invisible(query_sets(dir))
  quit(save = "no")
}
if (identical(mode, "op")) {
  dir <- args[2]; op <- args[3]; reps <- as.integer(args[4]); out <- args[5]
  if (!op %in% names(ops)) stop("unknown op ", op)
  q <- query_sets(dir)
  traits <- length(store_paths(dir))
  invisible(ops[[op]](dir, q))  # warm-up (page cache, lazy loading); not recorded
  rows <- list()
  for (r in seq_len(reps)) {
    gc(reset = TRUE)
    t0 <- proc.time()[["elapsed"]]
    v <- ops[[op]](dir, q)
    t <- proc.time()[["elapsed"]] - t0
    nout <- if (is.data.frame(v)) nrow(v) else if (is.list(v)) sum(vapply(v, NROW, 0)) else length(v)
    rows[[r]] <- cbind(csr_env_row(), data.table(
      op = op, n_rows = q$n, traits = traits, rep = r, wall_s = t,
      rows_out = nout, checksum = csr_checksum(v), r_max_mb = { g <- gc(); sum(g[, ncol(g)]) },
      threads = NA_integer_))
    rm(v)
  }
  emit(rows, out)
  quit(save = "no")
}
stop("mode must be build, op or list-ops")
