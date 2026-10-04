# read_candidates_batch() must return, for every store and in any batch
# composition, exactly what read_candidates() returns for that store alone.
#
# Regression: stage 2 kept one shared key slice per store in a list and
# assigned `shared[[i]] <- keys`.  For an identity group (stores with the same
# variant set) whose members had no candidate rows, `keys` was NULL, and
# `shared[[i]] <- NULL` deletes slot i, so every later store picked up another
# group's keys: its candidates were matched against the wrong rows and came
# back partial or empty, or the read failed with "subscript out of bounds".
# Batches of identical-panel stores (one group) never hit it; batches mixing
# variant sets with an all-empty group did (2,940 UKB-PPP stores: FIS1 1 of
# 367 flagged rows, EFNA4 0 of 246).

batch_rows_data <- function(set, seed, hits, n = 1200L) {
  set.seed(seed)
  half <- n / 2L
  # Each variant set has its own positions, so its own identity signature.
  offset <- match(set, c("x", "y", "z", "s1", "s2")) * 10000L
  pos <- rep(seq.int(100001L + offset, by = 3L, length.out = half), 2L)
  ref <- rep(c("C", "G", "T", "A"), length.out = n)
  alt <- rep(c("A", "C", "G", "T"), length.out = n)
  se <- exp(rnorm(n, log(0.02), 0.4))
  z <- pmax(pmin(rnorm(n), 3), -3)
  if (hits) z[sample.int(n, hits)] <- sample(c(-1, 1), hits, TRUE) * runif(hits, 6, 30)
  data.frame(
    chromosome = rep(c("1", "2"), each = half), base_pair_location = pos,
    reference_allele = ref, alternate_allele = alt,
    effect_allele = alt, other_allele = ref,
    beta = z * se, standard_error = se,
    effect_allele_frequency = runif(n, 0.05, 0.95), stringsAsFactors = FALSE
  )
}

batch_rows_stores <- function() {
  spec <- list(X1 = c("x", 0), X2 = c("x", 0), Y1 = c("y", 30), Y2 = c("y", 12),
               Z1 = c("z", 20), Z2 = c("z", 7), S1 = c("s1", 25), S2 = c("s2", 0))
  paths <- vapply(seq_along(spec), function(i) {
    s <- spec[[i]]
    out <- tempfile(paste0("batchrows-", names(spec)[i], "-"))
    compress_sumstats(batch_rows_data(s[[1]], 100L + i, as.integer(s[[2]])), out,
                      qc = "none", pvalue_flag = TRUE, pvalue_order = TRUE,
                      pvalue_order_threshold = 0.01, overwrite = TRUE)
    out
  }, character(1))
  names(paths) <- names(spec)
  # Two copies of the 0.4.5 legacy fixture: an older-format identity group
  # inside a mixed-format batch.  Each copy is a complete, self-contained store.
  legacy <- testthat::test_path("fixtures", "legacy-z9se6-a27b32d.cpr")
  for (nm in c("L1", "L2")) {
    dir <- tempfile(paste0("batchrows-", nm, "-"))
    dir.create(dir)
    file.copy(list.files(legacy, full.names = TRUE), dir)
    paths[[nm]] <- dir
  }
  paths
}

flag_hits <- function(path) {
  as.integer(open_compressor(path)$manifest$domains$pvalue_flag$hit_rows)
}

test_that("read_candidates_batch equals per-store reads in every composition", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  paths <- batch_rows_stores()
  sig <- vapply(paths, function(p) CompreSSoR:::pcodec_identity_signature(
    CompreSSoR:::pcodec_open_store_cached(p)), character(1))
  # The fixture really has distinct variant sets, shared within groups.
  expect_identical(sig[["X1"]], sig[["X2"]])
  expect_identical(sig[["L1"]], sig[["L2"]])
  expect_length(unique(sig), 6L)
  expect_identical(unname(vapply(paths, flag_hits, integer(1))),
                   c(0L, 0L, 30L, 12L, 20L, 7L, 25L, 0L, 75L, 75L))
  expect_identical(open_compressor(paths[["L1"]])$manifest$format_version,
                   "0.4.5-pcodec-native")
  expect_identical(open_compressor(paths[["X1"]])$manifest$format_version,
                   "0.4.6-pcodec-native")

  cols <- c("key", "p_value", "chromosome", "base_pair_location")
  single <- function(p, strategy, order) {
    read_candidates(p, 5e-8, columns = cols, order = order, strategy = strategy)
  }
  set.seed(1)
  compositions <- c(
    # Failed on main: the all-empty X group deletes a slot after the Y and Z
    # groups were assigned, shifting their keys onto later stores: the first
    # errors (subscript out of bounds), the second silently returns 0 of S1's
    # 25 flagged rows.
    list(c("Y1", "Z1", "X1", "X2", "S1", "Y2", "Z2"),
         c("Y1", "X1", "X2", "S1", "Z1", "Y2", "Z2"),
         c("L1", "X1", "X2", "S2", "L2", "S1"),
         names(paths)),
    replicate(6L, sample(names(paths), sample(2:10, 1L)), simplify = FALSE))
  for (comp in compositions) for (strategy in c("pvalue_flag", "auto")) {
    for (order in c("none", "exact")) for (threads in c(1L, 2L)) {
      got <- read_candidates_batch(as.list(paths[comp]), 5e-8, columns = cols,
                                   order = order, threads = threads,
                                   strategy = strategy)
      expect_identical(names(got), comp)
      for (nm in comp) {
        info <- paste(strategy, order, threads, paste(comp, collapse = ","), nm)
        expect_identical(got[[nm]], single(paths[[nm]], strategy, order), info = info)
        if (identical(strategy, "pvalue_flag")) {
          expect_identical(nrow(got[[nm]]), flag_hits(paths[[nm]]), info = info)
        }
      }
    }
  }
  # Per-store thresholds (one per store) follow the same path.
  comp <- c("Y1", "Z1", "X1", "X2", "S1", "Y2", "Z2")
  got <- read_candidates_batch(as.list(paths[comp]), rep(5e-8, length(comp)),
                               columns = cols, strategy = "pvalue_flag")
  for (nm in comp) expect_identical(got[[nm]], single(paths[[nm]], "pvalue_flag", "none"))
})

test_that("candidate reads refuse to return partial pvalue_flag data", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  path <- tempfile("batchrows-guard-")
  compress_sumstats(batch_rows_data("y", 7L, 30L), path, qc = "none",
                    pvalue_flag = TRUE, overwrite = TRUE)
  other <- tempfile("batchrows-guard-")
  compress_sumstats(batch_rows_data("z", 8L, 20L), other, qc = "none",
                    pvalue_flag = TRUE, overwrite = TRUE)
  ok <- read_candidates(path, 5e-8, strategy = "pvalue_flag")
  expect_identical(nrow(ok), 30L)

  # A flag stream that disagrees with the manifest's flagged-row count.
  real_flag_rows <- CompreSSoR:::candidates_read_flag_rows
  testthat::local_mocked_bindings(
    candidates_read_flag_rows = function(store, threads) {
      real_flag_rows(store, threads)[-1L]
    })
  expect_error(read_candidates(path, 5e-8, strategy = "pvalue_flag"),
               "manifest records 30 flagged rows")
  expect_error(read_candidates_batch(list(path, other), 5e-8, strategy = "pvalue_flag"),
               "manifest records")
})

test_that("candidate reads refuse to drop selected rows in the decode step", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  path <- tempfile("batchrows-guard-")
  compress_sumstats(batch_rows_data("y", 9L, 30L), path, qc = "none",
                    pvalue_flag = TRUE, overwrite = TRUE)
  real_fetch <- CompreSSoR:::candidates_fetch
  testthat::local_mocked_bindings(
    candidates_fetch = function(...) {
      x <- real_fetch(...)
      x[-nrow(x), , drop = FALSE]
    })
  expect_error(read_candidates(path, 5e-8, strategy = "pvalue_flag"), "rows lost")
  expect_error(read_candidates(path, 5e-8), "rows lost")
  # A region legitimately drops rows; that is not an error.
  expect_no_error(read_candidates(path, 5e-8, region = "chr1:1-200000000"))
})
