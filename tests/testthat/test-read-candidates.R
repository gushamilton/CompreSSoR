candidate_fixture <- function(n = 3000L, seed = 11L) {
  set.seed(seed)
  x <- make_fixture(n)
  x$chromosome <- rep(c("1", "2"), each = n / 2L)
  x$base_pair_location <- rep(seq.int(100001L, length.out = n / 2L), 2L)
  x$standard_error <- 0.02 + (seq_len(n) %% 17) / 1000
  z <- rnorm(n, sd = 1.4)
  z[sample.int(n, 120L)] <- sample(c(-1, 1), 120L, TRUE) * runif(120L, 3, 40)
  x$beta <- z * x$standard_error
  x$p_value <- NULL
  x$beta[c(5L, 1700L)] <- NA_real_
  x$variant_id <- x$rsid <- x$annotation <- NULL
  x
}

candidate_store <- function(flag, order = FALSE) {
  compress_sumstats(candidate_fixture(), tempfile("cand-"), overwrite = TRUE,
                    qc = "none", pvalue_flag = flag,
                    pvalue_order = order, pvalue_order_threshold = 0.05)
}

reference_filter <- function(store, t, region = NULL) {
  full <- read_sumstats(store, columns = c("chromosome", "base_pair_location",
    "effect_allele", "other_allele", "p_value"))
  keep <- !is.na(full$p_value) & full$p_value <= t
  if (!is.null(region)) {
    b <- CompreSSoR:::read_region_bounds(region)
    keep <- keep & full$chromosome == b$chromosome &
      full$base_pair_location >= b$start & full$base_pair_location <= b$end
  }
  out <- full[keep, , drop = FALSE]
  row.names(out) <- NULL
  drop_row(out)
}

drop_row <- function(x) {
  x <- as.data.frame(x)
  x$row <- NULL
  attr(x, "candidate_strategy") <- NULL
  attr(x, "candidate_threshold") <- NULL
  attr(x, "candidate_order") <- NULL
  attr(x, "source_bytes_read") <- NULL
  x
}

test_that("read_candidates equals the filtered full read for every strategy", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  store <- candidate_store(flag = FALSE)
  semantic <- CompreSSoR:::candidates_semantic(
    CompreSSoR:::pcodec_open_store_cached(store$path))
  edge <- CompreSSoR:::candidates_central_edge_p(semantic)$p
  thresholds <- c(0.01, 1e-3, edge, nextafter_down(edge), nextafter_up(edge),
                  edge * (1 - 2e-9), 4.65e-4, 5e-8, 1e-300, 0, 1)
  regions <- list(NULL, "chr1:100101-100900", "chr2:100001-101500")
  for (t in thresholds) for (region in regions) {
    got <- read_candidates(store, t, region = region)
    expect_identical(drop_row(got), reference_filter(store, t, region))
    expect_true(attr(got, "candidate_strategy") %in% c("z_exceptions", "z_stream"))
  }
  expect_identical(attr(read_candidates(store, 1e-3), "candidate_strategy"), "z_stream")
  expect_identical(attr(read_candidates(store, 5e-8), "candidate_strategy"), "z_exceptions")
  expect_identical(attr(read_candidates(store, edge * (1 - 2e-9)),
                        "candidate_strategy"), "z_exceptions")
  expect_identical(attr(read_candidates(store, nextafter_down(edge)),
                        "candidate_strategy"), "z_stream")
  expect_identical(attr(read_candidates(store, edge), "candidate_strategy"), "z_stream")
})

test_that("auto never uses the flag; explicit pvalue_flag is opt-in", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  store <- candidate_store(flag = TRUE)
  got <- read_candidates(store, 5e-8)
  expect_identical(attr(got, "candidate_strategy"), "z_exceptions")
  expect_identical(drop_row(got), reference_filter(store, 5e-8))
  expect_identical(attr(read_candidates(store, 5e-8, strategy = "pvalue_flag"),
                        "candidate_strategy"), "pvalue_flag")
  expect_error(read_candidates(store, 1e-6, strategy = "pvalue_flag"), "requires")
  expect_error(read_candidates(store, 0.5, strategy = "exceptions"), "requires")
})

test_that("flag threshold with disagreeing supplied p matches full read", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  set.seed(1); n <- 400
  d <- data.frame(chromosome = "1", base_pair_location = seq(1e5, by = 10, length.out = n),
    reference_allele = "A", alternate_allele = "G", effect_allele = "G",
    other_allele = "A", beta = rnorm(n, 0, 0.01), standard_error = 0.01,
    effect_allele_frequency = 0.3)
  d$p_value <- 2 * pnorm(-abs(d$beta / d$standard_error))
  d$p_value[1:3] <- 1e-10
  d$beta[4] <- NA; d$p_value[4] <- 1e-12
  d$beta[5] <- 0.1; d$p_value[5] <- 0.2
  P <- tempfile(); suppressMessages(capture.output(compress_sumstats(d, P)))
  full <- read_sumstats(P, columns = "p_value")
  for (th in c(5e-8, 4.99999e-8, 5.00001e-8)) {
    r <- read_candidates(P, th, columns = "p_value")
    expect_identical(r$row, which(!is.na(full$p_value) & full$p_value <= th) - 1L)
  }
})

test_that("empty results, ordering and errors behave", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  store <- candidate_store(flag = TRUE, order = TRUE)
  small <- compress_sumstats(make_fixture(64L), tempfile("cand-"), overwrite = TRUE,
                             qc = "none", pvalue_flag = FALSE)
  none <- read_candidates(small, 1e-300, columns = c("chromosome", "p_value"))
  expect_identical(attr(none, "candidate_strategy"), "z_exceptions")
  expect_identical(nrow(none), 0L)
  expect_identical(names(none), c("row", "chromosome", "p_value"))
  empty_region <- read_candidates(store, 0.01, region = "chr1:1-2")
  expect_identical(nrow(empty_region), 0L)
  expect_error(read_candidates(store, 2), "between 0 and 1")

  rec <- read_candidates(store, 0.01, order = "reconstructed")
  expect_false(is.unsorted(rec$p_value))
  expect_identical(attr(rec, "candidate_order"), "reconstructed_p_approximate")

  ex <- read_candidates(store, 0.01, order = "exact")
  expect_true("exact_rank" %in% names(ex))
  ranked <- ex$exact_rank[!is.na(ex$exact_rank)]
  expect_false(is.unsorted(ranked))
  expected <- read_pvalue_order(store)
  expect_identical(ex$row[seq_along(ranked)], expected[expected %in% ex$row][seq_along(ranked)])
  expect_error(read_candidates(store, 0.5, order = "exact"), "covers only")
  no_order <- candidate_store(flag = FALSE, order = FALSE)
  expect_error(read_candidates(no_order, 0.01, order = "exact"), "no exact")
})

test_that("read_candidates_batch equals per-store results", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  stores <- list(a = candidate_store(TRUE), b = candidate_store(FALSE))
  for (t in c(5e-8, 1e-3)) {
    listed <- read_candidates_batch(stores, t)
    for (nm in names(stores)) {
      expect_identical(drop_row(listed[[nm]]),
                       drop_row(read_candidates(stores[[nm]], t)))
    }
    bound <- read_candidates_batch(stores, t, bind = TRUE, threads = 2L)
    expect_identical(unique(bound$store), c("a", "b"))
    expect_identical(sum(bound$store == "a"), nrow(listed$a))
    expect_identical(names(attr(bound, "candidate_strategy")), c("a", "b"))
  }
})

test_that("all projectable columns match the full read", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  store <- candidate_store(flag = FALSE)
  cols <- c("global_position", "substitution", "chromosome", "base_pair_location",
            "reference_allele", "alternate_allele", "effect_allele", "other_allele",
            "z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  full <- read_sumstats(store, columns = cols)
  got <- read_candidates(store, 0.01, columns = cols, threads = 2L)
  expect_identical(got$row + 1L, which(!is.na(full$p_value) & full$p_value <= 0.01))
  want <- full[got$row + 1L, , drop = FALSE]
  row.names(want) <- NULL
  for (nm in cols) expect_equal(got[[nm]], want[[nm]], tolerance = 1e-12)
  expect_error(read_candidates(store, 0.01, columns = "nope"), "unknown output")
})

test_that("store caches are keyed on content, not mtime", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  A <- candidate_fixture(seed = 3L)
  P <- tempfile("stale-"); Q <- tempfile("stale-")
  compress_sumstats(A, P, overwrite = TRUE, qc = "none")
  rows <- 0:(nrow(A) - 1L)
  r1 <- read_sumstats(P, variants = rows, columns = c("beta", "standard_error"))
  mt <- file.mtime(list.files(P, full.names = TRUE))
  names(mt) <- list.files(P)
  B <- A; B$beta <- -B$beta
  compress_sumstats(B, Q, overwrite = TRUE, qc = "none")
  expected <- read_sumstats(Q, variants = rows, columns = c("beta", "standard_error"))
  # Rewrite P in place and restore every original mtime.
  compress_sumstats(B, P, overwrite = TRUE, qc = "none")
  for (f in list.files(P, full.names = TRUE)) {
    if (basename(f) %in% names(mt)) Sys.setFileTime(f, mt[[basename(f)]])
  }
  r2 <- read_sumstats(P, variants = rows, columns = c("beta", "standard_error"))
  expect_identical(r2, expected)
  expect_false(identical(r2, r1))
})
