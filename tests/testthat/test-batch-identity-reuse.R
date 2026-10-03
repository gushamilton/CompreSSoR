build_panel_store <- function(dir, name, n, seed, offset = 0L) {
  set.seed(seed)
  d <- make_fixture(n)
  d$base_pair_location <- as.integer(d$base_pair_location + offset)
  d$beta <- rnorm(n, 0, 0.05)
  d$standard_error <- runif(n, 0.01, 0.1)
  d$effect_allele_frequency <- runif(n, 0.05, 0.95)
  d$p_value <- 2 * pnorm(-abs(d$beta / d$standard_error))
  d$z <- d$beta / d$standard_error
  path <- file.path(dir, paste0(name, ".cpr"))
  compress_sumstats(d, path)
  list(path = path, data = d)
}

strip <- function(x) { attr(x, "source_bytes_read") <- NULL; x }

test_that("batched reads reuse identity resolution per panel group", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec unavailable")
  dir <- tempfile(); dir.create(dir)
  a1 <- build_panel_store(dir, "a1", 1200L, 1)
  a2 <- build_panel_store(dir, "a2", 1200L, 2)
  a3 <- build_panel_store(dir, "a3", 1200L, 3)
  b1 <- build_panel_store(dir, "b1", 1200L, 4, offset = 50000L)
  b2 <- build_panel_store(dir, "b2", 1200L, 5, offset = 50000L)
  paths <- c(a1 = a1$path, a2 = a2$path, b1 = b1$path, a3 = a3$path, b2 = b2$path)
  d <- a1$data
  keys <- compressor_variant_key(d$chromosome, d$base_pair_location,
                                 d$other_allele, d$effect_allele)
  mixed_keys <- c(keys[c(3, 50, 700, 3)], "1:5000:A:C")
  cols <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency")
  trace <- CompreSSoR:::.pcodec_batch_trace
  reference <- function(v = NULL, r = NULL, columns = cols)
    lapply(paths, function(p) strip(read_sumstats(p, variants = v, region = r, columns = columns)))

  check <- function(v = NULL, r = NULL, columns = cols, groups = 2L, threads = 1L,
                    share = TRUE) {
    old <- options(CompreSSoR.batch_share_panels = share)
    on.exit(options(old), add = TRUE)
    before <- trace$identity_resolutions
    got <- read_sumstats_batch(paths, v, columns = columns, threads = threads, region = r)
    expect_identical(trace$identity_resolutions - before, groups)
    expect_identical(lapply(got, strip), reference(v, r, columns))
    expect_identical(names(got), names(paths))
  }
  # Keys: absent + duplicate keys; panel A hits, panel B mostly misses.
  check(v = mixed_keys, groups = 2L)
  check(v = mixed_keys, threads = 2L)
  # Row ids.
  check(v = c(5L, 10L, 600L, 5L))
  # Region.
  check(r = "chr1:100100-100400")
  # Value-only columns, keys and region still resolved once per group.
  check(v = mixed_keys, columns = c("beta", "standard_error"))
  check(r = "chr1:100100-100400", columns = c("z", "p_value"))
  # Rows with value-only columns need no identity resolution at all.
  check(v = c(1L, 7L, 900L), columns = c("beta", "standard_error"), groups = 0L)
  # One pass per store (the default with several threads): one resolution per
  # store that needs identity work, same output.
  check(v = mixed_keys, groups = 5L, threads = 2L, share = FALSE)
  check(v = mixed_keys, groups = 5L, threads = 1L, share = FALSE)
  check(r = "chr1:100100-100400", groups = 5L, threads = 3L, share = FALSE)
  check(v = c(5L, 10L, 600L, 5L), groups = 5L, threads = 2L, share = FALSE)
  check(v = c(1L, 7L, 900L), columns = c("beta", "standard_error"), groups = 0L,
        threads = 2L, share = FALSE)
  check(v = mixed_keys, columns = c("z", "p_value"), groups = 5L, threads = 8L, share = FALSE)
  # Default strategy: panel sharing whatever the thread count.
  check(v = mixed_keys, groups = 2L, threads = 1L, share = NULL)
  check(v = mixed_keys, groups = 2L, threads = 2L, share = NULL)
  # A single store is a plain read with all threads.
  one <- read_sumstats_batch(paths[1], mixed_keys, columns = cols, threads = 2L)
  expect_identical(one[[1]], read_sumstats(paths[[1]], variants = mixed_keys, columns = cols))
  # Full reads.
  full <- read_sumstats_batch(paths[1:2], columns = cols, threads = 2L)
  expect_identical(lapply(full, strip),
                   lapply(paths[1:2], function(p) strip(read_sumstats(p, columns = cols))))
  # Errors raised while reading a store propagate in either strategy.
  bad_region <- tryCatch(read_sumstats(paths[[1]], region = "chrMT:1-10", columns = cols),
                         error = conditionMessage)
  expect_type(bad_region, "character")
  for (share in c(TRUE, FALSE)) {
    old <- options(CompreSSoR.batch_share_panels = share)
    expect_error(read_sumstats_batch(paths, columns = cols, region = "chrMT:1-10", threads = 2L),
                 bad_region, fixed = TRUE)
    options(old)
  }
  # Per-store request lists (different requests, same group still resolves each).
  per <- list(keys[1:3], keys[1:3], keys[5:6], keys[1:3], keys[5:6])
  got <- read_sumstats_batch(paths, per, columns = cols)
  expect_identical(lapply(got, strip),
                   Map(function(p, v) strip(read_sumstats(p, variants = v, columns = cols)),
                       paths, per))
  # Nothing matches.
  none <- read_sumstats_batch(paths[1:2], "1:5000:A:C", columns = cols)
  expect_identical(lapply(none, strip),
                   lapply(paths[1:2], function(p) strip(read_sumstats(p, variants = "1:5000:A:C", columns = cols))))
})
