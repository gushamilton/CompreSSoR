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
  b1 <- build_panel_store(dir, "b1", 1100L, 4, offset = 50000L)
  b2 <- build_panel_store(dir, "b2", 1100L, 5, offset = 50000L)
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
  # Default strategy. 5 stores on 2 panels (repeat fraction 0.6, below 0.9):
  # one pass per store.
  check(v = mixed_keys, groups = 5L, threads = 1L, share = NULL)
  check(v = mixed_keys, groups = 5L, threads = 2L, share = NULL)
  # 10 stores on one panel (repeat fraction 0.9), one shared key list: panel sharing.
  more <- lapply(4:10, function(i) build_panel_store(dir, paste0("a", i), 1200L, i + 2L)$path)
  apaths <- c(a1 = a1$path, a2 = a2$path, a3 = a3$path, setNames(unlist(more), paste0("a", 4:10)))
  old <- options(CompreSSoR.batch_share_panels = NULL)
  for (th in c(1L, 2L)) {
    before <- trace$identity_resolutions
    got <- read_sumstats_batch(apaths, mixed_keys, columns = cols, threads = th)
    expect_identical(trace$identity_resolutions - before, 1L)
    expect_identical(lapply(got, strip),
                     lapply(apaths, function(p) strip(read_sumstats(p, variants = mixed_keys, columns = cols))))
    # regions and row ids keep one pass per store
    before <- trace$identity_resolutions
    got <- read_sumstats_batch(apaths, region = "chr1:100100-100400", columns = cols, threads = th)
    expect_identical(trace$identity_resolutions - before, 10L)
  }
  options(old)
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

test_that("the default batch strategy rule", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec unavailable")
  dir <- tempfile(); dir.create(dir)
  a1 <- build_panel_store(dir, "a1", 600L, 1)
  a2 <- build_panel_store(dir, "a2", 600L, 2)
  b1 <- build_panel_store(dir, "b1", 500L, 4, offset = 50000L)
  rule <- CompreSSoR:::pcodec_batch_share_panels
  per <- CompreSSoR:::PCODEC_BATCH_SHARE_STORES_PER_THREAD
  old <- options(CompreSSoR.batch_share_panels = NULL)
  on.exit(options(old), add = TRUE)
  same <- as.list(rep(c(a1$path, a2$path), 5L))           # one panel, 10 stores
  mixed <- list(a1$path, b1$path)                         # two panels
  expect_equal(CompreSSoR:::pcodec_batch_panel_repeat_fraction(same), 0.9)
  expect_false(rule(4L, 1L, TRUE, same[1:4]))             # 0.75: too few repeats
  expect_equal(CompreSSoR:::pcodec_batch_panel_repeat_fraction(mixed), 0)
  expect_true(rule(10L, 1L, TRUE, same))
  expect_false(rule(10L, 1L, FALSE, same))                 # not one shared key list
  expect_false(rule(2L, 1L, TRUE, mixed))                 # panels not shared
  expect_false(rule(1L, 8L, TRUE, same[1L]))              # a single store
  many <- rep(same, length.out = per * 2L + 1L)
  expect_false(rule(length(many), 2L, TRUE, many))        # too many stores per thread
  expect_true(rule(length(many), 3L, TRUE, many))
  expect_false(rule(2L, 1L, TRUE, list("/nonexistent/a.cpr", "/nonexistent/a.cpr")))
  # an explicit option overrides the rule either way
  options(CompreSSoR.batch_share_panels = TRUE)
  expect_true(rule(2L, 1L, TRUE, mixed))
  options(CompreSSoR.batch_share_panels = FALSE)
  expect_false(rule(10L, 1L, TRUE, same))
})
