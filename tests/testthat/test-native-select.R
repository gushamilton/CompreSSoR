make_select_store <- function(dir, n = 3000L, seed = 11, name = "s.cpr") {
  set.seed(seed)
  d <- make_fixture(n)
  d$chromosome <- rep(c("1", "2", "X"), each = ceiling(n / 3))[seq_len(n)]
  d$base_pair_location <- as.integer(100000L + seq_len(n) * 7L)
  d$standard_error <- runif(n, 0.01, 0.1)
  d$beta <- rnorm(n, 0, 0.05)
  big <- sample.int(n, 120)
  d$beta[big] <- d$standard_error[big] * sample(c(-1, 1), 120, TRUE) * runif(120, 3.6, 9)
  d$effect_allele_frequency <- runif(n, 0.02, 0.98)
  d$effect_allele_frequency[sample.int(n, 40)] <- 0.5
  d$z <- d$beta / d$standard_error
  d$p_value <- 2 * pnorm(-abs(d$z))
  path <- file.path(dir, name)
  suppressWarnings(compress_sumstats(d, path))
  list(path = path, data = d)
}

test_that("native select is identical to the R path", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec unavailable")
  skip_if_not(is.loaded("compressor_read_pcodec_native_select", PACKAGE = "CompreSSoR"))
  dir <- tempfile(); dir.create(dir)
  st <- make_select_store(dir)
  d <- st$data
  n <- nrow(d)
  both <- function(..., threads = 1L) {
    op <- options(CompreSSoR.pcodec.native_select = TRUE); on.exit(options(op))
    a <- read_sumstats(st$path, ..., threads = threads)
    options(CompreSSoR.pcodec.native_select = FALSE)
    b <- read_sumstats(st$path, ..., threads = threads)
    expect_identical(a, b)
    a
  }
  cols <- list(NULL,
               c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error"),
               c("z", "p_value"), c("beta", "standard_error"), "standard_error",
               "effective_allele_frequency_placeholder")
  cols[[6]] <- "effect_allele_frequency"
  rows_sets <- list(
    c(0L, 5L, 2999L), sort(sample.int(n, 500) - 1L), rev(c(10L, 11L, 12L, 3L, 3L)),
    0:(n - 1L), integer())
  for (th in 1:2) for (cc in cols) {
    for (r in rows_sets) both(variants = r, columns = cc, threads = th)
  }
  # exception-bearing rows must be exercised
  exc_rows <- which(abs(d$beta / d$standard_error) >= 3.6) - 1L
  expect_gt(length(exc_rows), 20L)
  both(variants = exc_rows, columns = c("z", "beta", "standard_error", "effect_allele_frequency"))
  # regions: inclusive edges, empty, wide, other chromosomes
  p <- d$base_pair_location
  for (th in 1:2) for (cc in cols) {
    both(region = sprintf("chr1:%d-%d", p[10], p[40]), columns = cc, threads = th)
  }
  both(region = sprintf("chr1:%d-%d", p[10], p[10]), columns = "z")
  both(region = "chr1:1-5", columns = c("chromosome", "z"))
  both(region = "chr1:100000-100000000", columns = NULL)
  both(region = "chrX:1-100000000", columns = c("base_pair_location", "beta"))
  both(region = "chr2:1-100000000")
  # keys: duplicates, absent, mixed
  keys <- compressor_variant_key(d$chromosome, d$base_pair_location, d$other_allele, d$effect_allele)
  for (th in 1:2) for (cc in cols[1:4]) {
    both(variants = c(keys[c(3, 900, 1500, 3)], "1:5000:A:C"), columns = cc, threads = th)
  }
  both(variants = keys[exc_rows[1:15] + 1L], columns = c("z", "chromosome"))
  both(variants = "1:5000:A:C", columns = "z")
  both(variants = sample(keys, 800), columns = c("beta", "base_pair_location"))
})

test_that("native select handles NA values", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec unavailable")
  dir <- tempfile(); dir.create(dir)
  d <- make_fixture(1500L)
  d$beta[c(5, 600, 1400)] <- NA
  d$standard_error[c(7, 601)] <- NA
  d$effect_allele_frequency[c(9, 602)] <- NA
  d$z <- d$beta / d$standard_error
  d$p_value <- 2 * pnorm(-abs(d$z))
  path <- file.path(dir, "na.cpr")
  ok <- tryCatch({ suppressWarnings(compress_sumstats(d, path)); TRUE }, error = function(e) FALSE)
  skip_if_not(ok, "store with NA values cannot be written")
  rows <- c(4L, 5L, 6L, 8L, 599L, 600L, 601L, 1400L)
  op <- options(CompreSSoR.pcodec.native_select = TRUE); on.exit(options(op))
  a <- read_sumstats(path, variants = rows)
  options(CompreSSoR.pcodec.native_select = FALSE)
  b <- read_sumstats(path, variants = rows)
  expect_identical(a, b)
})

test_that("index cache invalidates when a store is rewritten in place", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec unavailable")
  dir <- tempfile(); dir.create(dir)
  s1 <- make_select_store(dir, n = 2500L, seed = 1)
  r1 <- read_sumstats(s1$path, variants = c(0L, 100L, 2000L), columns = c("z", "p_value"))
  r1b <- read_sumstats(s1$path, variants = c(0L, 100L, 2000L), columns = c("z", "p_value"))
  expect_identical(r1, r1b)
  unlink(s1$path, recursive = TRUE)
  Sys.sleep(0.05)
  s2 <- make_select_store(dir, n = 1800L, seed = 2)
  r2 <- read_sumstats(s2$path, variants = c(0L, 100L, 1700L), columns = c("z", "p_value"))
  op <- options(CompreSSoR.pcodec.native_select = FALSE); on.exit(options(op))
  expect_identical(r2, read_sumstats(s2$path, variants = c(0L, 100L, 1700L),
                                     columns = c("z", "p_value")))
  expect_error(read_sumstats(s2$path, variants = 2000L))
  expect_equal(nrow(read_sumstats(s2$path, columns = "z")), 1800L)
})

test_that("native select is identical across multiple blocks", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec unavailable")
  dir <- tempfile(); dir.create(dir)
  n <- 140000L
  st <- make_select_store(dir, n = n, seed = 5)
  d <- st$data
  both <- function(..., threads = 1L) {
    op <- options(CompreSSoR.pcodec.native_select = TRUE); on.exit(options(op))
    a <- read_sumstats(st$path, ..., threads = threads)
    options(CompreSSoR.pcodec.native_select = FALSE)
    expect_identical(a, read_sumstats(st$path, ..., threads = threads))
    invisible(a)
  }
  set.seed(3)
  rows <- sort(c(sample.int(n, 300) - 1L, 65535L, 65536L, 131071L, 131072L, 0L, n - 1L))
  keys <- compressor_variant_key(d$chromosome, d$base_pair_location, d$other_allele, d$effect_allele)
  p <- d$base_pair_location
  for (th in 1:2) {
    both(variants = rows, threads = th)
    both(variants = rows, columns = c("beta", "standard_error"), threads = th)
    both(variants = keys[rows + 1L], columns = c("chromosome", "base_pair_location", "z"), threads = th)
    both(region = sprintf("chr1:%d-%d", p[65000], p[66000]), threads = th)
    both(region = sprintf("chr1:%d-%d", p[65536], p[65537]), columns = "p_value", threads = th)
    both(region = sprintf("chr2:%d-%d", p[60000], p[100000]), columns = "z", threads = th)
  }
})
