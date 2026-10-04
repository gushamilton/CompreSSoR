# Every reader decodes a row to the same bits as the full read: row-ID, key,
# region and region+key reads, the per-block R path, threshold candidates
# (single and batched) and batched key reads all use the one native decoder.

identity_store <- local({
  path <- NULL
  function() {
    if (!is.null(path)) return(path)
    set.seed(5)
    n <- 200000L
    se <- exp(rnorm(n, log(0.02), 0.5))
    z <- rnorm(n) * 1.5
    z[sample.int(n, 2000L)] <- rnorm(2000L, 0, 9)
    eaf <- runif(n, 0.01, 0.99)
    eaf[sample.int(n, 4000L)] <- NA
    x <- data.frame(
      chromosome = rep(c("1", "2"), each = n / 2L),
      base_pair_location = rep(seq.int(10001L, by = 37L, length.out = n / 2L), 2L),
      reference_allele = "C", alternate_allele = "T", effect_allele = "T",
      other_allele = "C", beta = z * se, standard_error = se,
      effect_allele_frequency = eaf, stringsAsFactors = FALSE)
    path <<- tempfile("bit-identity-", fileext = ".cpr")
    compress_sumstats(x, path, overwrite = TRUE, threads = 2L)
    path
  }
})

bit_cols <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
              "z", "beta", "standard_error", "effect_allele_frequency", "p_value")

plain <- function(d, cols = setdiff(names(d), "row")) {
  d <- as.data.frame(d)[cols]
  attributes(d) <- list(names = cols, class = "data.frame",
                        row.names = .set_row_names(nrow(d)))
  d
}

test_that("row, key, region and per-block reads are bit-identical to the full read", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  p <- identity_store()
  full <- read_sumstats(p, columns = bit_cols)
  set.seed(9)
  rows <- sort(sample.int(nrow(full), 4000L)) - 1L
  keys <- paste(full$chromosome, full$base_pair_location, full$other_allele,
                full$effect_allele, sep = ":")
  ref <- plain(full[rows + 1L, ], bit_cols)
  for (th in c(1L, 3L)) {
    expect_identical(plain(read_sumstats(p, variants = rows, columns = bit_cols,
                                         threads = th), bit_cols), ref)
    expect_identical(plain(read_sumstats(p, variants = keys[rows + 1L],
                                         columns = bit_cols, threads = th), bit_cols), ref)
  }
  in_region <- full$chromosome == "2" & full$base_pair_location <= 900000
  expect_identical(plain(read_sumstats(p, region = "chr2:1-900000", columns = bit_cols),
                         bit_cols), plain(full[in_region, ], bit_cols))
  # region + keys takes the per-block R path
  both <- in_region & keys %in% keys[rows + 1L]
  expect_identical(plain(read_sumstats(p, region = "chr2:1-900000",
                                       variants = keys[rows + 1L], columns = bit_cols),
                         bit_cols), plain(full[both, ], bit_cols))
  old <- options(CompreSSoR.pcodec.native_select = FALSE)
  on.exit(options(old), add = TRUE)
  expect_identical(plain(read_sumstats(p, variants = rows, columns = bit_cols),
                         bit_cols), ref)
})

test_that("candidates and batched reads are bit-identical to the full read", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  p <- identity_store()
  full <- read_sumstats(p, columns = bit_cols)
  vcols <- c("z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  for (th in c(1e-3, 0.05)) {
    hit <- !is.na(full$p_value) & full$p_value <= th
    got <- read_candidates(p, th, columns = vcols)
    expect_identical(plain(got, vcols), plain(full[hit, ], vcols))
    reg <- read_candidates(p, th, region = "chr1:1-2000000", columns = vcols)
    expect_identical(plain(reg, vcols),
                     plain(full[hit & full$chromosome == "1" &
                                  full$base_pair_location <= 2000000, ], vcols))
    batch <- read_candidates_batch(list(p, p), th, columns = vcols, threads = 2L)
    expect_identical(plain(batch[[2]], vcols), plain(got, vcols))
  }
  keys <- paste(full$chromosome, full$base_pair_location, full$other_allele,
                full$effect_allele, sep = ":")[seq(1, nrow(full), by = 50)]
  b <- read_sumstats_batch(list(p, p), variants = keys, columns = bit_cols, threads = 2L)
  expect_identical(plain(b[[1]], bit_cols),
                   plain(full[seq(1, nrow(full), by = 50), ], bit_cols))
  expect_identical(plain(b[[2]], bit_cols), plain(b[[1]], bit_cols))
})

test_that("missing-EAF exception records carry only the EAF and decode unchanged", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  set.seed(3)
  n <- 20000L
  se <- exp(rnorm(n, log(0.02), 0.3))
  z <- rnorm(n)
  eaf <- runif(n, 0.05, 0.95)
  eaf[seq(1L, n, by = 2L)] <- NA
  x <- data.frame(chromosome = "1", base_pair_location = seq.int(1000L, by = 11L, length.out = n),
                  reference_allele = "A", alternate_allele = "G", effect_allele = "G",
                  other_allele = "A", beta = z * se, standard_error = se,
                  effect_allele_frequency = eaf, stringsAsFactors = FALSE)
  p <- tempfile(fileext = ".cpr")
  store <- compress_sumstats(x, p, overwrite = TRUE)
  e <- CompreSSoR:::pcodec_native_read_all_exceptions(store)
  eaf_only <- e$flags == 4L
  expect_gt(sum(eaf_only), n / 3)
  expect_true(all(e$z[eaf_only] == 0 & e$log2se[eaf_only] == 0))
  got <- read_sumstats(p, columns = c("z", "standard_error", "effect_allele_frequency"))
  expect_identical(is.na(got$effect_allele_frequency), is.na(eaf))
  expect_true(all(abs(got$standard_error / se - 1) < 0.05))
  expect_true(all(abs(got$z - z) < 0.01))
  # About 4 bytes per missing-EAF row (mostly the row id) instead of ~9
  # when the unused Z and SE fields held each row's own values.
  expect_lt(file.size(file.path(p, "exceptions.bin")), 6 * sum(is.na(eaf)))
})
