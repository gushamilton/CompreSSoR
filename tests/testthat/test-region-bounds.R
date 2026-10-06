region_fixture <- function() {
  lengths <- CompreSSoR:::compressor_chromosome_lengths("GRCh38")
  # Rows at both ends of chromosomes 2, 3 and 4 so that an unclamped region
  # end (or start) on chromosome 3 would reach the neighbouring chromosome.
  chromosome <- c("2", "2", "3", "3", "3", "4", "4", "4")
  position <- c(1L, as.integer(lengths[["2"]]), 1L, 5000L,
                as.integer(lengths[["3"]]), 1L, 10253L, 20000L)
  n <- length(position)
  input <- make_fixture(n)
  input$chromosome <- chromosome
  input$base_pair_location <- position
  input$variant_id <- paste(chromosome, position, input$other_allele,
                            input$effect_allele, sep = ":")
  input
}

test_that("region ranges are clamped to the chromosome span", {
  lengths <- CompreSSoR:::compressor_chromosome_lengths("GRCh38")
  offsets <- CompreSSoR:::compressor_chromosome_offsets("GRCh38")
  range <- CompreSSoR:::pcodec_native_region_range
  chr3 <- as.numeric(offsets[["3"]])
  len3 <- as.numeric(lengths[["3"]])
  expect_identical(range("chr3:197590791-198590791", "GRCh38"),
                   c(chr3 + 197590790, chr3 + len3 - 1))
  expect_identical(range("chr3:0-10", "GRCh38"), c(chr3, chr3 + 9))
  expect_identical(range(c(3, -100, 10), "GRCh38"), c(chr3, chr3 + 9))
  expect_identical(range("chr3:1-1", "GRCh38"), c(chr3, chr3))
  expect_null(range("chr3:300000000-400000000", "GRCh38"))
  expect_null(range("chr3:200-100", "GRCh38"))
  expect_null(range(c(3, -10, 0), "GRCh38"))
  expect_error(range("chrZ:1-10", "GRCh38"), "unsupported region chromosome")
  # Every chromosome, including the last, is clamped to its own span.
  for (chromosome in names(lengths)) {
    got <- range(sprintf("chr%s:1-%.0f", chromosome, 1e12), "GRCh38")
    expect_identical(got, as.numeric(offsets[[chromosome]]) +
                       c(0, as.numeric(lengths[[chromosome]]) - 1))
  }
})

test_that("a region past the chromosome end never returns the next chromosome", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- region_fixture()
  store <- compress_sumstats(input, tempfile("region-bounds-"), overwrite = TRUE)
  len3 <- as.numeric(CompreSSoR:::compressor_chromosome_lengths("GRCh38")[["3"]])
  columns <- c("chromosome", "base_pair_location", "z")

  for (native in c(TRUE, FALSE)) {
    old <- options(CompreSSoR.pcodec.native_select = native)
    on.exit(options(old), add = TRUE)

    past_end <- read_sumstats(store, region = "chr3:197590791-198590791",
                              columns = columns)
    expect_identical(past_end$chromosome, "3")
    expect_identical(past_end$base_pair_location, as.integer(len3))

    far_past <- read_sumstats(store, region = "chr3:1-4000000000", columns = columns)
    expect_identical(far_past$chromosome, rep("3", 3L))

    beyond <- read_sumstats(store, region = "chr3:300000000-400000000",
                            columns = columns)
    expect_identical(nrow(beyond), 0L)

    from_zero <- read_sumstats(store, region = "chr3:0-5000", columns = columns)
    expect_identical(from_zero$chromosome, c("3", "3"))
    expect_identical(from_zero$base_pair_location, c(1L, 5000L))

    negative <- read_sumstats(store, region = c(3, -1000, 5000), columns = columns)
    expect_identical(negative$base_pair_location, c(1L, 5000L))
    expect_true(all(negative$chromosome == "3"))

    chr4_start <- read_sumstats(store, region = "chr4:0-10253", columns = columns)
    expect_identical(chr4_start$chromosome, c("4", "4"))

    last <- read_sumstats(store, region = "chr2:1-999999999", columns = columns)
    expect_identical(last$chromosome, c("2", "2"))
    options(old)
  }

  batch <- read_sumstats_batch(list(store, store), region = "chr3:197590791-198590791",
                               columns = columns)
  for (part in batch) expect_identical(part$chromosome, "3")

  candidates <- read_candidates(store, region = "chr3:197590791-198590791",
                                pvalue_threshold = 1)
  if (nrow(candidates)) expect_true(all(candidates$chromosome == "3"))
})

test_that("region strings accept commas, scientific notation and chromosomes 23/24", {
  bounds <- CompreSSoR:::read_region_bounds
  expect_identical(bounds("1:1-1e6")[c("start", "end")], list(start = 1, end = 1e6))
  expect_identical(bounds("chr1:1-1,000,000")$end, 1e6)
  expect_identical(bounds("chr1:1,000-2.5e3")[c("start", "end")],
                   list(start = 1000, end = 2500))
  expect_identical(bounds("23:1-100000")$chromosome, "X")
  expect_identical(bounds("chr24:1-5")$chromosome, "Y")
  expect_identical(bounds(c("chr23", "1", "1e3"))$chromosome, "X")
  expect_identical(bounds(c(24, 1, 5))$chromosome, "Y")
  expect_error(bounds("1:1-1.5"), "whole number")
  expect_error(bounds("1:1-1e-3"), "whole number")
  expect_error(bounds("1:1-1,00"), "whole number")
  expect_error(bounds("1:5--3"), "whole number")
  expect_error(bounds("1:a-3"), "whole number")
  expect_error(bounds("chr1-100"), "region must look like")

  range <- CompreSSoR:::pcodec_native_region_range
  offsets <- CompreSSoR:::compressor_chromosome_offsets("GRCh38")
  expect_identical(range("23:1-100000", "GRCh38"),
                   as.numeric(offsets[["X"]]) + c(0, 99999))
  expect_identical(range("chrX:1-100000", "GRCh38"), range("23:1-100000", "GRCh38"))
  expect_identical(range("24:1-10", "GRCh38"), range("chrY:1-10", "GRCh38"))
  expect_identical(range("1:1-1e6", "GRCh38"), range("1:1-1,000,000", "GRCh38"))
})

test_that("reads accept the extended region syntax", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(30L)
  input$chromosome <- rep(c("1", "23"), each = 15L)
  input$variant_id <- paste(input$chromosome, input$base_pair_location,
                            input$other_allele, input$effect_allele, sep = ":")
  store <- compress_sumstats(input, tempfile("region-syntax-"), overwrite = TRUE)
  columns <- c("chromosome", "base_pair_location")
  a <- read_sumstats(store, region = "1:1-1e6", columns = columns)
  b <- read_sumstats(store, region = "chr1:1-1,000,000", columns = columns)
  expect_identical(a, b)
  expect_identical(nrow(a), 15L)
  x <- read_sumstats(store, region = "23:1-1e6", columns = columns)
  expect_identical(nrow(x), 15L)
  expect_true(all(x$chromosome == "X"))
})
