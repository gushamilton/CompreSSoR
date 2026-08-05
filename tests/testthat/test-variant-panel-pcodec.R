make_variant_panel_fixture <- function() {
  data.frame(
    variant_id = c("2:100:G:A", "1:200:C:T", "1:100:A:G", "X:300:T:C"),
    hm3 = c(1L, 0L, 1L, 0L),
    stringsAsFactors = FALSE
  )
}

test_that("the bundled core/HM3 panel resolves without environment variables", {
  path <- CompreSSoR:::bundled_variant_panel_path("core")
  skip_if(!nzchar(path) || !dir.exists(path), "bundled real panel is not present")
  core_chr1 <- CompreSSoR:::read_variant_set("core", chromosomes = "1")
  hm3_chr1 <- CompreSSoR:::read_variant_set("hm3", chromosomes = "1")
  read_info <- attr(core_chr1, "variant_panel_read")
  expect_true(isTRUE(read_info$selective))
  expect_lt(read_info$payload_bytes_read, read_info$payload_bytes_total)
  expect_lt(length(read_info$key_blocks), read_info$key_blocks_total)
  expect_equal(nrow(core_chr1), 492028L)
  expect_equal(sum(core_chr1$hm3), 99416L)
  expect_equal(nrow(hm3_chr1), 99416L)
  expect_true(all(hm3_chr1$hm3 == 1L))
  input <- data.frame(
    chromosome = "1", base_pair_location = 10583L,
    reference_allele = "G", alternate_allele = "A",
    effect_allele = "A", other_allele = "G", beta = 0.1,
    standard_error = 0.1, effect_allele_frequency = 0.2,
    stringsAsFactors = FALSE
  )
  output <- tempfile("bundled-core-selection-")
  store <- compress_sumstats(input, output, selection = "core", overwrite = TRUE)
  expect_equal(store$manifest$n_rows, 1L)
})

test_that("native variant panels round-trip exact identity and HM3 flags", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  path <- tempfile("variant-panel-")
  panel <- write_variant_panel(make_variant_panel_fixture(), path,
                               overwrite = TRUE)
  expect_s3_class(panel, "compressor_variant_panel")
  expect_true(validate_variant_panel(path, full = TRUE)$valid)
  got <- read_variant_panel(path)
  expect_equal(got$variant_id,
               c("1:100:A:G", "1:200:C:T", "2:100:G:A", "X:300:T:C"))
  expect_equal(got$hm3, c(1L, 0L, 1L, 0L))
  expect_equal(panel$manifest$rows, 4L)
  expect_equal(panel$manifest$core_rows, 4L)
  expect_equal(panel$manifest$hm3_rows, 2L)
  expect_true(isTRUE(panel$manifest$identity$exact))
  expect_identical(panel$manifest$tolerances$hm3, "lossless_binary_flag")
})

test_that("variant panel row, key, chromosome and HM3 access are compatible", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  path <- tempfile("variant-panel-access-")
  write_variant_panel(make_variant_panel_fixture(), path, overwrite = TRUE)
  expect_equal(read_variant_panel(path, rows = c(1, 4))$variant_id,
               c("1:100:A:G", "X:300:T:C"))
  expect_equal(read_variant_panel(path, keys = "2:100:G:A")$hm3, 1L)
  expect_equal(nrow(read_variant_panel(path, chromosomes = "chr1")), 2L)
  expect_equal(read_variant_panel(path, hm3_only = TRUE)$variant_id,
               c("1:100:A:G", "2:100:G:A"))
})

test_that("panel core universe and membership drive core selection", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  path <- tempfile("variant-panel-selection-")
  write_variant_panel(make_variant_panel_fixture(), path, overwrite = TRUE)
  panel <- read_variant_set(path)
  expect_equal(panel$hm3, c(1L, 0L, 1L, 0L))
  sumstats <- data.frame(
    chromosome = c("1", "1", "2", "3"),
    base_pair_location = c(100L, 999L, 100L, 100L),
    reference_allele = c("A", "A", "G", "A"),
    alternate_allele = c("G", "T", "A", "G"),
    effect_allele = c("G", "T", "A", "G"),
    other_allele = c("A", "A", "G", "A"),
    beta = c(.1, .2, .3, .4), standard_error = rep(.1, 4),
    effect_allele_frequency = rep(.2, 4), stringsAsFactors = FALSE
  )
  expect_equal(variant_set_membership(sumstats, panel), c(TRUE, FALSE, TRUE, FALSE))
  output <- tempfile("variant-panel-filtered-")
  store <- compress_sumstats(sumstats, output, selection = "core",
                             variant_set = path, overwrite = TRUE)
  expect_equal(store$manifest$n_rows, 2L)
})

test_that("variant panel corruption is detected by payload checksums", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  path <- tempfile("variant-panel-corrupt-")
  write_variant_panel(make_variant_panel_fixture(), path, overwrite = TRUE)
  stream <- file.path(path, "hm3.pco")
  bytes <- readBin(stream, raw(), n = file.info(stream)$size)
  bytes[[length(bytes)]] <- as.raw(bitwXor(as.integer(bytes[[length(bytes)]]), 1L))
  writeBin(bytes, stream)
  expect_error(validate_variant_panel(path, full = TRUE), "checksum mismatch")
})
