test_that("exact conversion round-trips all supported columns", {
  skip_if_not_installed("arrow")
  input <- make_fixture(500L)
  path <- file.path(tempdir(), "roundtrip-exact.cpr")
  store <- compress_sumstats(input, path, profile = "exact", keep_extras = TRUE,
                             backend = "parquet", overwrite = TRUE)
  got <- decompress_sumstats(store)
  expect_true(validate_compressor(store)$valid)
  expect_equal(got[input_core <- c("chromosome", "base_pair_location", "effect_allele",
                                   "other_allele", "variant_id", "rsid", "beta",
                                   "standard_error", "effect_allele_frequency", "p_value",
                                   "annotation")], input[input_core])
})

test_that("standard conversion round-trips with measured codec error", {
  skip_if_not_installed("arrow")
  input <- make_fixture(5000L)
  path <- file.path(tempdir(), "roundtrip-standard.cpr")
  store <- compress_sumstats(input, path, profile = "standard",
                             backend = "parquet", overwrite = TRUE)
  got <- decompress_sumstats(store)
  expect_true(validate_compressor(store)$valid)
  expect_lt(max(abs(log10(got$p_value) - log10(input$p_value)), na.rm = TRUE), 0.08)
  expect_equal(got$effect_allele_frequency, input$effect_allele_frequency, tolerance = 0.006)
  expect_lt(max(abs(got$beta - input$beta), na.rm = TRUE), 0.01)
  expect_lt(max(abs(got$standard_error - input$standard_error), na.rm = TRUE), 0.001)
})

test_that("regional cache reads agree with the durable store and respect variants", {
  skip_if_not_installed("arrow")
  input <- make_fixture(1200L)
  path <- file.path(tempdir(), "roundtrip-cache.cpr")
  store <- compress_sumstats(input, path, profile = "standard",
                             backend = "parquet", overwrite = TRUE)
  build_cache(store, overwrite = TRUE, block_rows = 128L)
  direct <- read_sumstats(store, region = "chr1:100100-100250")
  cached <- read_sumstats(store, region = "chr1:100100-100250", use_cache = TRUE)
  cols <- c("row", "variant_id", "beta", "standard_error",
            "effect_allele_frequency", "p_value")
  expect_equal(nrow(cached), nrow(direct))
  expect_equal(cached$variant_id, direct$variant_id)
  expect_lt(max(abs(log10(cached$p_value) - log10(direct$p_value)), na.rm = TRUE), 0.2)
  expect_equal(cached$beta, direct$beta, tolerance = 0.02)
  selected <- read_sumstats(store, region = "chr1:100100-100250",
                            variants = input$variant_id[101:105], use_cache = TRUE)
  expect_equal(nrow(selected), 5L)
})

test_that("validation detects broken store files", {
  skip_if_not_installed("arrow")
  input <- make_fixture(20L)
  path <- file.path(tempdir(), "broken-store.cpr")
  store <- compress_sumstats(input, path, profile = "exact",
                             backend = "parquet", overwrite = TRUE)
  unlink(file.path(path, "values.parquet"))
  result <- validate_compressor(path)
  expect_false(result$valid)
  expect_true(any(grepl("values.parquet", result$errors)))
})
