test_that("q8 cache builds and reads regional values", {
  skip_if_not_installed("arrow")
  input <- make_fixture(1200L)
  path <- file.path(tempdir(), "cache-store.cpr")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             profile = "standard", backend = "parquet", overwrite = TRUE)
  cache_path <- build_cache(store, overwrite = TRUE, block_rows = 256L)
  expect_true(dir.exists(cache_path))
  got <- read_sumstats(store, region = "chr1:100100-100250", use_cache = TRUE)
  expect_equal(nrow(got), 151L)
  expect_true(all(is.finite(got$standard_error)))
  expect_true(all(c("beta", "standard_error", "effect_allele_frequency", "p_value") %in% names(got)))
  expect_false("p_value" %in% names(arrow::read_parquet(file.path(path, "values.parquet"))))
})
