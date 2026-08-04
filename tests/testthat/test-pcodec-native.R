withr::local_options(CompreSSoR.native_pcodec = TRUE)

test_that("native Pcodec is available and round trips integer streams", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  for (spec in list(
    list(dtype = "u8", values = c(0L, 1L, 7L, 127L, 255L)),
    list(dtype = "u16", values = c(0L, 1L, 255L, 4096L, 65535L)),
    list(dtype = "u32", values = c(0, 1, 2147483648, 4000000000, 4294967295))
  )) {
    encoded <- CompreSSoR:::pcodec_native_compress(spec$values, spec$dtype)
    decoded <- CompreSSoR:::pcodec_native_decompress(encoded, length(spec$values), spec$dtype)
    expect_identical(decoded, spec$values)
  }
})

test_that("native 0.4 stores are the default and support full, regional, key, and row reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(2000L)
  path <- tempfile("pcodec-native-")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  expect_equal(store$manifest$format_version, "0.4.0-pcodec-native")
  expect_equal(store$manifest$codec$name, "pcodec_native_standalone_z9_eaf8_se6")
  expect_false(file.exists(file.path(path, "variants.parquet")))
  expect_true(validate_compressor(store, full = TRUE)$valid)

  full <- read_sumstats(store)
  expect_equal(nrow(full), nrow(input))
  expect_false("row" %in% names(full))
  expect_false("variant_id" %in% names(full))
  expect_equal(full[c("chromosome", "base_pair_location", "effect_allele", "other_allele")],
               input[c("chromosome", "base_pair_location", "effect_allele", "other_allele")])
  expect_equal(full$effect_allele_frequency, input$effect_allele_frequency, tolerance = 0.006)
  expect_lt(max(abs(full$beta - input$beta), na.rm = TRUE), 0.02)
  expect_lt(max(abs(full$standard_error - input$standard_error), na.rm = TRUE), 0.01)

  regional <- read_sumstats(store, region = "chr1:100100-100150",
                            columns = c("chromosome", "base_pair_location", "z", "beta",
                                         "standard_error", "effect_allele_frequency"))
  expect_true(all(regional$base_pair_location >= 100100L &
                  regional$base_pair_location <= 100150L))
  key <- compressor_variant_key(input$chromosome[17L], input$base_pair_location[17L],
                                input$other_allele[17L], input$effect_allele[17L])
  by_key <- read_sumstats(store, variants = key,
                          columns = c("chromosome", "base_pair_location", "effect_allele",
                                      "other_allele", "z"))
  expect_equal(nrow(by_key), 1L)
  expect_equal(compressor_variant_key(by_key$chromosome, by_key$base_pair_location,
                                      by_key$other_allele, by_key$effect_allele), key)
  by_row <- read_sumstats(store, variants = 16L, columns = c("z", "standard_error"))
  expect_equal(nrow(by_row), 1L)
  expect_false("row" %in% names(by_row))
  projected <- read_sumstats(store, columns = "z")
  expect_identical(names(projected), "z")
  expect_equal(nrow(projected), nrow(input))
})

test_that("native Pcodec preserves exceptional values and batched reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(100L)
  input$beta[c(1, 2)] <- c(100000, -100)
  input$standard_error[c(1, 2)] <- c(1000, 1e-6)
  input$effect_allele_frequency[3] <- NA_real_
  input$z <- input$beta / input$standard_error
  path1 <- tempfile("pcodec-native-batch-a-")
  path2 <- tempfile("pcodec-native-batch-b-")
  compress_sumstats(input, path1, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  compress_sumstats(input, path2, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  expect_true(validate_compressor(path1, full = TRUE)$valid)
  keys <- compressor_variant_key(input$chromosome, input$base_pair_location,
                                 input$other_allele, input$effect_allele)
  observed <- read_sumstats_batch(c(path1, path2), list(keys[c(1, 3, 100)], keys[2]),
                                  columns = c("z", "standard_error", "effect_allele_frequency"),
                                  threads = 2L)
  expect_length(observed, 2L)
  expect_equal(nrow(observed[[1L]]), 3L)
  expect_equal(nrow(observed[[2L]]), 1L)
  expect_equal(observed[[1L]]$z[1], input$z[1], tolerance = 1e-5)
  expect_equal(observed[[1L]]$effect_allele_frequency[2], input$effect_allele_frequency[3])
})
