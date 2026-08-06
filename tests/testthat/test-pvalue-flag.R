test_that("native stores write and read the standard p-value flag domain", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(256L)
  input$p_value <- c(1e-9, 5e-8, 5.0001e-8,
                     rep(1e-3, nrow(input) - 3L))
  path <- tempfile("pvalue-flag-")
  store <- compress_sumstats(input, path, overwrite = TRUE, threads = 1L)

  domain <- store$manifest$domains$pvalue_flag
  expect_true(is.list(domain))
  expect_equal(domain$threshold, 5e-8)
  expect_identical(domain$operator, "<=")
  expect_identical(domain$dtype, "uint8")
  expect_identical(domain$rows, nrow(input))
  expect_true(file.exists(file.path(path, "pvalue_flag.pco")))

  identity <- CompreSSoR:::pcodec_native_identity(input)
  row_order <- order(identity$global_position, identity$substitution, method = "radix")
  expected <- which(input$p_value[row_order] <= 5e-8) - 1L
  expect_identical(read_pvalue_flag(store), as.integer(expected))
  expected_logical <- rep(FALSE, nrow(input))
  expected_logical[expected + 1L] <- TRUE
  expect_identical(read_pvalue_flag(store, as = "logical"), expected_logical)

  mr <- read_sumstats(
    store, variants = expected,
    columns = c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error",
                "effect_allele_frequency")
  )
  expect_equal(nrow(mr), length(expected))
  expect_equal(mr$base_pair_location, input$base_pair_location[row_order][expected + 1L])
  expect_true(validate_compressor(store, full = TRUE)$valid)
})

test_that("p-value flag threshold is configurable and remains separate from selection", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(64L)
  input$p_value <- rep(1e-3, nrow(input))
  input$p_value[c(2L, 7L)] <- 1e-4
  path <- tempfile("pvalue-flag-threshold-")
  store <- compress_sumstats(
    input, path, overwrite = TRUE, threads = 1L,
    pvalue_flag_threshold = 1e-4
  )
  expect_identical(read_pvalue_flag(store), as.integer(c(1L, 6L)))
  expect_equal(store$manifest$domains$pvalue_flag$threshold, 1e-4)
})

test_that("p-value flag derives p from Z when no source p-value is supplied", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(64L)
  input$p_value <- NULL
  path <- tempfile("pvalue-flag-derived-")
  store <- compress_sumstats(
    input, path, overwrite = TRUE, threads = 1L,
    pvalue_flag_threshold = 0.05
  )
  identity <- CompreSSoR:::pcodec_native_identity(input)
  row_order <- order(identity$global_position, identity$substitution, method = "radix")
  expected_p <- 2 * stats::pnorm(-abs(input$beta / input$standard_error))
  expected <- which(expected_p[row_order] <= 0.05) - 1L
  expect_identical(read_pvalue_flag(store), as.integer(expected))
  expect_identical(store$manifest$domains$pvalue_flag$source, "derived_from_z")
})

test_that("native p-value flag can be disabled", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  path <- tempfile("pvalue-flag-disabled-")
  store <- compress_sumstats(
    make_fixture(64L), path, overwrite = TRUE, threads = 1L,
    pvalue_flag = FALSE
  )
  expect_false(file.exists(file.path(path, "pvalue_flag.pco")))
  expect_null(store$manifest$domains$pvalue_flag)
  expect_error(read_pvalue_flag(store), "no p-value flag domain")
})

test_that("p-value flag arguments are validated", {
  expect_error(
    compress_sumstats(make_fixture(8L), tempfile("pvalue-flag-invalid-"),
                      pvalue_flag_threshold = 2),
    "pvalue_flag_threshold"
  )
  expect_error(
    compress_sumstats(make_fixture(8L), tempfile("pvalue-flag-invalid-"),
                      pvalue_flag = NA),
    "pvalue_flag"
  )
})
