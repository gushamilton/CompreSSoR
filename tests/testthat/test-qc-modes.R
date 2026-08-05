test_that("compact and none modes write the same trusted fixture", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(24L)
  compact_path <- tempfile("compact-mode-")
  none_path <- tempfile("none-mode-")

  compact <- compress_sumstats(input, compact_path, overwrite = TRUE,
                               threads = 2L)
  none <- compress_sumstats(input, none_path, overwrite = TRUE,
                            qc = "none", threads = 2L)

  expect_equal(compact$manifest$n_rows, 24L)
  expect_equal(none$manifest$n_rows, 24L)
  expect_identical(compact$manifest$qc$mode, "compact")
  expect_identical(compact$manifest$qc$structural_qc, "compact")
  expect_identical(none$manifest$qc$mode, "none")
  expect_identical(none$manifest$qc$structural_qc, "bypassed")
  expect_identical(none$manifest$row_policy, "not_applied")
  expect_identical(none$manifest$preparation$preparation$qc$mode, "none")
  expect_true(isTRUE(validate_compressor(none, full = TRUE)$valid))
})

test_that("none mode filters through the numeric identity panel", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(4L)
  panel <- data.frame(
    variant_id = input$variant_id[c(1L, 3L)],
    hm3 = c(1L, 0L),
    stringsAsFactors = FALSE
  )
  path <- tempfile("none-panel-")
  store <- compress_sumstats(input, path, qc = "none", selection = "core",
                             variant_set = panel, overwrite = TRUE)
  expect_identical(store$manifest$n_rows, 2L)
  expect_identical(store$manifest$selection$kept_rows, 2L)
  expect_identical(store$manifest$selection$dropped_rows, 2L)
})

test_that("none mode drops unsupported identity rows before canonicalization", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(7L)
  lengths <- CompreSSoR:::compressor_chromosome_lengths("GRCh38")
  input$chromosome[1L] <- "15"
  input$base_pair_location[1L] <- lengths[["15"]] + 1
  input$base_pair_location[2L] <- Inf
  input$base_pair_location[3L] <- input$base_pair_location[3L] + 0.5
  input$chromosome[4L] <- "MT"
  input$reference_allele[5L] <- "N"
  input$other_allele[5L] <- "N"
  input$reference_allele[6L] <- input$alternate_allele[6L]
  input$other_allele[6L] <- input$reference_allele[6L]

  store <- compress_sumstats(input, tempfile("none-identity-filter-"),
                             qc = "none", overwrite = TRUE)
  expect_equal(store$manifest$n_rows, 1L)

  preparation <- store$manifest$preparation$preparation
  safety <- preparation$identity_safety
  expect_equal(preparation$input_rows, 7L)
  expect_equal(preparation$rows_after_identity_safety, 1L)
  expect_equal(preparation$unsupported_identity_rows, 6L)
  expect_equal(preparation$dropped_unsupported_identity, 6L)
  expect_equal(safety$input_rows, 7L)
  expect_equal(safety$kept_rows, 1L)
  expect_equal(safety$dropped_rows, 6L)
  expect_equal(safety$counts$invalid_primary_chromosome, 1L)
  expect_equal(safety$counts$nonfinite_coordinate, 1L)
  expect_equal(safety$counts$noninteger_coordinate, 1L)
  expect_equal(safety$counts$coordinate_out_of_range, 1L)
  expect_equal(safety$counts$invalid_allele, 1L)
  expect_equal(safety$counts$same_alleles, 1L)

  expect_equal(store$manifest$selection$input_rows, 1L)
  expect_equal(store$manifest$selection$kept_rows, 1L)
  expect_equal(store$manifest$selection$dropped_rows, 0L)
  expect_equal(read_sumstats(store)$base_pair_location, input$base_pair_location[7L])
  expect_true(isTRUE(validate_compressor(store, full = TRUE)$valid))
})

test_that("HM3 selection honours flags in a combined custom panel", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(4L)
  panel <- data.frame(
    variant_id = input$variant_id,
    hm3 = c(1L, 0L, 1L, 0L),
    stringsAsFactors = FALSE
  )
  store <- compress_sumstats(
    input, tempfile("none-hm3-"), qc = "none", selection = "hm3",
    variant_set = panel, overwrite = TRUE
  )
  expect_identical(store$manifest$n_rows, 2L)
  expect_equal(store$manifest$selection$selection, "hm3")
  expect_equal(store$manifest$selection$kept_rows, 2L)
})

test_that("none mode requires exact canonical columns and refuses QC policies", {
  input <- make_fixture(2L)
  aliases <- input
  names(aliases)[names(aliases) == "reference_allele"] <- "REF"
  expect_error(
    compress_sumstats(aliases, tempfile("none-alias-"), qc = "none"),
    "already-canonical"
  )
  missing_beta <- input[setdiff(names(input), "beta")]
  expect_error(
    compress_sumstats(missing_beta, tempfile("none-missing-"), qc = "none"),
    "already-canonical.*beta"
  )
  expect_error(
    compress_sumstats(input, tempfile("none-strict-"), qc = "none", strict = TRUE),
    "cannot be combined"
  )
  expect_error(
    compress_sumstats(input, tempfile("none-error-policy-"), qc = "none",
                      row_policy = "error"),
    "cannot be combined"
  )

  malformed_beta <- input
  malformed_beta$beta[1L] <- NA_real_
  expect_error(
    compress_sumstats(malformed_beta, tempfile("none-malformed-beta-"), qc = "none"),
    "requires finite beta"
  )
  malformed_se <- input
  malformed_se$standard_error[1L] <- 0
  expect_error(
    compress_sumstats(malformed_se, tempfile("none-malformed-se-"), qc = "none"),
    "requires positive finite standard_error"
  )

  inconsistent_z <- input
  inconsistent_z$z <- inconsistent_z$beta / inconsistent_z$standard_error
  inconsistent_z$z[1L] <- inconsistent_z$z[1L] + 0.5
  expect_error(
    compress_sumstats(inconsistent_z, tempfile("none-inconsistent-z-"), qc = "none"),
    "supplied z is inconsistent"
  )
})

test_that("none mode retains native duplicate safety without structural QC", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(3L)
  input$reference_allele[1L] <- "AT"
  input$other_allele[1L] <- "AT"
  store <- compress_sumstats(input, tempfile("none-indel-"), qc = "none",
                             overwrite = TRUE)
  expect_equal(store$manifest$n_rows, 2L)
  expect_equal(store$manifest$preparation$preparation$identity_safety$dropped_rows, 1L)
  expect_equal(store$manifest$preparation$preparation$identity_safety$counts$invalid_allele, 1L)
  duplicate <- rbind(input[2L, ], input[2L, ])
  expect_error(
    compress_sumstats(duplicate, tempfile("none-duplicate-"), qc = "none"),
    "duplicate full REF/ALT identity keys"
  )
})

test_that("compact ingestion can omit synthetic variant IDs", {
  raw <- make_fixture(3L)
  imported <- CompreSSoR:::import_sumstats_impl(
    raw, project_columns = TRUE, core_only = TRUE, run_qc = FALSE,
    construct_variant_id = FALSE
  )
  expect_false("variant_id" %in% names(imported))
  expect_false("rsid" %in% names(imported))
})
