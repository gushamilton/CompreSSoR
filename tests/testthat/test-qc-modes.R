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
})

test_that("none mode retains codec identity safety without structural QC", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(3L)
  input$reference_allele[1L] <- "AT"
  input$other_allele[1L] <- "AT"
  expect_error(
    compress_sumstats(input, tempfile("none-indel-"), qc = "none"),
    "REF and ALT"
  )
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
