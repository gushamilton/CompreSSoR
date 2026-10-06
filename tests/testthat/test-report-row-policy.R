test_that("report mode drops invalid indels and records structural rejection", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(4L)
  input$reference_allele[1L] <- "AT"
  input$alternate_allele[1L] <- "G"
  input$effect_allele[1L] <- "G"
  input$other_allele[1L] <- "AT"
  path <- file.path(tempdir(), "report-invalid-indel.cpr")

  expect_warning(
    store <- compress_sumstats(input, path, row_policy = "report",
                               overwrite = TRUE),
    "dropped 1 of 4 input rows.*indel=1"
  )

  expect_identical(store$manifest$n_rows, 3L)
  qc <- store$manifest$preparation$preparation$structural_qc
  expect_identical(qc$accepted_rows, 3L)
  expect_identical(qc$rejected_rows, 1L)
  expect_identical(unname(qc$rejection_counts[["indel"]]), 1L)

  got <- read_sumstats(store, columns = c("chromosome", "base_pair_location",
                                           "reference_allele", "alternate_allele",
                                           "effect_allele", "other_allele"))
  expect_identical(nrow(got), 3L)
  expect_false(any(got$base_pair_location == input$base_pair_location[1L]))
})

test_that("error mode remains fail-closed for invalid indels", {
  input <- make_fixture(4L)
  input$reference_allele[1L] <- "AT"
  input$alternate_allele[1L] <- "G"
  input$effect_allele[1L] <- "G"
  input$other_allele[1L] <- "AT"

  expect_error(
    compress_sumstats(input, file.path(tempdir(), "error-invalid-indel.cpr"),
                      row_policy = "error", overwrite = TRUE),
    "explicit REF and ALT"
  )
})

test_that("report mode drops inconsistent orientation without flipping alleles", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(4L)
  input$effect_allele[1L] <- input$reference_allele[1L]
  input$other_allele[1L] <- input$alternate_allele[1L]
  path <- file.path(tempdir(), "report-inconsistent-orientation.cpr")

  expect_warning(
    store <- compress_sumstats(input, path, row_policy = "report",
                               overwrite = TRUE),
    "orientation_mismatch=1"
  )

  qc <- store$manifest$preparation$preparation$structural_qc
  expect_identical(qc$rejected_rows, 1L)
  expect_identical(unname(qc$rejection_counts[["orientation_mismatch"]]), 1L)
  got <- read_sumstats(store, columns = c("base_pair_location", "reference_allele",
                                           "alternate_allele", "effect_allele",
                                           "other_allele"))
  expect_false(any(got$base_pair_location == input$base_pair_location[1L]))

  expect_error(
    compress_sumstats(input, file.path(tempdir(), "error-inconsistent-orientation.cpr"),
                      row_policy = "error", overwrite = TRUE),
    "inconsistent with explicit REF/ALT"
  )
})

test_that("report mode still requires explicit REF/ALT and never silently flips", {
  input <- make_fixture(2L)
  missing_identity <- input[c("chromosome", "base_pair_location", "effect_allele",
                              "other_allele", "beta", "standard_error")]
  expect_error(
    compress_sumstats(missing_identity,
                      file.path(tempdir(), "report-missing-ref-alt.cpr"),
                      row_policy = "report"),
    "requires explicit REF and ALT"
  )

  flipped <- make_fixture(2L)
  flipped$effect_allele <- flipped$reference_allele
  flipped$other_allele <- flipped$alternate_allele
  expect_error(
    compress_sumstats(flipped,
                      file.path(tempdir(), "report-flipped.cpr"),
                      row_policy = "error"),
    "inconsistent with explicit REF/ALT"
  )
})

test_that("report mode fails clearly when every row is structurally rejected", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(3L)
  input$reference_allele[1L] <- "AT"
  input$alternate_allele[1L] <- "G"
  input$effect_allele[1L] <- "G"
  input$other_allele[1L] <- "AT"
  input$effect_allele[2L] <- input$reference_allele[2L]
  input$other_allele[2L] <- input$alternate_allele[2L]
  input$reference_allele[3L] <- "N"
  input$other_allele[3L] <- "N"
  path <- file.path(tempdir(), "report-all-structurally-invalid.cpr")

  expect_error(
    compress_sumstats(input, path, row_policy = "report", overwrite = TRUE),
    "structural QC rejected 3 row\\(s\\)"
  )
  expect_false(dir.exists(path))
})

test_that("report-mode structural filtering is shared by the Parquet backend", {
  skip_if_not_installed("arrow")

  input <- make_fixture(4L)
  input$reference_allele[1L] <- "AT"
  input$alternate_allele[1L] <- "G"
  input$effect_allele[1L] <- "G"
  input$other_allele[1L] <- "AT"
  path <- file.path(tempdir(), "report-invalid-indel.parquet.cpr")

  expect_warning(
    store <- compress_sumstats(input, path, backend = "parquet", profile = "exact",
                               row_policy = "report", overwrite = TRUE),
    "indel=1"
  )

  expect_identical(store$manifest$n_rows, 3L)
  qc <- store$manifest$preparation$preparation$structural_qc
  expect_identical(qc$accepted_rows, 3L)
  expect_identical(qc$rejected_rows, 1L)
  expect_identical(unname(qc$rejection_counts[["indel"]]), 1L)
  expect_identical(nrow(read_sumstats(store)), 3L)
})

test_that("report mode warns once with the total and per-reason counts", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(1000L)
  input$standard_error[2:1000] <- 0
  path <- file.path(tempdir(), "report-warns-bad-se.cpr")
  warnings <- character()
  store <- withCallingHandlers(
    compress_sumstats(input, path, row_policy = "report", overwrite = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(store$manifest$n_rows, 1L)
  expect_length(warnings, 1L)
  expect_match(warnings, "dropped 999 of 1,000 input rows", fixed = TRUE)
  expect_match(warnings, "invalid_standard_error=999", fixed = TRUE)
  expect_match(warnings, "row_policy = \"error\"", fixed = TRUE)

  # Nothing dropped: no warning. row_policy = "error" still stops instead.
  expect_no_warning(compress_sumstats(make_fixture(50L), path, overwrite = TRUE))
  expect_error(
    compress_sumstats(input, path, row_policy = "error", overwrite = TRUE),
    "structural QC rejected"
  )
})

test_that("the drop warning names duplicate keys with differing values", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(20L)
  same <- input[5L, ]
  differing <- input[9L, ]
  differing$beta <- differing$beta + 0.5
  input <- rbind(input, same, differing)
  path <- file.path(tempdir(), "report-warns-duplicates.cpr")
  expect_warning(
    store <- compress_sumstats(input, path, row_policy = "report",
                               overwrite = TRUE),
    "duplicate_variant=2.*1 duplicated key\\(s\\) had differing values"
  )
  expect_identical(store$manifest$n_rows, 20L)
  expect_identical(store$manifest$dropped_rows$input_rows, 22L)
})

test_that("finite effects that overflow on read are rejected at write time", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(10L)
  input$beta[2L] <- 1e308
  input$standard_error[2L] <- 1
  input$beta[3L] <- -1e308
  input$standard_error[3L] <- 1
  input$beta[4L] <- 1e300
  input$standard_error[4L] <- 1e-300
  path <- file.path(tempdir(), "report-overflow.cpr")
  expect_warning(
    store <- compress_sumstats(input, path, row_policy = "report",
                               overwrite = TRUE),
    "dropped 3 of 10 input rows.*non_finite_effect=2"
  )
  expect_identical(store$manifest$n_rows, 7L)
  got <- read_sumstats(store, columns = c("beta", "z", "standard_error"))
  expect_true(all(is.finite(got$beta)))
  expect_true(all(is.finite(got$z)))

  expect_error(
    compress_sumstats(input, path, row_policy = "error", overwrite = TRUE),
    "non_finite_effect=2"
  )
  # Large but representable effects are kept and read back finite.
  big <- make_fixture(10L)
  big$beta[2L] <- 1e30
  big$standard_error[2L] <- 1
  expect_no_warning(store <- compress_sumstats(big, path, overwrite = TRUE))
  expect_true(all(is.finite(read_sumstats(store, columns = "beta")$beta)))
})
