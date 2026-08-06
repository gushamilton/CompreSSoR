test_that("report mode drops invalid indels and records structural rejection", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")

  input <- make_fixture(4L)
  input$reference_allele[1L] <- "AT"
  input$alternate_allele[1L] <- "G"
  input$effect_allele[1L] <- "G"
  input$other_allele[1L] <- "AT"
  path <- file.path(tempdir(), "report-invalid-indel.cpr")

  store <- compress_sumstats(input, path, row_policy = "report",
                             overwrite = TRUE)

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

  store <- compress_sumstats(input, path, row_policy = "report",
                             overwrite = TRUE)

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

  store <- compress_sumstats(input, path, backend = "parquet", profile = "exact",
                             row_policy = "report", overwrite = TRUE)

  expect_identical(store$manifest$n_rows, 3L)
  qc <- store$manifest$preparation$preparation$structural_qc
  expect_identical(qc$accepted_rows, 3L)
  expect_identical(qc$rejected_rows, 1L)
  expect_identical(unname(qc$rejection_counts[["indel"]]), 1L)
  expect_identical(nrow(read_sumstats(store)), 3L)
})
