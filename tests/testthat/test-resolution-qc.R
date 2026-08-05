test_that("the resolution matrix is deterministic and beta plus SE is primary", {
  matrix <- CompreSSoR:::sumstats_resolution_matrix()
  expect_identical(matrix$field,
                   c("beta", "standard_error", "z", "odds_ratio", "p_value"))
  expect_match(matrix$primary_route[[1L]], "explicit beta")
  expect_match(matrix$strict_fallback[[3L]], "beta / standard_error")
  expect_match(matrix$strict_fallback[[5L]], "disabled")

  input <- data.frame(
    CHR = "1", GENPOS = 101L, REF = "A", ALT = "G",
    ALLELE1 = "G", ALLELE0 = "A", BETA = 0.2, SE = 0.1, Z = 2,
    stringsAsFactors = FALSE
  )
  imported <- import_sumstats(input, row_policy = "error")
  expect_equal(imported$z, 2)
  expect_identical(attr(imported, "resolution_provenance")$z_source,
                   "supplied_and_checked_when_beta_and_standard_error_are_finite")
  expect_identical(attr(imported, "resolution_provenance")$p_to_se$enabled, FALSE)
})

test_that("conflicting Z is rejected and OR conversion stays at the boundary", {
  conflicting <- data.frame(
    chr = "1", pos = 101L, REF = "A", ALT = "G", EA = "G", OA = "A",
    beta = 0.2, SE = 0.1, Z = 3, stringsAsFactors = FALSE
  )
  expect_error(import_sumstats(conflicting, row_policy = "error"),
               "beta, z and standard_error are inconsistent")

  ambiguous <- conflicting
  ambiguous$Z <- 2
  ambiguous$OR <- 1.5
  expect_error(import_sumstats(ambiguous, row_policy = "error"),
               "beta and odds_ratio are inconsistent")

  odds_ratio <- conflicting[c("chr", "pos", "REF", "ALT", "EA", "OA", "SE")]
  odds_ratio$SE <- 0.1
  odds_ratio$OR <- 1.2
  imported <- import_sumstats(odds_ratio, row_policy = "error")
  expect_equal(imported$beta, log(1.2))
  expect_identical(attr(imported, "resolution_provenance")$beta_route,
                   "positive_OR_to_log_OR_boundary_conversion")
})

test_that("direct compression accepts OR plus SE and records resolution provenance", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- data.frame(
    CHR = "1", GENPOS = 101L, REF = "A", ALT = "G",
    ALLELE1 = "G", ALLELE0 = "A", OR = 1.2, SE = 0.1,
    stringsAsFactors = FALSE
  )
  path <- tempfile("or-se-native-")
  store <- compress_sumstats(input, path, overwrite = TRUE, row_policy = "error")
  got <- read_sumstats(store, columns = c("beta", "standard_error", "z"))

  expect_equal(got$beta, log(1.2), tolerance = 0.02)
  expect_equal(got$standard_error, 0.1, tolerance = 0.01)
  expect_identical(store$manifest$source$resolution$beta_route,
                   "positive_OR_to_log_OR_boundary_conversion")
  expect_true(all(c("OR", "SE") %in% unlist(store$manifest$source_columns_read)))
  expect_true(validate_compressor(store, full = TRUE)$valid)
})

test_that("p-value to SE inference is explicit, bounded, and conflict checked", {
  input <- data.frame(
    chr = "1", pos = 101L, REF = "A", ALT = "G", EA = "G", OA = "A",
    beta = 0.2, SE = NA_real_, P = 0.0455, stringsAsFactors = FALSE
  )
  disabled <- import_sumstats(input)
  expect_true(is.na(disabled$standard_error))
  expect_identical(attr(disabled, "resolution_provenance")$p_to_se$enabled, FALSE)
  expect_equal(unname(attr(disabled, "structural_qc_report")$rejections[["missing_statistics"]]), 1L)

  enabled <- import_sumstats(input, allow_p_to_se = TRUE, row_policy = "error")
  expect_true(is.finite(enabled$standard_error))
  expect_identical(attr(enabled, "resolution_provenance")$p_to_se$rows, 1L)
  expect_identical(attr(enabled, "resolution_provenance")$p_to_se$example_rows, 1L)

  conflict <- input
  conflict$SE <- 0.25
  expect_error(import_sumstats(conflict, allow_p_to_se = TRUE),
               "standard_error conflicts")
})

test_that("explicit p-to-SE resolves an absent SE column and remains opt-in", {
  input <- data.frame(
    chr = "1", pos = 101L, REF = "A", ALT = "G", EA = "G", OA = "A",
    beta = 0.2, P = 0.0455, stringsAsFactors = FALSE
  )
  expect_error(import_sumstats(input, row_policy = "error"),
               "Missing required summary-statistics columns: standard_error")

  enabled <- import_sumstats(input, allow_p_to_se = TRUE, row_policy = "error")
  expect_true(is.finite(enabled$standard_error))
  expect_identical(attr(enabled, "resolution_provenance")$standard_error_route,
                   "explicit_opt_in_p_value_to_standard_error_conversion")

  disabled_projection <- CompreSSoR:::read_sumstats_input(
    input, project_columns = TRUE, core_only = TRUE
  )
  enabled_projection <- CompreSSoR:::read_sumstats_input(
    input, project_columns = TRUE, core_only = TRUE, allow_p_to_se = TRUE
  )
  expect_false(any(toupper(names(disabled_projection)) == "P"))
  expect_true(any(toupper(names(enabled_projection)) == "P"))
})

test_that("compact QC retains aggregates but not row-level audit vectors", {
  input <- import_sumstats(make_fixture(5000L))
  input[2L, ] <- input[1L, ]
  input$beta[3L] <- NA_real_

  compact <- CompreSSoR:::apply_structural_qc(
    input, row_policy = "report", detail = "compact"
  )
  report <- compact$report
  expect_equal(nrow(compact$data), 4998L)
  expect_equal(report$rejected_rows, 3L)
  expect_equal(report$dropped_rows, 2L)
  expect_equal(unname(report$rejection_counts[["duplicate_variant"]]), 2L)
  expect_equal(unname(report$rejection_counts[["missing_statistics"]]), 1L)
  expect_true(all(c("rejection_counts", "examples", "accepted_rows", "dropped_rows") %in%
                  names(report)))
  expect_false(any(c("row_status", "canonical_key", "invalid_rows", "duplicate_rows",
                     "structurally_valid_rows", "internal") %in% names(report)))

  full <- CompreSSoR:::apply_structural_qc(
    input, row_policy = "report", detail = "full"
  )
  expect_true(all(c("row_status", "canonical_key", "invalid_rows") %in% names(full$report)))
  message(sprintf(
    "compact QC evidence: rows=%d columns=%d retained_audit_fields=%d rejected=%d",
    nrow(input), ncol(input),
    sum(c("row_status", "canonical_key", "invalid_rows") %in% names(report)),
    report$rejected_rows
  ))
})

test_that("import timings and canonical output fields are recorded", {
  input <- make_fixture(2000L)
  imported <- import_sumstats(input)
  timings <- attr(imported, "phase_timings")
  expect_identical(timings$unit, "seconds")
  expect_true(all(c("read", "normalise", "qc") %in% names(timings$phases)))
  expect_true(all(vapply(timings$phases, function(value) is.finite(value) && value >= 0,
                         logical(1))))
  expect_true(all(CompreSSoR:::required_sumstats_columns() %in% names(imported)))
  expect_true(all(c("z", "effect_allele_frequency", "p_value") %in% names(imported)))
})

test_that("native manifests record compact phase timings", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  store <- compress_sumstats(make_fixture(128L), tempfile("timed-native-"),
                             overwrite = TRUE, threads = 2L)
  timings <- store$manifest$timings
  expect_identical(timings$unit, "seconds")
  expect_true(all(c("read", "normalise", "qc", "identity_sort", "encode", "commit") %in%
                  names(timings$phases)))
  expect_true(all(vapply(timings$phases, function(value) is.finite(value) && value >= 0,
                         logical(1))))
  expect_true(validate_compressor(store, full = TRUE)$valid)
})
