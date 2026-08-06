core_plus_test_data <- function(positions, z, p_value = NULL) {
  se <- rep(0.1, length(positions))
  out <- data.frame(
    chromosome = rep("1", length(positions)),
    base_pair_location = as.integer(positions),
    reference_allele = rep("A", length(positions)),
    alternate_allele = rep("G", length(positions)),
    effect_allele = rep("G", length(positions)),
    other_allele = rep("A", length(positions)),
    beta = as.numeric(z) * se,
    standard_error = se,
    z = as.numeric(z),
    effect_allele_frequency = rep(0.2, length(positions)),
    stringsAsFactors = FALSE
  )
  if (!is.null(p_value)) out$p_value <- as.numeric(p_value)
  out
}

test_that("core-plus uses an inclusive threshold and inclusive 10-kb window", {
  threshold <- 1e-5
  z_boundary <- stats::qnorm(threshold / 2, lower.tail = FALSE)
  threshold_at_boundary <- 2 * stats::pnorm(-abs(z_boundary))
  positions <- c(89999L, 90000L, 100000L, 110000L, 110001L, 300000L)
  z <- c(0, 0, z_boundary, 0, 0, 0)
  panel <- data.frame(variant_id = "1:300000:A:G", stringsAsFactors = FALSE)

  selected <- CompreSSoR:::select_core_plus_variants(
    core_plus_test_data(positions, z), variant_set = panel,
    pvalue_threshold = threshold_at_boundary, region_padding = 10000L
  )

  expect_equal(selected$keep, c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE))
  expect_equal(selected$regions$chromosome, "1")
  expect_equal(selected$regions$start, 90000L)
  expect_equal(selected$regions$end, 110000L)
  expect_equal(selected$metadata$seed_snps, 1L)
  expect_equal(selected$metadata$threshold_operator, "<=")
  expect_equal(selected$metadata$padding_bp, 10000L)
  expect_equal(selected$metadata$window_bp_each_side, 10000L)
  expect_equal(selected$metadata$window_boundary, "inclusive")
  expect_equal(selected$metadata$union, "core_or_pvalue_regions")
})

test_that("overlapping seed windows merge deterministically and preserve the core union", {
  threshold <- 1e-5
  z_signal <- stats::qnorm(threshold / 4, lower.tail = FALSE)
  positions <- c(89999L, 90000L, 100000L, 110000L, 115000L, 125000L,
                125001L, 300000L)
  z <- c(0, 0, z_signal, 0, z_signal, 0, 0, 0)
  panel <- data.frame(variant_id = "1:300000:A:G", stringsAsFactors = FALSE)

  selected <- CompreSSoR:::select_core_plus_variants(
    core_plus_test_data(positions, z), variant_set = panel,
    pvalue_threshold = threshold, region_padding = 10000L
  )

  expect_equal(nrow(selected$regions), 1L)
  expect_equal(selected$regions$start, 90000L)
  expect_equal(selected$regions$end, 125000L)
  expect_equal(selected$regions$seed_snps, 2L)
  expect_equal(selected$keep, c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE))
})

test_that("missing or non-finite Z values cannot create p-value seeds", {
  positions <- c(100000L, 200000L, 300000L, 400000L)
  z <- c(NA_real_, Inf, -Inf, NaN)
  selected <- CompreSSoR:::pvalue_region_selection(
    core_plus_test_data(positions, z), pvalue_threshold = 1e-5,
    region_padding = 10000L
  )

  expect_false(any(selected$keep))
  expect_equal(selected$metadata$seed_snps, 0L)
  expect_equal(nrow(selected$regions), 0L)
})

test_that("selection derives p from beta/SE and does not substitute an input p column", {
  threshold <- 1e-5
  z_signal <- stats::qnorm(threshold / 4, lower.tail = FALSE)
  raw <- core_plus_test_data(c(100000L, 200000L), c(z_signal, 0),
                             p_value = c(1, NA_real_))
  raw$z <- NULL
  resolved <- CompreSSoR:::normalise_sumstats_columns(raw, parse_policy = "error")
  expect_equal(resolved$z[[1L]], z_signal)

  selected <- CompreSSoR:::pvalue_region_selection(
    resolved, pvalue_threshold = threshold, region_padding = 10000L
  )
  expect_true(selected$keep[[1L]])
  expect_false(selected$keep[[2L]])
  expect_equal(selected$metadata$p_value_source, "derived_from_z")
  expect_equal(selected$metadata$threshold_source, "pre_encoding_prepared")
})

test_that("a Z-only table is rejected by the strict core contract", {
  input <- core_plus_test_data(100000L, 5)
  input$beta <- NULL
  input$standard_error <- NULL
  expect_error(
    CompreSSoR:::normalise_sumstats_columns(input, parse_policy = "error"),
    "Missing required summary-statistics columns: beta, standard_error"
  )
})

test_that("core-plus records the threshold/window contract in the manifest", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  threshold <- 1e-5
  z_signal <- stats::qnorm(threshold / 4, lower.tail = FALSE)
  input <- core_plus_test_data(c(100000L, 110000L, 300000L),
                               c(z_signal, 0, 0))
  panel <- data.frame(variant_id = "1:300000:A:G", stringsAsFactors = FALSE)
  path <- tempfile("core-plus-manifest-")
  store <- compress_sumstats(
    input, path, selection = "core_plus", variant_set = panel,
    pvalue_threshold = threshold, region_padding = 10000L,
    overwrite = TRUE
  )

  selection <- store$manifest$selection
  expect_equal(selection$selection, "core_plus")
  expect_equal(selection$pvalue_threshold, threshold)
  expect_equal(selection$padding_bp, 10000L)
  expect_equal(selection$window_bp_each_side, 10000L)
  expect_equal(selection$window_boundary, "inclusive")
  expect_equal(selection$threshold_operator, "<=")
  expect_equal(selection$union, "core_or_pvalue_regions")
  expect_equal(selection$p_value_source, "derived_from_z")

  sidecar <- jsonlite::read_json(
    file.path(path, selection$file), simplifyVector = TRUE
  )
  expect_equal(sidecar$pvalue_threshold, threshold)
  expect_equal(sidecar$window_bp_each_side, 10000L)
  expect_equal(sidecar$window_boundary, "inclusive")
  expect_equal(sidecar$threshold_operator, "<=")
  expect_equal(sidecar$union, "core_or_pvalue_regions")
})
