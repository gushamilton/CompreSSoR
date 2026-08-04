make_edge_case_sumstats <- function() {
  data.frame(
    chromosome = c("1", "1", "1", "1", "1", "2", "1"),
    base_pair_location = c(100L, 200L, 300L, 400L, 500L, 600L, 100L),
    effect_allele = c("A", "C", "C", "A", "T", "A", "A"),
    other_allele = c("C", "A", "T", "T", "G", "G", "C"),
    beta = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.1),
    standard_error = rep(0.1, 7),
    effect_allele_frequency = c(0.2, 0.3, 0.4, 0.5, 0.2, 0.3, 0.2),
    p_value = rep(0.01, 7),
    variant_id = c("v1", "v2", "v3", "v4", "v5", "v6", "v1"),
    stringsAsFactors = FALSE
  )
}

make_edge_case_reference <- function() {
  data.frame(
    chromosome = rep("1", 5),
    base_pair_location = c(100L, 200L, 300L, 400L, 500L),
    effect_allele = c("A", "A", "G", "A", "G"),
    other_allele = c("C", "C", "A", "T", "A"),
    effect_allele_frequency = c(0.2, 0.7, 0.4, 0.5, 0.2),
    rsid = paste0("rs", 1:5),
    stringsAsFactors = FALSE
  )
}

test_that("QC uses the canonical allele key and drops unresolved rows", {
  input <- make_edge_case_sumstats()
  reference <- make_edge_case_reference()
  got <- harmonise_sumstats(input, reference, mode = "qc")
  expect_equal(nrow(got), 2L)
  expect_true(all(grepl("^1:[0-9]+:[ACGT]:[ACGT]$", got$variant_id)))
  expect_equal(attr(got, "alignment_stats")$dropped_unmatched, 2L)
  expect_equal(attr(got, "alignment_stats")$dropped_duplicates, 2L)
  expect_error(harmonise_sumstats(input, reference, mode = "qc", strict = TRUE),
               "reference alignment failed")
})

test_that("convert is an explicit no-reference escape hatch", {
  input <- make_edge_case_sumstats()
  converted <- harmonise_sumstats(input, reference = NULL, mode = "convert")
  expect_equal(nrow(converted), nrow(input))
  expect_true(all(converted$harmonisation_status == "unreferenced"))
  expect_error(harmonise_sumstats(input, reference = NULL, mode = "qc"), "reference is required")
})

test_that("variant-set filtering works in convert mode without harmonisation", {
  input <- make_edge_case_sumstats()
  panel <- input[c(1L, 2L), c("chromosome", "base_pair_location",
                               "other_allele", "effect_allele")]

  converted <- harmonise_sumstats(input, reference = NULL, mode = "convert",
                                   variant_set = panel)

  # The no-reference conversion path keeps both input rows that carry the
  # selected identity; it does not silently deduplicate them.
  expect_equal(nrow(converted), 3L)
  expect_equal(converted$base_pair_location, c(100L, 200L, 100L))
  expect_true(all(converted$harmonisation_status == "unreferenced"))
  expect_equal(attr(converted, "alignment_stats")$variant_set$name, "convert")
  expect_equal(attr(converted, "alignment_stats")$variant_set$kept_rows, 3L)
  expect_equal(attr(converted, "alignment_stats")$variant_set$dropped_rows, 4L)

  wrong_allele_panel <- panel[1L, , drop = FALSE]
  wrong_allele_panel$effect_allele <- "G"
  expect_error(
    harmonise_sumstats(input, reference = NULL, mode = "convert",
                       variant_set = wrong_allele_panel),
    "retained no variants"
  )
})

test_that("named common and tag panels resolve through environment variables", {
  input <- make_edge_case_sumstats()
  panel_path <- tempfile("common-variants-", fileext = ".tsv.gz")
  panel <- data.frame(variant_id = compressor_variant_key(
    input$chromosome[1:2], input$base_pair_location[1:2],
    input$other_allele[1:2], input$effect_allele[1:2]
  ), stringsAsFactors = FALSE
  )
  data.table::fwrite(panel, panel_path, sep = "\t")

  old_common <- Sys.getenv("COMPRESSOR_COMMON_VARIANTS", unset = NA_character_)
  old_tag <- Sys.getenv("COMPRESSOR_TAG_VARIANTS", unset = NA_character_)
  on.exit({
    if (is.na(old_common)) Sys.unsetenv("COMPRESSOR_COMMON_VARIANTS") else
      Sys.setenv(COMPRESSOR_COMMON_VARIANTS = old_common)
    if (is.na(old_tag)) Sys.unsetenv("COMPRESSOR_TAG_VARIANTS") else
      Sys.setenv(COMPRESSOR_TAG_VARIANTS = old_tag)
  }, add = TRUE)
  Sys.setenv(COMPRESSOR_COMMON_VARIANTS = panel_path,
             COMPRESSOR_TAG_VARIANTS = panel_path)

  common <- harmonise_sumstats(input, reference = NULL, mode = "convert",
                               variant_set = "common")
  tag <- harmonise_sumstats(input, reference = NULL, mode = "convert",
                            variant_set = "tag")

  expect_equal(common$variant_id, tag$variant_id)
  expect_equal(attr(common, "alignment_stats")$variant_set$name, "common")
  expect_equal(attr(tag, "alignment_stats")$variant_set$name, "tag")
  expect_equal(attr(common, "alignment_stats")$variant_set$rows, 2L)
})

test_that("chromosome-parallel harmonisation matches serial output", {
  input <- make_edge_case_sumstats()
  reference <- make_edge_case_reference()
  serial <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 1L)
  parallel <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 2L)
  cols <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "z", "effect_allele_frequency", "variant_id")
  expect_equal(serial[cols], parallel[cols])
})
