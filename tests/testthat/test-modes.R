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

test_that("chromosome-parallel harmonisation matches serial output", {
  input <- make_edge_case_sumstats()
  reference <- make_edge_case_reference()
  serial <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 1L)
  parallel <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 2L)
  cols <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "z", "effect_allele_frequency", "variant_id")
  expect_equal(serial[cols], parallel[cols])
})
