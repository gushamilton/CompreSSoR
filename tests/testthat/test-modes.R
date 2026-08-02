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
    variant_id = paste0("v", 1:5),
    stringsAsFactors = FALSE
  )
}

test_that("default QC mode preserves unresolved variants and records status", {
  input <- make_edge_case_sumstats()
  reference <- make_edge_case_reference()
  got <- harmonise_sumstats(input, reference, mode = "qc")
  expect_equal(nrow(got), nrow(input))
  status_counts <- table(got$harmonisation_status)
  expect_equal(names(status_counts), c("aligned", "ambiguous", "incompatible", "unmatched"))
  expect_equal(as.integer(status_counts), c(4L, 1L, 1L, 1L))
  expect_equal(sum(got$input_duplicate), 1L)
  expect_true(got$harmonisation_flip[2])
  expect_equal(got$beta[2], -input$beta[2])
  expect_equal(got$effect_allele[3], "G")
  expect_equal(got$other_allele[3], "A")
  expect_equal(got$beta[3], input$beta[3])
  expect_equal(attr(got, "alignment_stats")$output_rows, nrow(input))
  expect_error(harmonise_sumstats(input[-7, , drop = FALSE], reference, mode = "qc", strict = TRUE),
               "reference alignment failed")
})

test_that("conversion, core and HM3 modes are explicit", {
  input <- make_edge_case_sumstats()
  reference <- make_edge_case_reference()
  converted <- harmonise_sumstats(input, reference, mode = "convert")
  expect_equal(nrow(converted), nrow(input))
  expect_true(all(converted$harmonisation_status == "unreferenced"))

  core <- harmonise_sumstats(input, reference, mode = "core",
                             variant_set = data.frame(variant_id = c("v1", "v2")))
  expect_equal(nrow(core), 3L)
  expect_true(all(core$variant_id %in% c("v1", "v2")))

  hm3_path <- tempfile(fileext = ".bim")
  writeLines(c(
    "1\tv1\t0\t100\tC\tA",
    "1\tv2\t0\t200\tC\tA"
  ), hm3_path)
  hm3 <- harmonise_sumstats(input, reference, mode = "hm3", variant_set = hm3_path)
  expect_equal(nrow(hm3), 3L)
  expect_true(all(hm3$variant_id %in% c("v1", "v2")))
})

test_that("chromosome-parallel harmonisation matches serial output", {
  input <- make_edge_case_sumstats()
  reference <- make_edge_case_reference()
  serial <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 1L)
  parallel <- harmonise_sumstats(input, reference, mode = "qc", chrom_threads = 2L)
  cols <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "effect_allele_frequency", "variant_id", "harmonisation_status",
            "harmonisation_flip", "input_duplicate")
  expect_equal(serial[cols], parallel[cols])
})
