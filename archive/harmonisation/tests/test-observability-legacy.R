test_that("harmonisation attaches compact deterministic phase timings", {
  reference <- data.frame(
    chromosome = c("1", "2"), base_pair_location = c(100L, 200L),
    reference_allele = c("A", "G"), alternate_allele = c("C", "T"),
    effect_allele_frequency = c(0.2, 0.3), stringsAsFactors = FALSE
  )
  input <- data.frame(
    chromosome = c("1", "2", "2"), base_pair_location = c(100L, 200L, 300L),
    effect_allele = c("C", "T", "A"), other_allele = c("A", "G", "C"),
    beta = c(0.2, 0.3, 0.1), standard_error = 0.1,
    effect_allele_frequency = c(0.2, 0.3, 0.4), stringsAsFactors = FALSE
  )

  got <- harmonise_sumstats(input, reference, chrom_threads = 1L)
  timings <- attr(got, "timings")
  expect_true(is.list(timings))
  expect_identical(timings$schema, "CompreSSoR.harmonisation.observability.v1")
  expect_identical(timings$clock, "proc.time.elapsed")
  expect_identical(timings$rows$input, 3L)
  expect_identical(timings$rows$output, 2L)
  expect_identical(timings$rows$dropped, 1L)
  expected <- c(
    "input_discovery_decompression_import", "structural_qc", "liftover",
    "reference_discovery_load_normalisation", "matching_alignment",
    "strand_handling", "final_qc_report_construction", "output_preparation"
  )
  expect_true(all(expected %in% names(timings$phases)))
  expect_true(is.finite(timings$elapsed_seconds) && timings$elapsed_seconds >= 0)
  expect_true(is.finite(timings$cpu_seconds) && timings$cpu_seconds >= 0)
  for (phase in timings$phases) {
    expect_true(is.finite(phase$elapsed_seconds) && phase$elapsed_seconds >= 0)
    expect_true(is.finite(phase$cpu_seconds) && phase$cpu_seconds >= 0)
  }
  expect_identical(attr(got, "audit")$timings, timings)
  expect_identical(attr(got, "alignment_stats")$timings, timings)

  skip_on_os("windows")
  parallel <- harmonise_sumstats(input, reference, chrom_threads = 2L)
  expect_equal(parallel[c("chromosome", "base_pair_location", "variant_id")],
               got[c("chromosome", "base_pair_location", "variant_id")])
  expect_true("chromosome_partitioning" %in% names(attr(parallel, "timings")$phases))
  expect_true("matching_alignment" %in% names(attr(parallel, "timings")$phases))
})

test_that("progress callbacks have a bounded metadata-only event schema", {
  reference <- data.frame(
    chromosome = c("1", "2"), base_pair_location = c(100L, 200L),
    reference_allele = c("A", "G"), alternate_allele = c("C", "T"),
    stringsAsFactors = FALSE
  )
  input <- data.frame(
    chromosome = c("1", "2"), base_pair_location = c(100L, 200L),
    effect_allele = c("C", "T"), other_allele = c("A", "G"),
    beta = c(0.2, 0.3), standard_error = 0.1,
    effect_allele_frequency = c(0.2, 0.3), stringsAsFactors = FALSE
  )
  events <- list()
  callback <- function(event) {
    events[[length(events) + 1L]] <<- event
  }
  harmonise_sumstats(
    input, reference,
    observability = list(level = "events", callback = callback, memory = FALSE)
  )
  expect_true(length(events) >= 2L)
  expect_true(any(vapply(events, function(x) identical(x$event, "finalized"), logical(1L))))
  forbidden <- c("chromosome", "base_pair_location", "variant_id", "effect_allele",
                 "other_allele", "beta", "standard_error", "effect_allele_frequency")
  for (event in events) {
    expect_identical(event$schema, "CompreSSoR.harmonisation.observability.v1")
    expect_true(event$event %in% c("phase_start", "phase_end", "progress", "finalized"))
    expect_true(all(c("elapsed_seconds", "cpu_seconds", "rows_in", "rows_out",
                      "rows_dropped", "workers") %in% names(event)))
    expect_true(!any(names(event) %in% forbidden))
  }
})

test_that("error finalizes observability without returning partial data", {
  reference <- data.frame(
    chromosome = "1", base_pair_location = 100L,
    reference_allele = "A", alternate_allele = "C",
    stringsAsFactors = FALSE
  )
  input <- data.frame(
    chromosome = c("1", "1"), base_pair_location = c(100L, 999L),
    effect_allele = c("C", "G"), other_allele = c("A", "A"),
    beta = c(0.2, 0.3), standard_error = 0.1,
    effect_allele_frequency = c(0.2, 0.3), stringsAsFactors = FALSE
  )
  events <- list()
  callback <- function(event) events[[length(events) + 1L]] <<- event
  expect_error(harmonise_sumstats(
    input, reference, strict = TRUE,
    observability = list(level = "events", callback = callback)
  ), "reference alignment failed")
  finalized <- Filter(function(x) identical(x$event, "finalized"), events)
  expect_length(finalized, 1L)
  expect_identical(finalized[[1L]]$status, "error")
  expect_true(all(vapply(events, function(x) {
    is.finite(x$elapsed_seconds) && x$elapsed_seconds >= 0
  }, logical(1L))))
})
