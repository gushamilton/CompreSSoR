empty_prepared_input <- function() {
  data.frame(
    chromosome = character(),
    base_pair_location = integer(),
    reference_allele = character(),
    alternate_allele = character(),
    effect_allele = character(),
    other_allele = character(),
    beta = numeric(),
    standard_error = numeric(),
    effect_allele_frequency = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
}

expect_empty_input_failure <- function(input, output, row_policy) {
  expect_error(
    compress_sumstats(input, output, row_policy = row_policy),
    "input contains zero rows"
  )
  expect_false(file.exists(output) || dir.exists(output))
  expect_length(list.files(dirname(output), all.files = TRUE, no.. = TRUE), 0L)
}

test_that("empty prepared data frames fail before staging a store", {
  parent <- tempfile("compressor-empty-data-frame-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)

  for (row_policy in c("report", "error")) {
    output <- file.path(parent, paste0("empty-", row_policy, ".cpr"))
    expect_empty_input_failure(empty_prepared_input(), output, row_policy)
  }
})

test_that("empty prepared TSV files fail before staging a store", {
  parent <- tempfile("compressor-empty-tsv-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  input <- file.path(parent, "empty.tsv")
  writeLines(
    paste(
      "chromosome", "base_pair_location", "reference_allele",
      "alternate_allele", "effect_allele", "other_allele", "beta",
      "standard_error", "effect_allele_frequency", "p_value", sep = "\t"
    ),
    input
  )

  expect_empty_input_failure(
    input, file.path(parent, "stores", "empty-tsv.cpr"), row_policy = "report"
  )
})
