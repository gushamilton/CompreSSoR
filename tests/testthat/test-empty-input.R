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

expect_empty_input_failure <- function(input, output, row_policy,
                                       backend = "pcodec", qc = "compact",
                                       profile = "standard",
                                       expect_parent_absent = FALSE) {
  args <- list(
    input = input,
    output = output,
    row_policy = row_policy,
    backend = backend,
    qc = qc,
    profile = profile
  )
  expect_error(
    do.call(compress_sumstats, args),
    "input contains zero rows"
  )
  expect_false(file.exists(output) || dir.exists(output))
  if (isTRUE(expect_parent_absent)) {
    expect_false(dir.exists(dirname(output)))
  } else {
    expect_length(list.files(dirname(output), all.files = TRUE, no.. = TRUE), 0L)
  }
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

  # The empty input is rejected before the output parent is created.
  output <- file.path(parent, "stores", "empty-tsv.cpr")
  expect_empty_input_failure(input, output, row_policy = "report",
                             expect_parent_absent = TRUE)
})

test_that("empty prepared TSV.GZ files use the same diagnostic", {
  parent <- tempfile("compressor-empty-tsvgz-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)
  input <- file.path(parent, "empty.tsv.gz")
  con <- gzfile(input, open = "wt")
  writeLines(
    paste(
      "chr", "pos", "ref", "alt", "a1", "a2", "beta", "se", "eaf", "p",
      sep = "\t"
    ),
    con
  )
  close(con)

  output <- file.path(parent, "stores", "empty-tsvgz.cpr")
  expect_empty_input_failure(input, output, row_policy = "error",
                             expect_parent_absent = TRUE)
})

test_that("empty canonical input is rejected consistently for both backends and QC modes", {
  parent <- tempfile("compressor-empty-modes-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)

  for (backend in c("pcodec", "parquet")) {
    if (identical(backend, "parquet")) skip_if_not_installed("arrow")
    profile <- if (identical(backend, "parquet")) "exact" else "standard"
    for (row_policy in c("report", "error")) {
      output <- file.path(parent, paste0(backend, "-compact-", row_policy, ".cpr"))
      expect_empty_input_failure(
        empty_prepared_input(), output, row_policy = row_policy,
        backend = backend, qc = "compact", profile = profile
      )
    }
  }

  output <- file.path(parent, "pcodec-none.cpr")
  expect_empty_input_failure(
    empty_prepared_input(), output, row_policy = "report",
    backend = "pcodec", qc = "none", profile = "standard"
  )
})

test_that("empty aliased input is rejected before alias/schema resolution", {
  input <- data.frame(
    chr = character(), pos = integer(), ref = character(), alt = character(),
    a1 = character(), a2 = character(), beta = numeric(), se = numeric(),
    eaf = numeric(), p = numeric(), stringsAsFactors = FALSE
  )

  for (row_policy in c("report", "error")) {
    expect_error(
      import_sumstats(input, row_policy = row_policy),
      "input contains zero rows"
    )
    expect_error(
      preflight_sumstats(input, row_policy = row_policy),
      "input contains zero rows"
    )
  }
})

test_that("fully rejected non-empty input remains a structural-QC diagnostic", {
  input <- make_fixture(1L)
  input$reference_allele <- "AT"
  input$alternate_allele <- "G"
  input$effect_allele <- "G"
  input$other_allele <- "AT"

  expect_error(
    compress_sumstats(input, tempfile("all-rejected-empty-regression-"),
                      row_policy = "report", overwrite = TRUE),
    "structural QC rejected 1 row\\(s\\)"
  )
  expect_error(
    compress_sumstats(input, tempfile("all-rejected-error-regression-"),
                      row_policy = "error", overwrite = TRUE),
    "structural QC rejected 1 row\\(s\\)"
  )
})
