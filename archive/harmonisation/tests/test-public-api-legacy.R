test_that("import_sumstats exposes canonical columns and preserves extras", {
  input <- data.frame(
    `#CHROM` = c("chr1", "2"), POS = c(10L, 20L),
    ALT = c("C", "T"), REF = c("A", "G"),
    OR = c(1.2, 0.8), SE = c(0.1, 0.2), AF = c(0.2, 0.7),
    P = c(0.05, 0.01), marker = c("keep", "me"),
    check.names = FALSE
  )
  got <- import_sumstats(input)
  expect_equal(got$chromosome, c("1", "2"))
  expect_equal(got$base_pair_location, c(10L, 20L))
  expect_equal(got$beta, log(c(1.2, 0.8)))
  expect_equal(got$z, got$beta / got$standard_error)
  expect_equal(got$marker, c("keep", "me"))
  expect_equal(attr(got, "source_columns"), names(input))
})

test_that("compression exposes orthogonal issue-18 arguments and manifest provenance", {
  skip_if_not_installed("arrow")
  input <- make_fixture(12L)
  path <- tempfile("issue18-api-")
  store <- compress_sumstats(
    input, path, reference = NULL, backend = "parquet", profile = "exact",
    store_build = "hg38", selection = "full", threads = 2L,
    row_policy = "report", overwrite = TRUE
  )
  expect_s3_class(store, "compressor_store")
  expect_equal(store$manifest$selection_scope, "full")
  expect_equal(store$manifest$input_build, "GRCh38")
  expect_equal(store$manifest$stored_build, "GRCh38")
  expect_equal(store$manifest$profile, "exact")
  expect_true(nzchar(store$manifest$codec$name))
  expect_true(isTRUE(store$manifest$tolerances$exact))
  expect_equal(store$manifest$row_policy, "report")
  expect_equal(store$manifest$threads$requested, 2L)
  expect_equal(store$manifest$threads$effective, 1L)
  expect_equal(store$manifest$provenance$input$kind, "data.frame")
  expect_true(validate_compressor(store)$valid)
})

test_that("threads is a positive scalar and row policy can reject rows", {
  skip_if_not_installed("arrow")
  input <- make_fixture(2L)
  expect_error(
    compress_sumstats(input, tempfile("bad-threads-"), reference = NULL,
                      backend = "parquet", threads = 0L),
    "threads"
  )
  bad_reference <- input[1L, c("chromosome", "base_pair_location",
                               "effect_allele", "other_allele"), drop = FALSE]
  bad_reference$base_pair_location <- bad_reference$base_pair_location + 100L
  expect_error(
    compress_sumstats(input, tempfile("row-policy-"), reference = bad_reference,
                      backend = "parquet", selection = "full",
                      row_policy = "error", overwrite = TRUE),
    "reference alignment failed"
  )
})

test_that("liftover_sumstats is a public no-op for GRCh38", {
  input <- make_fixture(5L)
  got <- liftover_sumstats(input, input_build = "GRCh38")
  expect_equal(nrow(got), nrow(input))
  expect_true(all(got$liftover_status == "not_needed"))
  expect_equal(attr(got, "liftover_stats")$not_needed, 5L)
  expect_equal(attr(got, "genome_build"), "GRCh38")
  expect_false(any(grepl("^\\.compressor_liftover", names(got))))
})

test_that("liftover_sumstats accepts gzipped chains and reports status", {
  skip_if_not_installed("rtracklayer")
  chain <- tempfile(fileext = ".chain.gz")
  con <- gzfile(chain, open = "wt")
  writeLines(c(
    "chain 1 chr1 1000 + 0 1000 chr1 1000 + 0 1000 1",
    "1000",
    ""
  ), con)
  close(con)
  input <- data.frame(
    chromosome = "1", base_pair_location = 10L,
    other_allele = "A", effect_allele = "C",
    beta = 0.2, standard_error = 0.1,
    effect_allele_frequency = 0.3
  )
  got <- liftover_sumstats(input, input_build = "GRCh37", chain = chain)
  expect_equal(got$base_pair_location, 10L)
  expect_equal(got$liftover_status, "mapped")
  expect_equal(attr(got, "liftover_stats")$mapped, 1L)
})

test_that("the shipped example imports without network access", {
  path <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
  if (!nzchar(path)) path <- file.path("inst", "extdata", "example-grch38.tsv")
  got <- import_sumstats(path)
  expect_equal(nrow(got), 10L)
  expect_true(all(c("z", "beta", "standard_error", "effect_allele_frequency") %in% names(got)))
})
