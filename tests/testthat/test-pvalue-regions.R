test_that("p-value regions merge overlapping padded windows", {
  input <- data.frame(
    chromosome = rep("1", 8),
    base_pair_location = c(50000L, 100000L, 150000L, 200000L,
                           250000L, 300000L, 350000L, 400000L),
    z = c(0, 5, 5, 0, 0, 5, 0, 0),
    stringsAsFactors = FALSE
  )
  selected <- CompreSSoR:::pvalue_region_selection(
    input, pvalue_threshold = 1e-5, region_padding = 50000L
  )

  expect_equal(selected$metadata$seed_snps, 3L)
  expect_equal(selected$metadata$regions, 2L)
  expect_equal(selected$metadata$kept_rows, 7L)
  expect_equal(selected$regions$start, c(50000L, 250000L))
  expect_equal(selected$regions$end, c(200000L, 350000L))
  expect_equal(input$base_pair_location[selected$keep],
               c(50000L, 100000L, 150000L, 200000L, 250000L, 300000L, 350000L))
})

test_that("p-value thresholds use the pre-encoding harmonised statistic", {
  input <- data.frame(
    chromosome = c("1", "1"), base_pair_location = c(100L, 200L),
    z = c(5, 0), p_value = c(0.9, 1e-20), stringsAsFactors = FALSE
  )
  selected <- CompreSSoR:::pvalue_region_selection(
    input, pvalue_threshold = 1e-5, region_padding = 0L
  )

  expect_equal(selected$keep, c(TRUE, FALSE))
  expect_equal(selected$metadata$threshold_source, "pre_encoding_harmonised")
  expect_equal(selected$metadata$threshold_statistic, "p_value_from_harmonised_z")
  expect_equal(selected$metadata$threshold_operator, "<")
  expect_equal(selected$metadata$threshold_semantics$encoding, "not_encoded")
  expect_equal(selected$metadata$threshold_semantics$derivation,
               "2 * pnorm(-abs(z))")
})

test_that("p-value regions are selected independently across chromosomes", {
  input <- data.frame(
    chromosome = c("2", "1", "X", "Y", "2"),
    base_pair_location = c(120000L, 200000L, 300000L, 400000L, 160000L),
    z = c(0, 5, 0, 5, 5),
    stringsAsFactors = FALSE
  )
  selected <- CompreSSoR:::pvalue_region_selection(
    input, pvalue_threshold = 1e-5, region_padding = 50000L
  )

  expect_equal(selected$metadata$seed_snps, 3L)
  expect_equal(selected$metadata$regions, 3L)
  expect_equal(selected$keep, c(TRUE, TRUE, FALSE, TRUE, TRUE))
})

test_that("p-value region selection records a compact storage tag", {
  input <- data.frame(
    chromosome = rep("1", 8),
    base_pair_location = c(50000L, 100000L, 150000L, 200000L,
                           250000L, 300000L, 350000L, 400000L),
    effect_allele = rep("A", 8), other_allele = rep("C", 8),
    z = c(0, 5, 5, 0, 0, 5, 0, 0),
    beta = c(0, 0.5, 0.5, 0, 0, 0.5, 0, 0),
    standard_error = rep(0.1, 8),
    effect_allele_frequency = rep(0.3, 8),
    stringsAsFactors = FALSE
  )
  skip_if_not_installed("arrow")
  path <- tempfile("pvalue-regions-", fileext = ".cpr")
  store <- CompreSSoR::compress_sumstats(
    input, path, reference = NULL, mode = "pvalue_regions",
    backend = "parquet", overwrite = TRUE
  )

  expect_equal(store$manifest$selection$tag, "core")
  expect_equal(store$manifest$selection$pvalue_threshold, 1e-5)
  expect_equal(store$manifest$selection$padding_bp, 50000L)
  expect_equal(store$manifest$selection$seed_snps, 3L)
  expect_equal(store$manifest$selection$regions, 2L)
  expect_true(file.exists(file.path(path, "core_regions.json")))
  expect_equal(nrow(CompreSSoR::read_sumstats(store)), 7L)
  expect_true(CompreSSoR::validate_compressor(store)$valid)
})

test_that("p-value region arguments are validated", {
  input <- data.frame(
    chromosome = "1", base_pair_location = 100L, z = 5,
    stringsAsFactors = FALSE
  )
  expect_error(CompreSSoR:::pvalue_region_selection(input, pvalue_threshold = 0),
               "pvalue_threshold")
  expect_error(CompreSSoR:::pvalue_region_selection(input, region_padding = 1.5),
               "region_padding")
})

test_that("core_plus stores the core panel and significant regions", {
  input <- data.frame(
    chromosome = rep("1", 8),
    base_pair_location = c(50000L, 100000L, 150000L, 200000L,
                           250000L, 300000L, 350000L, 400000L),
    effect_allele = rep("A", 8), other_allele = rep("C", 8),
    z = c(0, 5, 5, 0, 0, 5, 0, 0),
    beta = c(0, 0.5, 0.5, 0, 0, 0.5, 0, 0),
    standard_error = rep(0.1, 8),
    effect_allele_frequency = rep(0.3, 8),
    stringsAsFactors = FALSE
  )
  skip_if_not_installed("arrow")
  reference <- input[c("chromosome", "base_pair_location",
                       "effect_allele", "other_allele")]
  panel <- input[8L, c("chromosome", "base_pair_location",
                       "effect_allele", "other_allele")]
  path <- tempfile("core-plus-", fileext = ".cpr")
  store <- CompreSSoR::compress_sumstats(
    input, path, reference = reference, mode = "core_plus",
    variant_set = panel, backend = "parquet", overwrite = TRUE
  )

  expect_equal(store$manifest$selection$tag, "core_plus")
  expect_equal(store$manifest$selection$core_variant_rows, 1L)
  expect_equal(store$manifest$selection$seed_snps, 3L)
  expect_equal(store$manifest$selection$kept_rows, 8L)
  expect_true(file.exists(file.path(path, "core_plus_regions.json")))
  expect_equal(nrow(CompreSSoR::read_sumstats(store)), 8L)
  expect_true(CompreSSoR::validate_compressor(store)$valid)
})

test_that("core-plus can filter an already harmonised object without a second reference pass", {
  input <- data.frame(
    chromosome = rep("1", 8),
    base_pair_location = c(50000L, 100000L, 150000L, 200000L,
                           250000L, 300000L, 350000L, 400000L),
    effect_allele = rep("A", 8), other_allele = rep("C", 8),
    z = c(0, 5, 5, 0, 0, 5, 0, 0),
    beta = c(0, 0.5, 0.5, 0, 0, 0.5, 0, 0),
    standard_error = rep(0.1, 8),
    effect_allele_frequency = rep(0.3, 8),
    stringsAsFactors = FALSE
  )
  reference <- input[c("chromosome", "base_pair_location",
                       "effect_allele", "other_allele")]
  panel <- input[8L, c("chromosome", "base_pair_location",
                       "effect_allele", "other_allele")]
  harmonised <- CompreSSoR::harmonise_sumstats(input, reference, mode = "qc")
  skip_if_not_installed("arrow")
  path <- tempfile("core-plus-pre-harmonised-", fileext = ".cpr")
  store <- CompreSSoR::compress_sumstats(
    harmonised, path, reference = NULL, mode = "core_plus",
    variant_set = panel, backend = "parquet", overwrite = TRUE
  )

  expect_equal(store$manifest$harmonisation$method, "pre_harmonised")
  expect_equal(store$manifest$selection$tag, "core_plus")
  expect_equal(store$manifest$selection$kept_rows, 8L)
  expect_true("variant_set" %in% names(store$manifest$harmonisation$alignment))
  expect_equal(nrow(CompreSSoR::read_sumstats(store)), 8L)
})

test_that("canonical panel membership is allele-aware", {
  data <- data.frame(
    variant_id = c("1:100:A:C", "1:100:A:G", "1:101:C:T"),
    chromosome = "1", base_pair_location = c(100L, 100L, 101L),
    other_allele = c("A", "A", "C"),
    effect_allele = c("C", "G", "T"),
    rsid = c("rs1", "rs1", "rs2"), stringsAsFactors = FALSE
  )
  panel <- data.frame(variant_id = "chr1:100:a:c", stringsAsFactors = FALSE)
  expect_identical(
    CompreSSoR:::variant_set_membership(data, CompreSSoR:::normalise_variant_set_columns(panel)),
    c(TRUE, FALSE, FALSE)
  )
})

test_that("compressed panel readers select identity columns before normalization", {
  skip_if_not_installed("data.table")
  panel_path <- tempfile("panel-", fileext = ".tsv.gz")
  panel_data <- data.frame(
    vid = c("1:100:A:C", "1:200:G:T"),
    chromosome = c("1", "1"), base_pair_location = c(100L, 200L),
    other_allele = c("A", "G"), effect_allele = c("C", "T"),
    large_unused_annotation = c("x", "y"), stringsAsFactors = FALSE
  )
  con <- gzfile(panel_path, open = "wt")
  write.table(panel_data, con, sep = "\t", row.names = FALSE,
              quote = FALSE)
  close(con)
  panel <- CompreSSoR:::read_variant_set(panel_path)
  expect_equal(attr(panel, "variant_set_metadata")$columns_read,
               "vid")
  expect_false("large_unused_annotation" %in% names(panel))
  expect_false("vid" %in% names(panel))
  expect_equal(nrow(panel), 2L)
})
