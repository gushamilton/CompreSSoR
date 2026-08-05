test_that("standard stores use semantic codes and omit P", {
  skip_if_not_installed("arrow")
  input <- make_fixture(120L)
  # This test targets storage semantics, not palindromic harmonisation.
  input$effect_allele <- "C"
  input$other_allele <- "A"
  reference <- input[c("chromosome", "base_pair_location", "effect_allele",
                       "other_allele", "rsid", "effect_allele_frequency")]
  path <- file.path(tempdir(), "standard.cpr")
  store <- compress_sumstats(input, path, reference = reference,
                             keep_extras = TRUE, backend = "parquet", overwrite = TRUE)
  expect_s3_class(store, "compressor_store")
  expect_equal(store$manifest$codec$name, "semantic_z9_eaf8_se6")
  expect_equal(unlist(store$manifest$logical_columns), c("z", "standard_error", "effect_allele_frequency"))
  expect_equal(store$manifest$derived_columns$p_value, "2 * pnorm(-abs(z))")
  expect_true(validate_compressor(store)$valid)
  values <- arrow::read_parquet(file.path(path, "values.parquet"))
  expect_false("p_value" %in% names(values))
  expect_false("beta" %in% names(values))
  got <- read_sumstats(store, region = "chr1:100100-100150",
                       columns = c("variant_id", "rsid", "z", "beta", "standard_error", "p_value", "annotation"))
  expect_equal(nrow(got), 21L)
  expect_true(all(c("z", "beta", "standard_error", "p_value", "annotation") %in% names(got)))
  expect_true(all(grepl("^1:[0-9]+:[ACGT]:[ACGT]$", got$variant_id)))
  light <- read_sumstats(store, region = "chr1:100100-100150",
                         columns = c("variant_id", "z", "standard_error",
                                     "effect_allele_frequency"))
  expect_equal(names(light), c("variant_id", "z", "standard_error",
                               "effect_allele_frequency"))
  expect_false("beta" %in% names(light))
  expect_false("p_value" %in% names(light))
})

test_that("Parquet stores retain inline identity after reference-backed QC", {
  skip_if_not_installed("arrow")
  input <- make_fixture(40L)
  input <- input[seq(1L, nrow(input), by = 4L), , drop = FALSE]
  reference <- input[c("chromosome", "base_pair_location", "effect_allele",
                       "other_allele", "rsid")]
  reference_path <- tempfile(fileext = ".parquet")
  build_canonical_reference(reference, reference_path, overwrite = TRUE)
  path <- tempfile("shared-reference-")
  store <- compress_sumstats(input, path, reference = reference_path,
                             profile = "exact", backend = "parquet", overwrite = TRUE)
  index <- arrow::read_parquet(file.path(path, "variants.parquet"))
  expect_equal(store$manifest$variant_storage, "inline")
  expect_true(all(c("row", "chromosome", "base_pair_location",
                    "effect_allele", "other_allele") %in% names(index)))
  expect_false("reference_index" %in% names(index))
  expect_equal(read_sumstats(store)$variant_id,
               paste(reference$chromosome, reference$base_pair_location,
                     reference$other_allele, reference$effect_allele, sep = ":"))
  expect_true(validate_compressor(store)$valid)
})

test_that("exact profile round-trips Z and derives P without storing it", {
  skip_if_not_installed("arrow")
  input <- make_fixture(100L)
  path <- file.path(tempdir(), "exact.cpr")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             profile = "exact", backend = "parquet", overwrite = TRUE)
  got <- read_sumstats(store)
  expect_equal(got$z, input$beta / input$standard_error, tolerance = 0)
  expect_equal(got$beta, input$beta, tolerance = 1e-15)
  expect_equal(got$p_value, input$p_value, tolerance = 0)
  expect_false("p_value" %in% names(arrow::read_parquet(file.path(path, "values.parquet"))))
})

test_that("reverse and complementary alleles are harmonised", {
  input <- data.frame(chromosome = c("1", "1"), base_pair_location = c(10L, 20L),
                      effect_allele = c("G", "T"), other_allele = c("C", "G"),
                      beta = c(.2, .3), standard_error = c(.1, .1),
                      effect_allele_frequency = c(.2, .2), rsid = c("rs1", "rs2"))
  reference <- data.frame(chromosome = c("1", "1"), base_pair_location = c(10L, 20L),
                          effect_allele = c("C", "A"), other_allele = c("G", "C"),
                          rsid = c("rs1", "rs2"))
  got <- harmonise_sumstats(input, reference)
  expect_equal(nrow(got), 1L)
  expect_equal(got$beta, .3)
  expect_equal(got$effect_allele, "A")
  expect_equal(got$other_allele, "C")
})

test_that("validation detects broken store files", {
  skip_if_not_installed("arrow")
  input <- make_fixture(20L)
  path <- file.path(tempdir(), "broken-store.cpr")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             profile = "exact", backend = "parquet", overwrite = TRUE)
  unlink(file.path(path, "values.parquet"))
  result <- validate_compressor(path)
  expect_false(result$valid)
  expect_true(any(grepl("values.parquet", result$errors)))
})
