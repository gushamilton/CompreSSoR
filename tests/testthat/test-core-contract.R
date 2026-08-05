test_that("the strict core requires explicit REF, ALT and prepared orientation", {
  input <- make_fixture(8L)
  missing_identity <- input[c("chromosome", "base_pair_location", "effect_allele",
                              "other_allele", "beta", "standard_error")]
  expect_error(
    compress_sumstats(missing_identity, tempfile("missing-ref-alt-")),
    "requires explicit REF and ALT"
  )

  flipped <- input
  flipped$effect_allele <- flipped$reference_allele
  flipped$other_allele <- flipped$alternate_allele
  expect_error(
    compress_sumstats(flipped, tempfile("flipped-")),
    "inconsistent with explicit REF/ALT"
  )
})

test_that("build conversion is outside the compression core", {
  expect_error(
    compress_sumstats(make_fixture(4L), tempfile("cross-build-"),
                      input_build = "GRCh37", store_build = "GRCh38"),
    "must be the same build"
  )
})

test_that("GRCh37 is a first-class native identity build", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(128L)
  path <- tempfile("grch37-native-")
  store <- compress_sumstats(input, path, input_build = "hg19",
                             store_build = "GRCh37", threads = 2L,
                             overwrite = TRUE)
  expect_identical(store$manifest$genome_build, "GRCh37")
  expect_identical(store$manifest$identity$chromosome_table$id,
                   "grch37_primary_1_22_X_Y")
  expect_identical(store$manifest$identity$input_build, "GRCh37")
  expect_true(validate_compressor(store, full = TRUE)$valid)
  got <- read_sumstats(store, columns = c("chromosome", "base_pair_location",
                                           "reference_allele", "alternate_allele"))
  expect_equal(got$base_pair_location, sort(input$base_pair_location))
  expect_equal(got$reference_allele,
               input$reference_allele[order(input$base_pair_location)])
  expect_equal(got$alternate_allele,
               input$alternate_allele[order(input$base_pair_location)])
})

test_that("core manifests make external preparation explicit", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  store <- compress_sumstats(make_fixture(6L), tempfile("manifest-contract-"),
                             input_build = "GRCh38", store_build = "GRCh38",
                             overwrite = TRUE)
  expect_identical(store$manifest$reference$status, "not_used")
  expect_identical(store$manifest$identity$chain$status, "not_used")
  expect_identical(store$manifest$preparation$method, "strict_prepared_input")
  expect_false("harmonisation" %in% names(store$manifest$threads))
})
