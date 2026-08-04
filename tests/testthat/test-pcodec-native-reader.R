test_that("native stream reader returns compact keys without identity expansion", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(5000L)
  path <- tempfile("pcodec-native-reader-")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)

  compact <- read_sumstats(
    store,
    columns = c("global_position", "substitution", "z", "standard_error",
                "effect_allele_frequency")
  )
  expect_identical(names(compact), c("global_position", "substitution", "z",
                                     "standard_error", "effect_allele_frequency"))
  expect_equal(nrow(compact), nrow(input))
  expect_true(is.double(compact$global_position))
  expect_true(is.integer(compact$substitution))
  identity <- CompreSSoR:::pcodec_native_identity(input)
  order <- order(identity$global_position, identity$substitution, method = "radix")
  expect_equal(compact$global_position, identity$global_position[order])
  expect_equal(compact$substitution, identity$substitution[order])

  old <- getOption("CompreSSoR.pcodec.native_stream_reader")
  options(CompreSSoR.pcodec.native_stream_reader = FALSE)
  fallback <- read_sumstats(
    store,
    columns = c("z", "standard_error", "effect_allele_frequency")
  )
  options(CompreSSoR.pcodec.native_stream_reader = old)
  expect_equal(compact$z, fallback$z)
  expect_equal(compact$standard_error, fallback$standard_error)
  expect_equal(compact$effect_allele_frequency,
               fallback$effect_allele_frequency)
})
