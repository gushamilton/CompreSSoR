test_that("native full-scan decoder matches the validated R decoder", {
  skip_if_not(is.loaded("compressor_decode_native", PACKAGE = "CompreSSoR"))
  n <- 257L
  beta <- c(seq(-4, 4, length.out = 40), sin(seq_len(n - 40) / 7) / 3)
  se <- 0.02 + (seq_len(n) %% 19) / 1000
  eaf <- (seq_len(n) %% 97) / 100
  eaf[c(1, 3, 121)] <- c(1, NA_real_, NA_real_)
  encoded <- CompreSSoR:::q_encode(beta, se, eaf, block_rows = 64L)
  old <- getOption("CompreSSoR.native_decode", TRUE)
  on.exit(options(CompreSSoR.native_decode = old), add = TRUE)

  options(CompreSSoR.native_decode = FALSE)
  fallback <- CompreSSoR:::q_decode(encoded$main, encoded$exceptions,
                                    encoded$metadata, include_beta = TRUE,
                                    include_p = TRUE)
  options(CompreSSoR.native_decode = TRUE)
  native <- CompreSSoR:::q_decode(encoded$main, encoded$exceptions,
                                  encoded$metadata, include_beta = TRUE,
                                  include_p = TRUE)

  expect_equal(native$z, fallback$z, tolerance = 1e-12)
  expect_equal(native$standard_error, fallback$standard_error, tolerance = 1e-12)
  expect_equal(native$effect_allele_frequency, fallback$effect_allele_frequency,
               tolerance = 1e-12)
  expect_equal(native$beta, fallback$beta, tolerance = 1e-12)
  expect_equal(native$p_value, fallback$p_value, tolerance = 1e-12)
})

test_that("native decoder handles malformed sentinel codes without indexing out of bounds", {
  skip_if_not(is.loaded("compressor_decode_native", PACKAGE = "CompreSSoR"))
  encoded <- CompreSSoR:::q_encode(
    beta = c(0.1, 0.2, 0.3, 0.4),
    se = c(0.02, 0.03, 0.04, 0.05),
    eaf = c(0.2, 0.3, 0.4, 0.5),
    block_rows = 2L
  )
  encoded$main$z_code <- c(-1L, encoded$metadata$z_count, 2^encoded$metadata$z_bits, 0L)
  encoded$main$se_code <- c(-1L, 0L, encoded$metadata$se_count, 0L)
  encoded$main$eaf_code <- c(-1L, 0L, encoded$metadata$eaf_count + 1L, 0L)
  empty_exceptions <- data.frame(
    row = integer(), z_value = numeric(), se_value = numeric(),
    eaf_value = numeric(), flags = integer()
  )
  old <- getOption("CompreSSoR.native_decode", TRUE)
  on.exit(options(CompreSSoR.native_decode = old), add = TRUE)
  options(CompreSSoR.native_decode = TRUE)
  decoded <- CompreSSoR:::q_decode(encoded$main, empty_exceptions,
                                   encoded$metadata, include_beta = TRUE,
                                   include_p = TRUE)
  expect_length(decoded$z, 4L)
  expect_true(is.na(decoded$z[1]))
  expect_true(is.na(decoded$z[2]))
  expect_true(is.na(decoded$standard_error[1]))
  expect_true(is.na(decoded$effect_allele_frequency[1]))
})

test_that("native decoder rejects unsupported code domains", {
  skip_if_not(is.loaded("compressor_decode_native", PACKAGE = "CompreSSoR"))
  encoded <- CompreSSoR:::q_encode(
    beta = c(0.1, 0.2), se = c(0.02, 0.03), eaf = c(0.2, 0.3),
    block_rows = 2L
  )
  encoded$metadata$z_bits <- 17L
  expect_error(
    CompreSSoR:::q_decode(encoded$main, encoded$exceptions,
                           encoded$metadata, include_beta = FALSE,
                           include_p = FALSE),
    "16-bit limit"
  )
})
