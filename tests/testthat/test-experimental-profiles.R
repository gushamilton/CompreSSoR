# Experimental quantisation sweep options (bench/quant-frontier). Opt-in only;
# the standard grammar and the default store are unchanged.

test_that("experimental profiles need the opt-in option", {
  expect_error(CompreSSoR:::pcodec_native_profile("z8/eaf8/se6"), "support")
  expect_error(CompreSSoR:::pcodec_native_profile("z10/eaf6/se8+xse"), "support")
  expect_error(CompreSSoR:::pcodec_native_profile("z10/eaf8/se8+xse@5"), "support")
  old <- options(CompreSSoR.experimental_profiles = TRUE)
  on.exit(options(old))
  p <- CompreSSoR:::pcodec_native_profile("z14/eaf6/se10+xse@5")
  expect_identical(c(p$z_bits, p$eaf_bits, p$se_bits), c(14L, 6L, 10L))
  expect_identical(p$z_range, c(-5, 5))
  expect_identical(p$eaf_count, 63L)
  expect_identical(p$se_physical_dtype, "uint16")
  expect_error(CompreSSoR:::pcodec_native_profile("z16/eaf8/se8"), "support")
  # The standard profiles parse exactly as before.
  d <- CompreSSoR:::pcodec_native_profile("z10/eaf8/se8+xse")
  expect_identical(d$codec_name, "pcodec_native_standalone_z10_eaf8_se8_xse_zstd_exceptions")
  expect_identical(d$se_physical_dtype, "uint8")
  expect_identical(d$eaf_count, 255L)
})

test_that("experimental profiles round-trip through full and keyed reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  input <- make_fixture(3000L)
  old <- options(CompreSSoR.experimental_profiles = TRUE)
  on.exit(options(old))
  cols <- c("chromosome", "base_pair_location", "other_allele", "effect_allele",
            "beta", "standard_error", "effect_allele_frequency", "p_value")
  for (name in c("z14/eaf6/se10+xse@5", "z6/eaf4/se4", "z8/eaf8/se12+xse@2.5")) {
    options(CompreSSoR.native_profile = name)
    path <- tempfile("experimental-profile-")
    store <- compress_sumstats(input, path, overwrite = TRUE)
    options(CompreSSoR.native_profile = NULL)
    p <- CompreSSoR:::pcodec_native_profile(name)
    expect_identical(store$manifest$semantic_codec$name, name)
    expect_identical(store$manifest$semantic_codec$se_physical_dtype, p$se_physical_dtype)
    full <- read_sumstats(path, columns = cols)
    z <- input$beta / input$standard_error
    key_in <- paste(input$chromosome, input$base_pair_location, input$other_allele, input$effect_allele)
    key_out <- paste(full$chromosome, full$base_pair_location, full$other_allele, full$effect_allele)
    m <- match(key_out, key_in)
    expect_false(anyNA(m))
    zerr <- abs(full$beta / full$standard_error - z[m])
    central <- abs(z[m]) < p$z_range[2]
    half_step <- diff(p$z_range) / (2 * p$z_count)
    expect_true(max(zerr[central]) <= half_step * (1 + 1e-6) + 1e-12)
    keys <- compressor_variant_key(full$chromosome[1:50], full$base_pair_location[1:50],
                                   full$other_allele[1:50], full$effect_allele[1:50])
    sel <- read_sumstats(path, variants = keys, columns = cols)
    expect_identical(sel$standard_error, full$standard_error[match(
      paste(sel$chromosome, sel$base_pair_location), paste(full$chromosome, full$base_pair_location))])
  }
})

test_that("experimental lossless value streams are exact and leave standard reads alone", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  input <- make_fixture(2000L)
  a <- tempfile("lossless-off-"); b <- tempfile("lossless-on-")
  compress_sumstats(input, a, overwrite = TRUE)
  old <- options(CompreSSoR.experimental_lossless_values = TRUE)
  compress_sumstats(input, b, overwrite = TRUE)
  options(old)
  for (f in c("z.pco", "se.pco", "eaf.pco", "position.pco", "substitution.pco", "exceptions.bin")) {
    expect_identical(unname(tools::md5sum(file.path(a, f))), unname(tools::md5sum(file.path(b, f))))
  }
  expect_identical(as.list(read_sumstats(a)), as.list(read_sumstats(b)))
  lv <- CompreSSoR:::pcodec_native_read_lossless_values(b)
  expect_identical(sort(lv$beta), sort(input$beta))
  expect_identical(sort(lv$standard_error), sort(input$standard_error))
  expect_error(CompreSSoR:::pcodec_native_read_lossless_values(a), "no experimental lossless")
})
