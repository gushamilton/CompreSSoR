test_that("native Pcodec is available and round trips integer streams", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  for (spec in list(
    list(dtype = "u8", values = c(0L, 1L, 7L, 127L, 255L)),
    list(dtype = "u16", values = c(0L, 1L, 255L, 4096L, 65535L)),
    list(dtype = "u32", values = c(0, 1, 2147483648, 4000000000, 4294967295))
  )) {
    encoded <- CompreSSoR:::pcodec_native_compress(spec$values, spec$dtype)
    decoded <- CompreSSoR:::pcodec_native_decompress(encoded, length(spec$values), spec$dtype)
    expect_identical(decoded, spec$values)
  }
})

test_that("native Pcodec selects access-appropriate thread defaults", {
  expect_equal(CompreSSoR:::pcodec_native_default_threads(), 4L)
  expect_equal(CompreSSoR:::pcodec_native_default_threads(region = "chr1:1-10"), 1L)
  expect_equal(CompreSSoR:::pcodec_native_default_threads(variants = 0L), 1L)
  expect_equal(CompreSSoR:::pcodec_native_default_threads(threads = 2L), 2L)
  old <- getOption("CompreSSoR.pcodec.threads")
  on.exit(options(CompreSSoR.pcodec.threads = old), add = TRUE)
  options(CompreSSoR.pcodec.threads = 3L)
  expect_equal(CompreSSoR:::pcodec_native_default_threads(), 3L)
})

test_that("legacy Z9/SE6 constants are unchanged", {
  # Readers fall back to these for manifests that predate explicit counts.
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_PROFILE, "z9/eaf8/se6")
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_BITS, 6L)
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_COUNT, 62L)
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_MISSING_CODE, 62L)
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_EXCEPTION_CODE, 63L)
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_PHYSICAL_DTYPE, "uint8")
  expect_identical(CompreSSoR:::PCODEC_NATIVE_SE_PHYSICAL_BITS, 8L)
  expect_identical(CompreSSoR:::PCODEC_NATIVE_Z_BITS_LEGACY, 9L)
  expect_identical(
    CompreSSoR:::PCODEC_NATIVE_CODEC_NAME,
    "pcodec_native_standalone_z9_eaf8_se6_zstd_exceptions"
  )
  legacy <- CompreSSoR:::pcodec_native_profile("z9/eaf8/se6")
  expect_identical(legacy$z_count, 510L)
  expect_identical(legacy$se_count, 62L)
  expect_identical(legacy$codec_name, CompreSSoR:::PCODEC_NATIVE_CODEC_NAME)
  params <- CompreSSoR:::pcodec_native_semantic_params(list(
    name = "z9/eaf8/se6", z_bits = 9L, se_bits = 6L, se_count = 62L))
  expect_identical(params$z_count, 510L)
  expect_identical(params$se_count, 62L)
  expect_identical(params$eaf_count, 255L)
  expect_error(CompreSSoR:::pcodec_native_profile("z8/eaf8/se6"), "support")
  expect_error(CompreSSoR:::pcodec_native_profile("z10/eaf8/se9"), "support")
  expect_error(CompreSSoR:::pcodec_native_profile("standard"), "unknown")
  expect_error(CompreSSoR:::pcodec_native_profile("z10/eaf8/se8+foo"), "unknown")
  xse <- CompreSSoR:::pcodec_native_profile("z9/eaf8/se6+xse")
  expect_identical(xse$exception_se, "exact")
  expect_identical(xse$codec_name, "pcodec_native_standalone_z9_eaf8_se6_xse_zstd_exceptions")
})

test_that("native format and manifest expose the semantic profile contract", {
  expect_identical(CompreSSoR:::PCODEC_NATIVE_FORMAT, "0.4.6-pcodec-native")
  expect_identical(CompreSSoR:::PCODEC_NATIVE_DEFAULT_PROFILE, "z10/eaf8/se8+xse")
  default <- CompreSSoR:::pcodec_native_profile()
  expect_identical(default$name, "z10/eaf8/se8+xse")
  expect_identical(default$exception_se, "exact")
  expect_identical(CompreSSoR:::pcodec_native_profile("z9/eaf8/se6")$exception_se, "quantised")
  expect_identical(default$z_count, 1022L)
  expect_identical(default$se_count, 254L)

  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(32L)
  for (name in c(CompreSSoR:::PCODEC_NATIVE_DEFAULT_PROFILE, "z9/eaf8/se6", "z9/eaf8/se7")) {
    old <- options(CompreSSoR.native_profile = name)
    store <- compress_sumstats(input, tempfile("pcodec-native-provenance-"),
                               overwrite = TRUE)
    options(old)
    profile <- CompreSSoR:::pcodec_native_profile(name)
    semantic <- store$manifest$semantic_codec
    expect_identical(store$manifest$format_version, CompreSSoR:::PCODEC_NATIVE_FORMAT)
    expect_identical(semantic$name, name)
    expect_identical(semantic$z_bits, profile$z_bits)
    expect_identical(semantic$z_count, profile$z_count)
    expect_identical(semantic$z_missing, profile$z_count)
    expect_identical(semantic$z_exception, profile$z_count + 1L)
    expect_identical(semantic$se_bits, profile$se_bits)
    expect_identical(semantic$se_count, profile$se_count)
    expect_identical(semantic$se_missing, profile$se_count)
    expect_identical(semantic$se_exception, profile$se_count + 1L)
    expect_identical(semantic$eaf_count, 255L)
    expect_identical(semantic$exception_se, profile$exception_se)
    expect_identical(semantic$se_physical_dtype, CompreSSoR:::PCODEC_NATIVE_SE_PHYSICAL_DTYPE)
    expect_identical(semantic$se_physical_bits, CompreSSoR:::PCODEC_NATIVE_SE_PHYSICAL_BITS)
    expect_identical(store$manifest$codec$name, profile$codec_name)
    expect_identical(store$manifest$codec$se_bits, profile$se_bits)
    expect_identical(store$manifest$tolerances$se_profile, name)
    expect_identical(store$manifest$tolerances$se_physical_storage,
                     paste0("uint8 container for semantic SE", profile$se_bits, " codes"))
    expect_equal(store$manifest$tolerances$z_abs_max_central, 7 / (2 * profile$z_count))
    expect_equal(store$manifest$tolerances$se_relative_max,
                 2^(4 / profile$se_count) - 1)
  }
})

test_that("every reader honours the stored Z/SE profile", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- legacy_store_fixture_data(4000L)
  z_true <- input$beta / input$standard_error
  reads <- list()
  for (name in c("z9/eaf8/se6", "z9/eaf8/se8", "z10/eaf8/se7+xse", "z12/eaf8/se8")) {
    profile <- CompreSSoR:::pcodec_native_profile(name)
    old <- options(CompreSSoR.native_profile = name)
    path <- tempfile("profile-read-")
    compress_sumstats(input, path, qc = "none", pvalue_flag = TRUE,
                      pvalue_order = TRUE, pvalue_order_threshold = 0.05,
                      overwrite = TRUE)
    options(old)
    got <- legacy_store_reads(path)
    full <- got$full
    key_in <- paste(input$chromosome, input$base_pair_location, input$other_allele)
    at <- match(paste(full$chromosome, full$base_pair_location, full$other_allele), key_in)
    central <- abs(z_true[at]) < 3.5
    # Central Z within half a bin; SE within half a residual bin (plus EAF8).
    expect_lte(max(abs(full$z - z_true[at])[central]), 7 / (2 * profile$z_count) + 1e-12)
    rel_se <- abs(full$standard_error / input$standard_error[at] - 1)
    expect_lte(stats::median(rel_se), 2^(1 / profile$se_count) - 1)
    if (identical(profile$exception_se, "exact")) {
      # Z-exception rows take SE from their float32 exception record.
      expect_lt(max(rel_se[!central]), 1e-6)
    } else {
      expect_gt(max(rel_se[!central]), 1e-4)
    }
    expect_equal(full$p_value, 2 * pnorm(-abs(full$z)), tolerance = 1e-12)
    # Every other reader agrees with the full read.
    expect_identical(got$full_threads, full)
    region <- full$chromosome == "1" & full$base_pair_location >= 100201 &
      full$base_pair_location <= 101400
    expect_equal(as.data.frame(got$region)[names(full)],
                 as.data.frame(full[region, ]), ignore_attr = TRUE)
    for (t in c("candidates_1e3", "candidates_1e5")) {
      threshold <- if (t == "candidates_1e3") 1e-3 else 1e-5
      expect_identical(sort(got[[t]]$base_pair_location + 1e7 * (got[[t]]$chromosome == "2")),
                       sort((full$base_pair_location + 1e7 * (full$chromosome == "2"))[
                         !is.na(full$p_value) & full$p_value <= threshold]),
                       label = paste(name, t))
    }
    reads[[name]] <- got
  }
  # p-value domains come from the exact input, not the stored codes.
  for (name in names(reads)[-1]) {
    expect_identical(reads[[name]]$flag, reads[[1]]$flag)
    expect_identical(reads[[name]]$order, reads[[1]]$order)
    expect_identical(reads[[name]]$candidates_flag$base_pair_location,
                     reads[[1]]$candidates_flag$base_pair_location)
  }
})

test_that("native 0.4 stores are the default and support full, regional, key, and row reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(2000L)
  path <- tempfile("pcodec-native-")
  store <- compress_sumstats(input, path, overwrite = TRUE)
  expect_equal(store$manifest$format_version, "0.4.6-pcodec-native")
  expect_equal(store$manifest$codec$name, "pcodec_native_standalone_z10_eaf8_se8_xse_zstd_exceptions")
  expect_equal(store$manifest$key_block_rows, 131072L)
  expect_equal(store$manifest$value_block_rows, 65536L)
  expect_equal(store$manifest$codec$page_rows, 131072L)
  expect_equal(store$manifest$tolerances$eaf_abs_max, 0.004)
  expect_true(is.finite(store$manifest$tolerances$se_relative_max))
  expect_false(file.exists(file.path(path, "variants.parquet")))
  expect_true(validate_compressor(store, full = TRUE)$valid)

  full <- read_sumstats(store)
  full_parallel <- read_sumstats(store, threads = 2L)
  expect_equal(nrow(full), nrow(input))
  expect_equal(full_parallel, full)
  expect_false("row" %in% names(full))
  expect_false("variant_id" %in% names(full))
  expect_equal(full[c("chromosome", "base_pair_location", "effect_allele", "other_allele")],
               input[c("chromosome", "base_pair_location", "effect_allele", "other_allele")])
  expect_equal(full$effect_allele_frequency, input$effect_allele_frequency, tolerance = 0.006)
  expect_lt(max(abs(full$beta - input$beta), na.rm = TRUE), 0.02)
  expect_lt(max(abs(full$standard_error - input$standard_error), na.rm = TRUE), 0.01)

  regional <- read_sumstats(store, region = "chr1:100100-100150",
                            columns = c("chromosome", "base_pair_location", "z", "beta",
                                         "standard_error", "effect_allele_frequency"))
  expect_true(all(regional$base_pair_location >= 100100L &
                  regional$base_pair_location <= 100150L))
  key <- compressor_variant_key(input$chromosome[17L], input$base_pair_location[17L],
                                input$other_allele[17L], input$effect_allele[17L])
  by_key <- read_sumstats(store, variants = key,
                          columns = c("chromosome", "base_pair_location", "effect_allele",
                                      "other_allele", "z"))
  expect_equal(nrow(by_key), 1L)
  expect_equal(compressor_variant_key(by_key$chromosome, by_key$base_pair_location,
                                      by_key$other_allele, by_key$effect_allele), key)
  by_row <- read_sumstats(store, variants = 16L, columns = c("z", "standard_error"))
  expect_equal(nrow(by_row), 1L)
  expect_false("row" %in% names(by_row))
  projected <- read_sumstats(store, columns = "z")
  expect_identical(names(projected), "z")
  expect_equal(nrow(projected), nrow(input))
})

test_that("native Pcodec preserves exceptional values and batched reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(100L)
  input$beta[c(1, 2)] <- c(100000, -100)
  input$standard_error[c(1, 2)] <- c(1000, 1e-6)
  input$effect_allele_frequency[3] <- NA_real_
  input$z <- input$beta / input$standard_error
  path1 <- tempfile("pcodec-native-batch-a-")
  path2 <- tempfile("pcodec-native-batch-b-")
  first_store <- compress_sumstats(input, path1, overwrite = TRUE)
  second_store <- compress_sumstats(input, path2, overwrite = TRUE)
  expect_identical(first_store$manifest$integrity$payload_sha256,
                   second_store$manifest$integrity$payload_sha256)
  expect_identical(first_store$manifest$integrity$canonical_sha256,
                   second_store$manifest$integrity$canonical_sha256)
  expect_identical(
    CompreSSoR:::pcodec_canonical_manifest_sha256(first_store$manifest),
    first_store$manifest$integrity$canonical_sha256
  )
  expect_true(validate_compressor(path1, full = TRUE)$valid)
  keys <- compressor_variant_key(input$chromosome, input$base_pair_location,
                                 input$other_allele, input$effect_allele)
  observed <- read_sumstats_batch(c(path1, path2), list(keys[c(1, 3, 100)], keys[2]),
                                  columns = c("z", "standard_error", "effect_allele_frequency"),
                                  threads = 2L)
  serial <- read_sumstats_batch(c(path1, path2), list(keys[c(1, 3, 100)], keys[2]),
                                columns = c("z", "standard_error", "effect_allele_frequency"),
                                threads = 1L)
  expect_length(observed, 2L)
  expect_equal(observed, serial)
  expect_equal(nrow(observed[[1L]]), 3L)
  expect_equal(nrow(observed[[2L]]), 1L)
  expect_equal(observed[[1L]]$z[1], input$z[1], tolerance = 1e-5)
  expect_equal(observed[[1L]]$effect_allele_frequency[2], input$effect_allele_frequency[3])
})

test_that("a native Pcodec store can be used as an identity-only panel", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(128L)
  panel_input <- input[c(1L, 17L, 99L), , drop = FALSE]
  panel_path <- tempfile("pcodec-variant-panel-")
  compress_sumstats(panel_input, panel_path, overwrite = TRUE)

  panel <- CompreSSoR:::read_variant_set(panel_path)
  expect_equal(nrow(panel), 3L)
  expect_true(all(grepl("^[0-9XY]+:[0-9]+:[ACGT]:[ACGT]$", panel$variant_id)))
  expect_equal(attr(panel, "variant_set_metadata")$rows, 3L)

  output_path <- tempfile("pcodec-filtered-")
  filtered_store <- compress_sumstats(
    input, output_path, selection = "core", variant_set = panel_path,
    overwrite = TRUE
  )
  expect_equal(filtered_store$manifest$n_rows, 3L)
  expect_equal(nrow(read_sumstats(filtered_store)), 3L)
})

test_that("native Pcodec keeps configurable stream frames aligned", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(70000L)
  path <- tempfile("pcodec-native-frame-")
  store <- compress_sumstats(input, path, block_rows = 65536L,
                             threads = 4L,
                             overwrite = TRUE)
  expect_equal(store$manifest$block_rows, 65536L)
  expect_equal(store$manifest$writer$requested_workers, 4L)
  expect_true(store$manifest$writer$effective_workers >= 1L)
  expect_lte(store$manifest$writer$effective_workers, 2L)
  expect_identical(store$manifest$writer$partition,
                   "deterministic_contiguous_blocks")
  expect_true(validate_compressor(store, full = TRUE)$valid)
  observed <- read_sumstats(store, columns = c("z", "standard_error",
                                                "effect_allele_frequency"))
  expect_equal(nrow(observed), nrow(input))
  expect_lt(max(abs(observed$z - input$beta / input$standard_error), na.rm = TRUE), 0.02)

  keys <- compressor_variant_key(
    input$chromosome[seq(1L, nrow(input), by = 5000L)],
    input$base_pair_location[seq(1L, nrow(input), by = 5000L)],
    input$other_allele[seq(1L, nrow(input), by = 5000L)],
    input$effect_allele[seq(1L, nrow(input), by = 5000L)]
  )
  serial <- read_sumstats(store, variants = keys,
                          columns = c("chromosome", "base_pair_location", "z",
                                      "standard_error", "effect_allele_frequency"),
                          threads = 1L)
  parallel <- read_sumstats(store, variants = keys,
                            columns = c("chromosome", "base_pair_location", "z",
                                        "standard_error", "effect_allele_frequency"),
                            threads = 2L)
  expect_equal(parallel, serial)
})

test_that("native exception frames use Zstandard and round trip", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  raw <- as.raw(rep(c(0L, 1L, 2L, 3L, 255L), 2000L))
  encoded <- CompreSSoR:::pcodec_native_zstd_compress(raw, level = 19L)
  expect_lt(length(encoded), length(raw))
  expect_identical(CompreSSoR:::pcodec_native_zstd_decompress(encoded, length(raw)), raw)
})

test_that("native writer validates and bounds worker counts", {
  expect_equal(CompreSSoR:::pcodec_native_validate_worker_count(1), 1L)
  for (value in list(0, 1.5, NA_real_, Inf, c(1, 2))) {
    expect_error(
      CompreSSoR:::pcodec_native_validate_worker_count(value),
      "positive integer"
    )
  }
  check_limited <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) &&
    tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_")) != "false"
  expected_workers <- if (check_limited) 2L else 3L
  expect_equal(CompreSSoR:::pcodec_native_effective_workers(8L, 3L), expected_workers)
  expect_equal(CompreSSoR:::pcodec_native_effective_workers(8L, 0L), 0L)
  expect_error(
    CompreSSoR:::pcodec_native_writer_workers(list(threads = 0L)),
    "positive integer"
  )
})

test_that("native stream code validation preserves reserved missing codes", {
  codes <- list(
    z = c(0L, 510L, 511L), se = c(0L, 62L, 63L), eaf = c(0L, 255L),
    position = c(0, 1), substitution = c(1L, 2L)
  )
  expect_silent(CompreSSoR:::pcodec_native_validate_code_domains(
    codes, c("z", "se", "eaf", "position", "substitution"),
    list(z_count = 510L, se_count = 62L, eaf_count = 255L)
  ))
  codes$se[1] <- 64L
  expect_error(
    CompreSSoR:::pcodec_native_validate_code_domains(
      codes, "se", list(se_count = 62L)
    ),
    "out-of-domain"
  )
})

test_that("native index validation rejects non-contiguous partitions", {
  expect_error(
    CompreSSoR:::pcodec_native_validate_index_partition(
      list(list(row_start = 1, row_stop = 2, values = 1)), 2,
      "native test blocks"
    ),
    "contiguous"
  )
})

test_that("native stream block compression merges deterministically", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  values <- as.integer(seq_len(7000L) %% 65536L)
  serial_path <- tempfile("pcodec-native-serial-")
  parallel_path <- tempfile("pcodec-native-parallel-")
  serial <- CompreSSoR:::pcodec_native_append_stream(
    values, serial_path, "u16", block_rows = 1024L, workers = 1L
  )
  parallel <- CompreSSoR:::pcodec_native_append_stream(
    values, parallel_path, "u16", block_rows = 1024L, workers = 4L
  )
  expect_identical(readBin(serial_path, raw(), file.info(serial_path)$size),
                   readBin(parallel_path, raw(), file.info(parallel_path)$size))
  expect_identical(serial$blocks, parallel$blocks)
  expect_equal(serial$effective_workers, 1L)
  check_limited <- nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_")) &&
    tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_")) != "false"
  expected_workers <- if (check_limited) 2L else 4L
  expect_equal(parallel$effective_workers, expected_workers)
})

test_that("native position gaps reset at block boundaries and reject unsorted input", {
  expect_identical(
    CompreSSoR:::pcodec_native_position_gaps(c(0, 4294967295, 4294967295), 2L),
    c(0, 4294967295, 0)
  )
  expect_error(
    CompreSSoR:::pcodec_native_position_gaps(c(2, 1), 2L),
    "sorted"
  )
})

test_that("native empty integer streams are valid empty round trips", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  encoded <- CompreSSoR:::pcodec_native_compress(integer(), "u8")
  expect_length(encoded, 0L)
  expect_identical(CompreSSoR:::pcodec_native_decompress(encoded, 0L, "u8"),
                   integer())
})
