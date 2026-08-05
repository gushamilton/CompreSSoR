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

test_that("native 0.4 stores are the default and support full, regional, key, and row reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(2000L)
  path <- tempfile("pcodec-native-")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  expect_equal(store$manifest$format_version, "0.4.5-pcodec-native")
  expect_equal(store$manifest$codec$name, "pcodec_native_standalone_z9_eaf8_se6_zstd_exceptions")
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
  compress_sumstats(input, path1, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  compress_sumstats(input, path2, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
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
  compress_sumstats(panel_input, panel_path, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)

  panel <- CompreSSoR:::read_variant_set(panel_path)
  expect_equal(nrow(panel), 3L)
  expect_true(all(grepl("^[0-9XY]+:[0-9]+:[ACGT]:[ACGT]$", panel$variant_id)))
  expect_equal(attr(panel, "variant_set_metadata")$rows, 3L)

  filtered <- harmonise_sumstats(input, reference = NULL, mode = "convert",
                                 variant_set = panel_path)
  expect_equal(nrow(filtered), 3L)
  expect_equal(filtered$base_pair_location, panel$base_pair_location)

  output_path <- tempfile("pcodec-filtered-")
  filtered_store <- compress_sumstats(
    input, output_path, reference = NULL, mode = "convert",
    variant_set = panel_path, assume_grch38_ref_alt = TRUE, overwrite = TRUE
  )
  expect_equal(filtered_store$manifest$n_rows, 3L)
  expect_equal(nrow(read_sumstats(filtered_store)), 3L)
})

test_that("native Pcodec keeps configurable stream frames aligned", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(70000L)
  path <- tempfile("pcodec-native-frame-")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, block_rows = 65536L,
                             chrom_threads = 4L,
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
