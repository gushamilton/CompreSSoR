# These tests exercise the maintained 0.2/0.3 Python format. The native 0.4
# format has its own tests in test-pcodec-native.R; keeping the legacy suite
# explicitly on the Python path prevents a native installation from changing
# the expected on-disk contract of these compatibility tests.
withr::local_options(CompreSSoR.native_pcodec = FALSE)

pcodec_test_python <- function() {
  candidates <- c(
    getOption("CompreSSoR.python", NULL),
    Sys.getenv("COMPRESSOR_PYTHON", unset = ""),
    Sys.which("python3"), Sys.which("python")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    resolved <- if (grepl("[/\\\\]", candidate)) candidate else Sys.which(candidate)
    if (!nzchar(resolved) || !file.exists(resolved)) next
    old <- getOption("CompreSSoR.python")
    options(CompreSSoR.python = resolved)
    ok <- tryCatch({
      CompreSSoR:::pcodec_invocation(refresh = TRUE)
      TRUE
    }, error = function(e) FALSE)
    options(CompreSSoR.python = old)
    if (ok) return(resolved)
  }
  NULL
}

test_that("Pcodec is the default self-contained backend", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(2000L)
  path <- file.path(tempdir(), "pcodec-default.cpr")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  expect_equal(store$manifest$backend, "pcodec")
  expect_equal(store$manifest$variant_storage, "self_contained_identity_key")
  expect_equal(store$manifest$format_version, "0.3.0-pcodec")
  expect_equal(store$manifest$key_block_rows, 4096L)
  expect_equal(store$manifest$value_block_rows, 4096L)
  expect_equal(store$manifest$semantic_codec$se_center_block_rows, 65536L)
  expect_false(file.exists(file.path(path, "variants.parquet")))
  expect_true(validate_compressor(store)$valid)
  expect_true(validate_compressor(store, full = TRUE)$valid)

  got <- read_sumstats(store)
  expect_equal(nrow(got), nrow(input))
  expect_false("variant_id" %in% names(got))
  expect_equal(got[c("chromosome", "base_pair_location", "effect_allele", "other_allele")],
               input[c("chromosome", "base_pair_location", "effect_allele", "other_allele")])
  expect_equal(got$effect_allele_frequency, input$effect_allele_frequency, tolerance = 0.006)
  expect_lt(max(abs(got$beta - input$beta), na.rm = TRUE), 0.02)
  expect_lt(max(abs(got$standard_error - input$standard_error), na.rm = TRUE), 0.01)

  old_native <- getOption("CompreSSoR.native_bridge")
  options(CompreSSoR.native_bridge = FALSE)
  fallback <- read_sumstats(store)
  options(CompreSSoR.native_bridge = old_native)
  expect_equal(got, fallback, tolerance = 1e-12)

  regional <- read_sumstats(store, region = "chr1:100100-100150",
                            columns = c("chromosome", "base_pair_location", "z", "beta",
                                         "standard_error", "effect_allele_frequency", "p_value"))
  expected <- input[input$chromosome == "1" & input$base_pair_location >= 100100L &
                      input$base_pair_location <= 100150L, , drop = FALSE]
  expect_equal(nrow(regional), nrow(expected))
  expect_true(all(regional$base_pair_location >= 100100L & regional$base_pair_location <= 100150L))

  selected <- read_sumstats(store, variants = c(0L, nrow(input) - 1L),
                            columns = c("base_pair_location", "z", "beta"))
  expect_equal(nrow(selected), 2L)
  keys <- compressor_variant_key(
    input$chromosome, input$base_pair_location,
    input$other_allele, input$effect_allele
  )
  selected_by_key <- read_sumstats(
    store,
    variants = c(keys[c(17L, 3L)], "2:200000000:A:C"),
    columns = c("chromosome", "base_pair_location", "effect_allele", "other_allele", "z")
  )
  selected_keys <- compressor_variant_key(
    selected_by_key$chromosome, selected_by_key$base_pair_location,
    selected_by_key$other_allele, selected_by_key$effect_allele
  )
  expect_setequal(selected_keys, keys[c(3L, 17L)])
  expect_error(read_sumstats(store, variants = "rs123"),
               "chromosome:position:REF:ALT")
  identity_only <- read_sumstats(store, variants = 0L, columns = "chromosome")
  expect_identical(names(identity_only), "chromosome")
  expect_error(read_sumstats(store, variants = 1.5), "whole-number")
  expect_error(read_sumstats(store, columns = character()), "at least one")
})

test_that("canonical key construction validates the identity contract", {
  expect_identical(
    compressor_variant_key(character(), numeric(), character(), character()),
    character()
  )
  expect_identical(compressor_variant_key("chr1", 123, "a", "g"), "1:123:A:G")
  expect_identical(
    compressor_variant_key(c("1", "X"), c(1, 2), "A", c("C", "T")),
    c("1:1:A:C", "X:2:A:T")
  )
  expect_error(compressor_variant_key("MT", 1, "A", "C"), "1-22")
  expect_error(compressor_variant_key("1", 1.5, "A", "C"), "whole-number")
  expect_error(compressor_variant_key("1", 248956423, "A", "C"), "outside its GRCh38")
  expect_identical(compressor_variant_key("X", 156040895, "A", "G"),
                   "X:156040895:A:G")
  expect_error(compressor_variant_key("1", 1, "A", "A"), "distinct")
})

test_that("batched canonical reads equal independent reads", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(200L)
  path <- tempfile("pcodec-batch-")
  compress_sumstats(input, path, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  keys <- compressor_variant_key(
    input$chromosome, input$base_pair_location,
    input$other_allele, input$effect_allele
  )
  variants <- list(keys[c(2L, 50L, 199L)], keys[17L])
  columns <- c("chromosome", "base_pair_location", "effect_allele",
               "other_allele", "beta", "standard_error")
  observed <- read_sumstats_batch(
    c(first = path, second = path), variants, columns = columns, threads = 2L
  )
  expected <- lapply(variants, function(selected) {
    read_sumstats(path, variants = selected, columns = columns)
  })
  expect_named(observed, c("first", "second"))
  expect_equal(unname(observed), expected, tolerance = 1e-12)
  old_report <- getOption("CompreSSoR.report_source_bytes")
  options(CompreSSoR.report_source_bytes = TRUE)
  reported <- read_sumstats_batch(
    path, variants[1L], columns = columns, threads = 1L
  )
  options(CompreSSoR.report_source_bytes = old_report)
  expect_gt(attr(reported, "source_bytes_read"), 0)
  old_coalesce <- getOption("CompreSSoR.coalesce_batch_reads")
  options(CompreSSoR.coalesce_batch_reads = FALSE)
  uncoalesced <- read_sumstats_batch(
    c(first = path, second = path), rep(variants[1L], 2L),
    columns = columns, threads = 2L
  )
  options(CompreSSoR.coalesce_batch_reads = old_coalesce)
  expect_equal(unname(uncoalesced), rep(expected[1L], 2L), tolerance = 1e-12)

  options(CompreSSoR.coalesce_batch_reads = FALSE)
  large_reply <- read_sumstats_batch(
    rep(path, 100L), rep(variants[2L], 100L),
    columns = columns, threads = 1L
  )
  options(CompreSSoR.coalesce_batch_reads = old_coalesce)
  expect_length(large_reply, 100L)
  expect_equal(unname(large_reply), rep(expected[2L], 100L), tolerance = 1e-12)
  expect_error(read_sumstats_batch(path, keys[1L], threads = 0L),
               "positive integer")
})

test_that("persistent Pcodec worker is reused and restarts safely", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- options(
    CompreSSoR.python = python,
    CompreSSoR.persistent_worker = TRUE
  )
  CompreSSoR:::pcodec_stop_worker()
  on.exit({
    CompreSSoR:::pcodec_stop_worker()
    options(old)
  }, add = TRUE)

  input <- make_fixture(100L)
  path <- tempfile("pcodec-worker-")
  compress_sumstats(input, path, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  columns <- c("chromosome", "base_pair_location", "beta", "standard_error")
  first <- read_sumstats(path, variants = 0:9, columns = columns)
  worker <- get(".pcodec_state", asNamespace("CompreSSoR"))$worker
  expect_true(worker$is_alive())
  first_pid <- worker$get_pid()

  second <- read_sumstats(path, variants = 0:9, columns = columns)
  expect_equal(second, first, tolerance = 1e-12)
  expect_identical(
    get(".pcodec_state", asNamespace("CompreSSoR"))$worker$get_pid(),
    first_pid
  )

  worker$kill()
  restarted <- read_sumstats(path, variants = 0:9, columns = columns)
  expect_equal(restarted, first, tolerance = 1e-12)
  expect_false(identical(
    get(".pcodec_state", asNamespace("CompreSSoR"))$worker$get_pid(),
    first_pid
  ))

  options(CompreSSoR.persistent_worker = FALSE)
  fallback <- read_sumstats(path, variants = 0:9, columns = columns)
  expect_equal(fallback, first, tolerance = 1e-12)
})

test_that("Pcodec convert mode requires an explicit REF/ALT assertion", {
  input <- make_fixture(5L)
  expect_error(
    compress_sumstats(input, tempfile(), reference = NULL, mode = "convert",
                      overwrite = TRUE),
    "assume_grch38_ref_alt"
  )
  expect_error(
    compress_sumstats(input, tempfile(), reference = input, mode = "qc",
                      drop_unresolved = FALSE, overwrite = TRUE),
    "cannot retain unresolved"
  )
})

test_that("Pcodec validation detects corruption in a payload", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)
  input <- make_fixture(1000L)
  path <- file.path(tempdir(), "pcodec-corrupt.cpr")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  z_path <- file.path(path, store$manifest$files$z)
  con <- file(z_path, open = "r+b")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  offset <- floor(file.info(z_path)$size / 2)
  seek(con, where = offset, origin = "start")
  byte <- readBin(con, "raw", n = 1L)
  seek(con, where = offset, origin = "start")
  writeBin(as.raw(bitwXor(as.integer(byte), 1L)), con)
  close(con)
  result <- validate_compressor(path)
  expect_false(result$valid)
  expect_match(result$errors, "checksum")
  expect_error(read_sumstats(path, columns = "z"), "checksum mismatch|corrupt")
})

test_that("Pcodec rejects unsupported exact and extras requests clearly", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)
  input <- make_fixture(20L)
  expect_error(
    compress_sumstats(input, tempfile(), reference = NULL, mode = "convert",
                      profile = "exact", assume_grch38_ref_alt = TRUE,
                      overwrite = TRUE),
    "backend='parquet'"
  )
  expect_error(
    compress_sumstats(input, tempfile(), reference = NULL, mode = "convert",
                      keep_extras = TRUE, assume_grch38_ref_alt = TRUE,
                      overwrite = TRUE),
    "core summary-statistics"
  )
})

test_that("Pcodec drops and audits identities outside its SNV contract", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(3L)
  input$chromosome <- c("1", "1", "MT")
  input$effect_allele <- c("A", "AT", "G")
  input$other_allele <- c("C", "A", "T")
  path <- tempfile("pcodec-supported-identities-")
  store <- compress_sumstats(
    input, path, reference = NULL, mode = "convert",
    assume_grch38_ref_alt = TRUE, overwrite = TRUE
  )
  expect_equal(store$manifest$n_rows, 1L)
  expect_equal(store$manifest$harmonisation$alignment$unsupported_pcodec_rows, 2L)
  expect_equal(store$manifest$harmonisation$alignment$dropped_unsupported_pcodec, 2L)
  expect_equal(read_sumstats(store)$base_pair_location, input$base_pair_location[1L])

  expect_error(
    compress_sumstats(
      input, tempfile(), reference = NULL, mode = "convert", strict = TRUE,
      assume_grch38_ref_alt = TRUE, overwrite = TRUE
    ),
    "unsupported row"
  )
})

test_that("Pcodec handles spaces and failed overwrites atomically", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(100L)
  path <- file.path(tempdir(), "pcodec store with spaces.cpr")
  unlink(path, recursive = TRUE, force = TRUE)
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  expect_true(validate_compressor(store)$valid)
  manifest_before <- readLines(file.path(path, "manifest.json"), warn = FALSE)

  broken <- rbind(input, input[1L, , drop = FALSE])
  expect_error(
    compress_sumstats(broken, path, reference = NULL, mode = "convert",
                      assume_grch38_ref_alt = TRUE, overwrite = TRUE),
    "duplicate full REF/ALT identity keys"
  )
  expect_identical(readLines(file.path(path, "manifest.json"), warn = FALSE),
                   manifest_before)
  expect_true(validate_compressor(path)$valid)
})

test_that("Pcodec rejects altered manifests before native decoding", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(100L)
  path <- tempfile("pcodec-manifest-contract-")
  compress_sumstats(input, path, reference = NULL, mode = "convert",
                    assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  manifest_path <- file.path(path, "manifest.json")
  manifest <- CompreSSoR:::read_manifest(manifest_path)
  manifest$semantic_codec$eaf_bits <- 1L
  CompreSSoR:::write_manifest(manifest, manifest_path)
  CompreSSoR:::seal_pcodec_manifest(manifest_path)

  result <- validate_compressor(path, full = TRUE)
  expect_false(result$valid)
  expect_match(result$errors, "semantic codec constants")

  manifest$semantic_codec$eaf_bits <- 8L
  manifest$semantic_codec$z_range <- c(-35, 35)
  CompreSSoR:::write_manifest(manifest, manifest_path)
  expect_error(open_compressor(path), "manifest checksum mismatch")
})

test_that("projected full reads avoid unrelated durable streams", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(500L)
  path <- tempfile("pcodec-projection-")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  unrelated <- file.path(
    path, unlist(store$manifest$files[c("position", "eaf", "se")], use.names = FALSE)
  )
  hidden <- paste0(unrelated, ".hidden")
  expect_true(all(file.rename(unrelated, hidden)))
  on.exit({
    restore <- file.exists(hidden)
    if (any(restore)) file.rename(hidden[restore], unrelated[restore])
  }, add = TRUE)

  z <- read_sumstats(store, columns = "z")
  expect_identical(names(z), "z")
  expect_equal(nrow(z), nrow(input))
})

test_that("whole, fallback, regional and sparse exception reads agree", {
  python <- pcodec_test_python()
  skip_if(is.null(python), "Pcodec Python dependencies are not available")
  old <- getOption("CompreSSoR.python")
  options(CompreSSoR.python = python)
  on.exit(options(CompreSSoR.python = old), add = TRUE)

  input <- make_fixture(20L)
  input$z <- input$beta / input$standard_error
  input$standard_error[10L] <- 1000
  input$z[10L] <- 100
  input$beta[10L] <- 100000
  path <- tempfile("pcodec-access-parity-")
  store <- compress_sumstats(input, path, reference = NULL, mode = "convert",
                             assume_grch38_ref_alt = TRUE, overwrite = TRUE)
  columns <- c("z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  whole <- read_sumstats(store, columns = columns)
  old_native <- getOption("CompreSSoR.native_bridge")
  options(CompreSSoR.native_bridge = FALSE)
  fallback <- read_sumstats(store, columns = columns)
  options(CompreSSoR.native_bridge = old_native)
  regional <- read_sumstats(
    store,
    region = paste0("chr1:", input$base_pair_location[10L], "-", input$base_pair_location[10L]),
    columns = columns
  )
  sparse <- read_sumstats(store, variants = 9L, columns = columns)
  expected <- whole[10L, , drop = FALSE]
  row.names(expected) <- NULL

  expect_equal(whole, fallback, tolerance = 1e-12)
  expect_equal(regional, expected, tolerance = 1e-12)
  expect_equal(sparse, expected, tolerance = 1e-12)
})

test_that("native bridge fails closed on hostile codec domains", {
  bridge <- tempfile("pcodec-hostile-bridge-")
  dir.create(bridge)
  writeBin(as.raw(0L), file.path(bridge, "se_code.bin"))
  writeBin(as.raw(255L), file.path(bridge, "eaf_code.bin"))
  file.create(file.path(bridge, "exception_row.bin"),
              file.path(bridge, "exception_z.bin"),
              file.path(bridge, "exception_se.bin"),
              file.path(bridge, "exception_eaf.bin"),
              file.path(bridge, "exception_flags.bin"))
  hostile_call <- function() .Call(
    "compressor_read_pcodec_bridge", normalizePath(bridge), "standard_error",
    1, 0, -3.5, 3.5, 510L, 62L, 1L, 65536L, numeric(),
    PACKAGE = "CompreSSoR"
  )
  for (i in seq_len(20L)) {
    expect_error(hostile_call(), "codec constants are invalid")
  }
  malformed_count <- function(value) .Call(
    "compressor_read_pcodec_bridge", normalizePath(bridge), "standard_error",
    value, 0, -3.5, 3.5, 510L, 62L, 255L, 65536L, numeric(),
    PACKAGE = "CompreSSoR"
  )
  expect_error(malformed_count(NA_real_), "row count is invalid")
  expect_error(malformed_count(Inf), "row count is invalid")
  expect_error(malformed_count(1.5), "row count is invalid")
  expect_error(.Call(
    "compressor_read_pcodec_bridge", normalizePath(bridge), NA_character_,
    1, 0, -3.5, 3.5, 510L, 62L, 255L, 65536L, numeric(),
    PACKAGE = "CompreSSoR"
  ), "requested columns contain NA")
})
