# Integrity hardening: allele orientation under qc = "none", parallel worker
# failures, checksum verification, batch identity sharing, out-of-range rows,
# native error handling, crash-safe manifests and small API edge cases.

skip_unless_native <- function() {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
}

hardening_data <- function(n = 2000L, seed = 1L, chromosome = "1") {
  set.seed(seed)
  ref <- rep(c("C", "G", "T", "A"), length.out = n)
  alt <- rep(c("A", "C", "G", "T"), length.out = n)
  se <- exp(rnorm(n, log(0.02), 0.3))
  z <- rnorm(n)
  z[seq(7L, n, by = 97L)] <- 8
  data.frame(
    chromosome = rep(chromosome, n),
    base_pair_location = seq.int(100001L, by = 7L, length.out = n),
    reference_allele = ref, alternate_allele = alt,
    effect_allele = alt, other_allele = ref,
    beta = z * se, standard_error = se,
    effect_allele_frequency = runif(n, 0.05, 0.95),
    stringsAsFactors = FALSE
  )
}

write_hardening_store <- function(data = hardening_data(), ...) {
  path <- tempfile("hardening-", fileext = ".cpr")
  compress_sumstats(data, path, overwrite = TRUE, ...)
  normalizePath(path)
}

copy_store <- function(path) {
  out <- tempfile("hardening-copy-", fileext = ".cpr")
  dir.create(out)
  file.copy(list.files(path, full.names = TRUE, all.files = TRUE, no.. = TRUE), out)
  normalizePath(out)
}

flip_byte <- function(file, at = NULL) {
  bytes <- readBin(file, raw(), n = file.info(file)$size)
  at <- at %||% max(1L, length(bytes) %/% 2L)
  bytes[at] <- xor(bytes[at], as.raw(0x5a))
  writeBin(bytes, file)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------- (5) alleles

test_that("qc = 'none' drops effect = REF rows like compact QC, never flipping beta", {
  skip_unless_native()
  x <- data.frame(chromosome = "1", base_pair_location = c(1000L, 2000L, 3000L),
                  reference_allele = "T", alternate_allele = "A",
                  effect_allele = c("T", "A", "a "), other_allele = c("A", "T", "t"),
                  beta = 0.028, standard_error = 0.01,
                  effect_allele_frequency = 0.3, stringsAsFactors = FALSE)
  none <- compress_sumstats(x, tempfile(), qc = "none", overwrite = TRUE)
  compact <- compress_sumstats(x, tempfile(), qc = "compact", overwrite = TRUE)
  cols <- c("base_pair_location", "effect_allele", "other_allele", "beta")
  expect_identical(read_sumstats(none, columns = cols),
                   read_sumstats(compact, columns = cols))
  got <- read_sumstats(none, columns = cols)
  # Row 1 (effect = REF) is dropped; rows 2-3 (effect = ALT, any case or
  # spacing) are kept with their beta unchanged.
  expect_identical(got$base_pair_location, c(2000L, 3000L))
  expect_true(all(got$beta > 0))
  safety <- none$manifest$preparation$preparation$identity_safety
  expect_equal(safety$counts$orientation_mismatch, 1L)
  expect_equal(none$manifest$dropped_rows$reasons$orientation_mismatch, 1L)
  expect_equal(compact$manifest$dropped_rows$reasons$orientation_mismatch, 1L)
})

test_that("qc = 'none' keeps payloads unchanged for correctly oriented input", {
  skip_unless_native()
  x <- hardening_data(500L)
  a <- compress_sumstats(x, tempfile(), qc = "none", overwrite = TRUE)
  b <- compress_sumstats(x, tempfile(), qc = "compact", overwrite = TRUE)
  expect_identical(a$manifest$integrity$payload_sha256,
                   b$manifest$integrity$payload_sha256)
  expect_equal(a$manifest$preparation$preparation$identity_safety$counts$orientation_mismatch, 0L)
})

# ------------------------------------------------------- (4) parallel workers

test_that("pcodec_parallel_lapply re-raises worker errors with their label", {
  f <- CompreSSoR:::pcodec_parallel_lapply
  for (threads in c(1L, 2L)) {
    expect_error(
      f(1:4, function(i) if (i == 3L) stop("boom ", i) else i, threads = threads,
        labels = paste0("store-", 1:4), what = "test read"),
      "test read failed for 'store-3': boom 3")
    expect_error(f(1:4, function(i) if (i == 3L) stop("plain") else i, threads = threads),
                 "plain")
  }
  expect_identical(f(1:3, function(i) i * 2L, threads = 2L), list(2L, 4L, 6L))
})

test_that("a dead forked worker is an error, never missing data", {
  skip_on_os("windows")
  main <- Sys.getpid()
  expect_error(
    CompreSSoR:::pcodec_parallel_lapply(1:4, function(i) {
      if (i == 2L && Sys.getpid() != main) tools::pskill(Sys.getpid(), 9L)
      i
    }, threads = 2L, labels = letters[1:4]),
    "returned no result")
  # NULL from FUN is indistinguishable from a dead worker and is rejected.
  expect_error(CompreSSoR:::pcodec_parallel_lapply(1:2, function(i) NULL),
               "returned no result")
})

test_that("region and row candidate reads stop when a decode worker dies", {
  skip_unless_native()
  skip_on_os("windows")
  # Two 131,072-row key blocks so the key decode is split across workers.
  n <- 140000L
  data <- hardening_data(n)
  path <- write_hardening_store(data, pvalue_order = FALSE)
  reg <- "chr1:1-248000000"
  ref <- read_candidates(path, 1e-3, region = reg, columns = c("key", "p_value"),
                         threads = 2L)
  expect_gt(nrow(ref), 0L)
  main <- Sys.getpid()
  original <- CompreSSoR:::candidates_stream_reader
  local_mocked_bindings(candidates_stream_reader = function(store) {
    io <- original(store)
    open <- io$open
    io$open <- function(spec) {
      if (Sys.getpid() != main && identical(spec$file, "position.pco")) {
        tools::pskill(Sys.getpid(), 9L)
      }
      open(spec)
    }
    io
  }, .package = "CompreSSoR")
  expect_error(read_candidates(path, 1e-3, region = reg,
                               columns = c("key", "p_value"), threads = 2L),
               "returned no result|failed")
  expect_error(read_candidates(path, 1e-3, columns = c("key", "p_value"), threads = 2L),
               "returned no result|failed")
})

# ---------------------------------------------------------- (2) checksums

test_that("validate_compressor(full = TRUE) detects a corrupted byte in every stream", {
  skip_unless_native()
  base <- write_hardening_store()
  expect_true(validate_compressor(base, full = TRUE)$valid)
  store <- open_compressor(base)
  streams <- setdiff(names(store$manifest$integrity$files), "native.index.json")
  expect_true(all(c("position.pco", "substitution.pco", "z.pco", "eaf.pco",
                    "se.pco", "exceptions.bin") %in% streams))
  for (file in streams) {
    copy <- copy_store(base)
    flip_byte(file.path(copy, file))
    res <- validate_compressor(copy, full = TRUE)
    expect_false(res$valid, label = file)
    expect_true(any(grepl(file, res$errors, fixed = TRUE)), label = file)
  }
  # A truncated stream is reported by size.
  copy <- copy_store(base)
  bytes <- readBin(file.path(copy, "z.pco"), raw(), n = 1e7)
  writeBin(bytes[-length(bytes)], file.path(copy, "z.pco"))
  res <- validate_compressor(copy, full = TRUE)
  expect_false(res$valid)
  expect_true(any(grepl("z.pco: .*bytes on disk", res$errors)))
})

test_that("a corrupted native index is rejected when it is parsed", {
  skip_unless_native()
  base <- write_hardening_store()
  copy <- copy_store(base)
  index <- file.path(copy, "native.index.json")
  text <- readLines(index)
  hit <- grep("first_position", text)[2L]
  digits <- regmatches(text[hit], regexpr("[0-9]+", text[hit]))
  text[hit] <- sub(digits, as.character(as.numeric(digits) + 1), text[hit], fixed = TRUE)
  writeLines(text, index)
  expect_error(read_sumstats(copy, region = "chr1:100000-100200"),
               "index checksum mismatch")
  res <- validate_compressor(copy, full = FALSE)
  expect_false(res$valid)
  expect_match(paste(res$errors, collapse = " "), "index checksum mismatch")
})

test_that("payload_sha256 must agree with the per-file record", {
  skip_unless_native()
  store <- open_compressor(write_hardening_store())
  expect_length(CompreSSoR:::pcodec_integrity_problems(store), 0L)
  store$manifest$integrity$payload_sha256 <- strrep("0", 64)
  expect_match(CompreSSoR:::pcodec_integrity_problems(store), "payload_sha256",
               all = FALSE)
})

test_that("native file hashing equals digest()", {
  skip_unless_native()
  path <- write_hardening_store()
  files <- list.files(path, full.names = TRUE)
  expect_identical(
    CompreSSoR:::pcodec_sha256_files(files, threads = 3L),
    vapply(files, digest::digest, character(1), algo = "sha256", file = TRUE,
           USE.NAMES = FALSE))
  expect_true(is.na(CompreSSoR:::pcodec_sha256_files(file.path(path, "absent"))))
})

# --------------------------------------------------- (3) batch identity sharing

test_that("batch identity sharing refuses a store whose identity streams are corrupt", {
  skip_unless_native()
  data <- hardening_data(3000L)
  a <- write_hardening_store(data)
  b <- copy_store(write_hardening_store(data))
  keys <- read_sumstats(a, columns = c("chromosome", "base_pair_location",
                                       "other_allele", "effect_allele"))
  kv <- paste(keys$chromosome, keys$base_pair_location, keys$other_allele,
              keys$effect_allele, sep = ":")[seq(1, 3000, by = 3)]
  ok <- read_sumstats_batch(list(a, b), variants = kv, threads = 1L)
  expect_identical(nrow(ok[[2]]), length(kv))
  # Half-copied store: truncated position stream, manifest intact.
  bytes <- readBin(file.path(b, "position.pco"), raw(), n = 1e7)
  writeBin(bytes[seq_len(length(bytes) %/% 2L)], file.path(b, "position.pco"))
  old <- options(CompreSSoR.batch_share_panels = TRUE)
  on.exit(options(old), add = TRUE)
  expect_error(read_sumstats_batch(list(a, b), variants = kv, threads = 1L),
               "failed identity verification")
  expect_error(read_candidates_batch(list(a, b), 5e-8, columns = c("key", "p_value")),
               "failed identity verification")
  # A store whose identity bytes flip without a size change is also caught.
  c2 <- copy_store(a)
  flip_byte(file.path(c2, "substitution.pco"))
  expect_error(read_sumstats_batch(list(a, c2), variants = kv, threads = 1L),
               "failed identity verification")
})

test_that("identity verification names the corrupt store wherever it is in the set", {
  skip_unless_native()
  data <- hardening_data(1500L)
  good <- write_hardening_store(data)
  bad <- copy_store(write_hardening_store(data))
  other <- copy_store(good)
  flip_byte(file.path(bad, "position.pco"))
  for (order in list(c(bad, good, other), c(good, bad, other), c(good, other, bad))) {
    stores <- lapply(order, open_compressor)
    rm(list = ls(CompreSSoR:::.pcodec_identity_verified),
       envir = CompreSSoR:::.pcodec_identity_verified)
    expect_error(CompreSSoR:::pcodec_verify_identity_files(stores, threads = 2L),
                 basename(bad), fixed = TRUE)
  }
  ok <- lapply(c(good, other), open_compressor)
  expect_true(CompreSSoR:::pcodec_verify_identity_files(ok, threads = 2L))
  files <- file.path(c(good, good, bad), "position.pco")
  expect_identical(CompreSSoR:::pcodec_files_equal(files[c(1, 1)], files[c(2, 3)], threads = 2L),
                   c(TRUE, FALSE))
})

test_that("identity verification is cached per file stamp", {
  skip_unless_native()
  data <- hardening_data(1000L)
  a <- write_hardening_store(data)
  b <- write_hardening_store(data)
  stores <- lapply(c(a, b), open_compressor)
  CompreSSoR:::pcodec_verify_identity_files(stores)
  env <- CompreSSoR:::.pcodec_identity_verified
  expect_true(all(c(stores[[1]]$path, stores[[2]]$path) %in% ls(env)))
  calls <- 0L
  local_mocked_bindings(pcodec_sha256_files = function(paths, threads = 1L) {
    calls <<- calls + 1L
    character()
  }, .package = "CompreSSoR")
  CompreSSoR:::pcodec_verify_identity_files(stores)
  expect_identical(calls, 0L)
})

# ------------------------------------------------ (13) rows past chromosome end

test_that("rows past a chromosome end warn, are counted, and stop above the limit", {
  skip_unless_native()
  data <- hardening_data(1000L)
  end21 <- CompreSSoR:::compressor_chromosome_lengths("GRCh38")[["1"]]
  few <- data
  few$base_pair_location[1:5] <- end21 + 1:5
  for (qc in c("compact", "none")) {
    expect_warning(store <- compress_sumstats(few, tempfile(), qc = qc, overwrite = TRUE),
                   "5 of 1,000 input rows")
    expect_identical(store$manifest$dropped_rows$past_chromosome_end, 5L)
    expect_identical(store$manifest$n_rows, 995L)
  }
  many <- data
  many$base_pair_location[1:50] <- end21 + 1:50
  for (qc in c("compact", "none")) {
    expect_error(compress_sumstats(many, tempfile(), qc = qc, overwrite = TRUE),
                 "another genome build")
  }
  old <- options(CompreSSoR.max_out_of_range_fraction = 0.1)
  on.exit(options(old), add = TRUE)
  expect_warning(store <- compress_sumstats(many, tempfile(), overwrite = TRUE),
                 "50 of 1,000")
  expect_identical(store$manifest$dropped_rows$past_chromosome_end, 50L)
})

# ------------------------------------------------- (10) native error handling

test_that("native failures surface as R errors repeatedly without corrupting state", {
  skip_unless_native()
  garbage <- as.raw(seq_len(64) %% 256)
  for (i in 1:200) {
    expect_error(CompreSSoR:::pcodec_native_decompress(garbage, 1000L, "u16"))
    expect_error(.Call("compressor_pcodec_compress_blocks", c(1L, -1L), "u8", 2L, 8L,
                       131072L, 1L, PACKAGE = "CompreSSoR"), "out-of-range")
  }
  bad <- tempfile()
  expect_error(.Call("compressor_read_pcodec_native_codes", rep(bad, 6),
                     matrix(0, 0, 7), matrix(0, 0, 6), matrix(0, 0, 6),
                     matrix(0, 0, 6), matrix(0, 0, 6), matrix(0, 0, 4), -1,
                     "z", "zstd", 1L, PACKAGE = "CompreSSoR"))
  values <- c(1L, 2L, 3L)
  blob <- CompreSSoR:::pcodec_native_compress(values, "u8")
  expect_identical(CompreSSoR:::pcodec_native_decompress(blob, 3L, "u8"), values)
})

# ------------------------------------------------- (12) manifests and threads

test_that("store bytes and hashes do not depend on the thread count", {
  skip_unless_native()
  data <- hardening_data(3000L)
  one <- open_compressor(write_hardening_store(data, threads = 1L))
  four <- open_compressor(write_hardening_store(data, threads = 4L))
  expect_identical(one$manifest$integrity$payload_sha256,
                   four$manifest$integrity$payload_sha256)
  expect_identical(one$manifest$integrity$canonical_sha256,
                   four$manifest$integrity$canonical_sha256)
  index <- jsonlite::read_json(file.path(one$path, "native.index.json"))
  expect_null(index$writer)
  expect_null(index$streams$z$requested_workers)
  expect_null(index$exceptions$effective_workers)
})

test_that("the committed manifest is complete and no temporary files remain", {
  skip_unless_native()
  path <- write_hardening_store()
  expect_setequal(list.files(path, all.files = TRUE, no.. = TRUE),
                  c(names(open_compressor(path)$manifest$integrity$files),
                    "manifest.json", "manifest.sha256"))
  expect_true(is.finite(open_compressor(path)$manifest$timings$phases$commit))
  expect_true(validate_compressor(path, full = TRUE)$valid)
})

test_that("overwrite = FALSE is re-checked at commit time", {
  target <- tempfile("commit-race-")
  transaction <- CompreSSoR:::stage_store_output(target, overwrite = FALSE)
  dir.create(target)  # another writer created the destination meanwhile
  expect_error(CompreSSoR:::commit_store_output(transaction), "already exists")
  expect_true(dir.exists(transaction$staging))
  unlink(c(target, transaction$staging), recursive = TRUE)
})

test_that("stores written before this change still validate in full", {
  skip_unless_native()
  fixture <- test_path("fixtures", "legacy-z9se6-a27b32d.cpr")
  skip_if_not(dir.exists(fixture))
  copy <- copy_store(fixture)
  expect_true(validate_compressor(copy, full = TRUE)$valid)
})

# ------------------------------------------------------------- (16) API edges

test_that("regions parse with or without a chr prefix", {
  skip_unless_native()
  path <- write_hardening_store()
  a <- read_sumstats(path, region = "chr1:100001-101000")
  expect_gt(nrow(a), 0L)
  expect_identical(read_sumstats(path, region = "1:100001-101000"), a)
  expect_identical(read_sumstats(path, region = c("1", "100001", "101000")), a)
  expect_error(read_sumstats(path, region = "ch1:100001-101000"),
               "unsupported region chromosome")
})

test_that("non-integer row IDs are rejected instead of truncated", {
  skip_unless_native()
  path <- write_hardening_store()
  expect_error(read_sumstats(path, variants = c(0, 1.5)), "whole numbers")
  expect_error(read_sumstats_batch(list(path, path), variants = c(0, 1.5)),
               "whole numbers")
  expect_identical(nrow(read_sumstats(path, variants = c(0, 2))), 2L)
})
