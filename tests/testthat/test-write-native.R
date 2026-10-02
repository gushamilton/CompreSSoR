test_that("streamed object SHA-256 equals digest::digest", {
  skip_if_not(is.loaded("compressor_serialize_sha256", PACKAGE = "CompreSSoR"))
  set.seed(7)
  n <- 5000L
  frame <- data.frame(
    chromosome = sample(c(as.character(1:22), "X", NA), n, TRUE),
    base_pair_location = sample.int(1e8, n),
    reference_allele = sample(c("A", "C", "G", "T", "AT"), n, TRUE),
    beta = c(rnorm(n - 3L), NA, NaN, Inf),
    flag = factor(sample(c("x", "y", NA), n, TRUE)),
    stringsAsFactors = FALSE
  )
  attr(frame, "extra") <- list(a = 1, b = "two")
  objects <- list(frame, frame[0, ], letters, 1:10, list(), NULL,
                  c(a = 1.5, b = NA), paste0("é", 1:100))
  big <- data.frame(x = rep(frame$chromosome, 300), y = rnorm(n * 300))
  for (object in c(objects, list(big))) {
    expected <- digest::digest(object, algo = "sha256", serialize = TRUE)
    expect_identical(CompreSSoR:::object_sha256(object), expected)
    expect_identical(CompreSSoR:::object_sha256(object, threads = 2L), expected)
  }
  for (version in 2:3) {
    for (threads in 1:2) {
      expect_identical(
        .Call("compressor_serialize_sha256", frame, version, 14, threads,
              PACKAGE = "CompreSSoR"),
        digest::digest(frame, algo = "sha256", serializeVersion = version)
      )
    }
  }
})

test_that("map_unique_character and map_unique_values equal the direct transform", {
  x <- c(rep(c(" a", "c ", NA, "", ".", "na", "G"), 40), paste0("v", 1:30))
  f <- function(v) {
    out <- toupper(trimws(v))
    out[out %in% c("", ".", "NA", "N/A")] <- NA_character_
    out
  }
  expect_identical(CompreSSoR:::map_unique_character(x, f), f(x))
  expect_identical(CompreSSoR:::map_unique_character(factor(x), f), f(as.character(factor(x))))
  g <- function(v) {
    out <- trimws(as.character(v))
    out[is.na(v) | out %in% c("", ".", "NA", "N/A")] <- NA_character_
    out
  }
  numeric_x <- c(rep(c(1, NA, NaN, 2.5), 30))
  expect_identical(CompreSSoR:::map_unique_values(x, g), g(x))
  expect_identical(CompreSSoR:::map_unique_values(factor(x), g), g(factor(x)))
  expect_identical(CompreSSoR:::map_unique_values(numeric_x, g), g(numeric_x))
  expect_identical(CompreSSoR:::map_unique_values(rep(NA_character_, 100), g),
                   g(rep(NA_character_, 100)))
  expect_identical(CompreSSoR:::row_index_mask(10L, c(2L, 5L, NA, 11L, 0L)),
                   seq_len(10L) %in% c(2L, 5L, NA, 11L, 0L))
  expect_identical(CompreSSoR:::row_index_mask(5L, integer()), logical(5L))
})

test_that("native block compression equals per-block compression", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  set.seed(11)
  specs <- list(
    list(values = sample.int(256L, 5000L, TRUE) - 1L, dtype = "u8"),
    list(values = sample.int(512L, 5000L, TRUE) - 1L, dtype = "u16"),
    list(values = c(0, 4294967295, as.numeric(sample.int(1e6, 4998L))), dtype = "u32")
  )
  for (spec in specs) {
    for (threads in c(1L, 3L)) {
      blobs <- .Call("compressor_pcodec_compress_blocks",
                     if (spec$dtype == "u32") as.numeric(spec$values) else as.integer(spec$values),
                     spec$dtype, 1024L, CompreSSoR:::PCODEC_NATIVE_LEVEL,
                     CompreSSoR:::PCODEC_NATIVE_PAGE_ROWS, threads, PACKAGE = "CompreSSoR")
      starts <- seq.int(1L, length(spec$values), by = 1024L)
      expected <- lapply(starts, function(start) {
        CompreSSoR:::pcodec_native_compress(
          spec$values[start:min(length(spec$values), start + 1023L)], spec$dtype)
      })
      expect_identical(blobs, expected)
    }
  }
  expect_error(.Call("compressor_pcodec_compress_blocks", c(1L, NA), "u8", 1024L, 8L, 0L, 1L,
                     PACKAGE = "CompreSSoR"), "out-of-range")
})

test_that("native exception frames equal the R record layout", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  rows <- sort(sample.int(9000L, 700L)) - 1L
  exceptions <- data.frame(
    row = as.integer(rows), z = c(rnorm(698L) * 10, NA, Inf),
    log2se = rnorm(700L), eaf = c(NA, runif(699L)),
    flags = sample(1:7, 700L, TRUE)
  )
  blocks <- CompreSSoR:::pcodec_native_block_template(as.numeric(seq_len(9000L)), 1024L)
  out <- tempfile("exceptions-")
  dir.create(out)
  written <- CompreSSoR:::pcodec_native_write_exceptions(exceptions, out, blocks, workers = 3L)
  stops <- vapply(blocks, function(b) as.numeric(b$row_stop), numeric(1))
  member <- findInterval(exceptions$row, stops) + 1L
  expected <- unlist(lapply(seq_along(blocks), function(b) {
    raw_blob <- CompreSSoR:::pcodec_native_exception_bytes(exceptions[member == b, , drop = FALSE])
    if (length(raw_blob)) CompreSSoR:::pcodec_native_zstd_compress(raw_blob, 19L) else raw()
  }))
  path <- file.path(out, "exceptions.bin")
  expect_identical(readBin(path, raw(), file.info(path)$size), expected)
  expect_identical(vapply(written$blocks, `[[`, integer(1), "count"),
                   as.integer(tabulate(member, length(blocks))))
  expect_identical(vapply(written$blocks, `[[`, integer(1), "raw_length"),
                   as.integer(tabulate(member, length(blocks)) * 17L))
})
