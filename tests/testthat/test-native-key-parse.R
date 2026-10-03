r_target_keys <- function(variants, build = "GRCh38", trim = FALSE) {
  if (trim) variants <- unique(trimws(variants))
  parsed <- CompreSSoR:::pcodec_native_target_keys_r(variants, build = build)
  codes <- sort(unique(as.numeric(parsed$position) * 16 + as.numeric(parsed$substitution)))
  list(position = floor(codes / 16), substitution = as.integer(codes %% 16), codes = codes)
}

r_error <- function(variants, build = "GRCh38", trim = FALSE) {
  tryCatch({ r_target_keys(variants, build, trim); NA_character_ },
           error = function(e) conditionMessage(e))
}

test_that("native key parsing equals the R parser on canonical keys", {
  skip_if_not(is.loaded("compressor_parse_variant_keys", PACKAGE = "CompreSSoR"))
  parse <- CompreSSoR:::pcodec_native_target_keys
  set.seed(3)
  n <- 20000L
  bases <- c("A", "C", "G", "T")
  ref <- sample(4L, n, TRUE)
  alt <- (ref + sample(3L, n, TRUE) - 1L) %% 4L + 1L
  keys <- paste(sample(c(1:22, "X", "Y"), n, TRUE), sample(4e7, n, TRUE),
                bases[ref], bases[alt], sep = ":")
  keys <- c(keys, keys[1:50])                     # duplicates
  for (build in c("GRCh38", "GRCh37")) {
    expect_identical(parse(keys, build), r_target_keys(keys, build))
    expect_identical(parse(keys, build, trim = TRUE), r_target_keys(keys, build, TRUE))
  }
  odd <- c(" chr1:100:a:g\t", "CHRX:5:C:T", "23:5:C:T", "x:6:c:t", "24:10:A:C",
           "Y:11:A:C", "chr22:50818468:G:A", "1:0100:A:G", "1:248956422:T:A")
  expect_identical(parse(odd, trim = TRUE), r_target_keys(odd, trim = TRUE))
  expect_identical(parse(odd[-1]), r_target_keys(odd[-1]))
  expect_identical(parse(character()), r_target_keys(character()))
})

test_that("malformed keys keep the R parser's validation and messages", {
  parse <- CompreSSoR:::pcodec_native_target_keys
  good <- "1:100:A:G"
  bad <- list(NA_character_, "", "1:100:A", "1:100:A:G:T", "1:100:A:A", "1:100:A:N",
              "1:100:AC:G", "1:0:A:G", "1:-5:A:G", "22:50818469:A:G", "25:1:A:G",
              "MT:1:A:G", "01:100:A:G", "chrchr1:1:A:G", "1:12345678901:A:G",
              "1:abc:A:G", "chr:1:A:G")
  for (b in bad) for (trim in c(FALSE, TRUE)) {
    v <- c(good, b)
    expected <- r_error(v, trim = trim)
    expect_false(is.na(expected), info = b)
    expect_error(parse(v, trim = trim), expected, fixed = TRUE, info = b)
  }
  # Spellings the strict native parser leaves to the R parser, which accepts
  # them: the result is still the R parser's.
  lenient <- c("1:1e3:A:G", "1: 100:A:G", "1:100.5:A:G", "1:100:A:G:", "chr1 :7:A:C")
  for (b in lenient) {
    expect_identical(parse(c(good, b), trim = TRUE), r_target_keys(c(good, b), trim = TRUE),
                     info = b)
  }
  # GRCh37 chromosome 22 is longer than GRCh38's
  expect_identical(parse("22:50818469:A:G", "GRCh37"),
                   r_target_keys("22:50818469:A:G", "GRCh37"))
})

test_that("read_sumstats reports malformed keys exactly as before", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  store <- compress_sumstats(make_fixture(600L), tempfile("keyparse-"), overwrite = TRUE)
  for (b in c("1:100:A", "1:100001:A:A", "MT:1:A:G", "1:0:A:G", "22:50818469:A:G")) {
    expected <- r_error(c("1:100001:C:A", b), trim = TRUE)
    expect_error(read_sumstats(store, variants = c("1:100001:C:A", b)), expected,
                 fixed = TRUE, info = b)
    expect_error(read_sumstats_batch(list(store, store), c("1:100001:C:A", b), threads = 2L),
                 expected, fixed = TRUE, info = b)
  }
  expect_error(read_sumstats(store, variants = NA_character_), "invalid canonical variant key",
               fixed = TRUE)
  # whitespace, case and the chr prefix are accepted as before
  keys <- c(" chr1:100001:c:a ", "1:100002:G:C", "1:100002:G:C", "1:1:A:C")
  ref <- pcodec_native_read_store_reference(store, variants = keys, threads = 1L)
  for (th in c(1L, 2L)) {
    expect_identical(CompreSSoR:::pcodec_native_read_store(store, variants = keys, threads = th), ref)
  }
})

test_that("merge join matches the per-block R reader for dense and sparse keys", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  input <- make_fixture(300000L)
  store <- compress_sumstats(input, tempfile("mergejoin-"), overwrite = TRUE)
  cols <- c("chromosome", "base_pair_location", "z", "standard_error")
  set.seed(9)
  sets <- list(all = input$variant_id, every_other = input$variant_id[c(TRUE, FALSE)],
               sparse = sample(input$variant_id, 37L),
               edges = input$variant_id[c(1L, 131072L, 131073L, 300000L)],
               with_absent = c(sample(input$variant_id, 500L), "1:5:A:C", "1:200000000:A:C"))
  for (nm in names(sets)) for (th in c(1L, 3L)) {
    expect_identical(
      CompreSSoR:::pcodec_native_read_store(store, variants = sets[[nm]], columns = cols, threads = th),
      pcodec_native_read_store_reference(store, variants = sets[[nm]], columns = cols, threads = th),
      info = nm)
  }
})

test_that("selective reads pick threads from the touched key blocks", {
  auto <- CompreSSoR:::pcodec_native_select_auto_threads
  # 10 key blocks of 100 rows covering positions 0..9999 (column 6 first, 7 last)
  starts <- seq(0, 900, by = 100)
  mats <- list(position = cbind(0, 0, 100, starts, starts + 100, starts * 10, starts * 10 + 999))
  keys <- function(pos) list(position = pos, substitution = rep(1L, length(pos)),
                             codes = pos * 16 + 1)
  expect_identical(auto(mats, keys(c(5, 6)), NULL, NULL, NULL), 1L)
  expect_identical(auto(mats, keys(c(5, 1500)), NULL, NULL, NULL), 1L)
  expect_identical(auto(mats, keys(c(5, 1500, 2500)), NULL, NULL, NULL), 3L)
  expect_identical(auto(mats, keys(seq(0, 9999, by = 50)), NULL, NULL, NULL),
                   CompreSSoR:::PCODEC_NATIVE_SELECT_MAX_THREADS)
  expect_identical(auto(mats, NULL, NULL, 0, 999), 1L)
  expect_identical(auto(mats, NULL, NULL, 0, 3500), 1L)   # regions stay serial
  expect_identical(auto(mats, NULL, c(0L, 99L, 150L), NULL, NULL), 1L)
  expect_identical(auto(mats, NULL, c(0L, 150L, 250L, 350L), NULL, NULL), 4L)
  # an explicit thread count or option is never overridden
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  store <- compress_sumstats(make_fixture(600L), tempfile("autothreads-"), overwrite = TRUE)
  k <- make_fixture(600L)$variant_id[c(1, 300, 600)]
  expect_identical(read_sumstats(store, variants = k),
                   read_sumstats(store, variants = k, threads = 1L))
})
