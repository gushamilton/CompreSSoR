test_that("native key columns match the R definition", {
  skip_if_not(is.loaded("compressor_pcodec_key_columns", PACKAGE = "CompreSSoR"),
              "native key column builder is not built")
  offsets <- CompreSSoR:::compressor_chromosome_offsets("GRCh38")
  lengths <- CompreSSoR:::compressor_chromosome_lengths("GRCh38")
  set.seed(3)
  # Sorted positions across every chromosome, both ends of each chromosome,
  # and an unsorted tail so the binary-search fallback is exercised too.
  ends <- c(as.numeric(offsets), as.numeric(offsets) + as.numeric(lengths) - 1)
  position <- c(sort(c(ends, runif(5000, 0, sum(as.numeric(lengths)) - 1))),
                rev(ends), sample(ends))
  position <- floor(position)
  substitution <- sample(setdiff(0:15, c(0L, 5L, 10L, 15L)), length(position), TRUE)
  native <- CompreSSoR:::pcodec_native_key_columns(position, substitution)
  old <- options(CompreSSoR.native_key_columns = FALSE)
  on.exit(options(old), add = TRUE)
  reference <- CompreSSoR:::pcodec_native_key_columns(position, substitution)
  expect_identical(native, reference)
  expect_identical(names(native), CompreSSoR:::PCODEC_NATIVE_KEY_COLUMNS)
  options(old)

  subset <- CompreSSoR:::pcodec_native_key_columns(
    position, substitution, columns = c("other_allele", "chromosome"))
  expect_identical(names(subset), c("chromosome", "other_allele"))
  expect_identical(subset$other_allele, reference$reference_allele)

  # Integer positions, NA and out-of-table codes follow the R fallback.
  odd_position <- c(0L, NA_integer_, 248956421L, 248956422L)
  odd_substitution <- c(1L, NA_integer_, 16L, -1L)
  native <- CompreSSoR:::pcodec_native_key_columns(odd_position, odd_substitution)
  options(CompreSSoR.native_key_columns = FALSE)
  expect_identical(native,
                   CompreSSoR:::pcodec_native_key_columns(odd_position, odd_substitution))
})

test_that("full reads with native key columns equal the R fallback", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(3000L)
  input$chromosome <- rep(c("1", "2", "X"), each = 1000L)
  store <- compress_sumstats(input, tempfile("key-columns-"), overwrite = TRUE)
  columns <- list(NULL,
                  c("chromosome", "base_pair_location", "other_allele",
                    "effect_allele", "beta", "standard_error",
                    "effect_allele_frequency", "p_value"),
                  c("reference_allele", "alternate_allele", "global_position"),
                  c("beta", "standard_error"))
  for (cols in columns) {
    native <- read_sumstats(store, columns = cols)
    old <- options(CompreSSoR.native_key_columns = FALSE)
    reference <- read_sumstats(store, columns = cols)
    options(old)
    expect_identical(native, reference)
    expect_s3_class(native, "data.frame")
    expect_identical(attr(native, "row.names"), seq_len(3000L))
  }
  full <- read_sumstats(store)
  expect_identical(full$chromosome, input$chromosome)
  expect_identical(full$effect_allele, input$effect_allele)
  expect_identical(full$other_allele, input$other_allele)
})

test_that("native code validation agrees with the R checks", {
  skip_if_not(is.loaded("compressor_validate_native_codes", PACKAGE = "CompreSSoR"),
              "native code validation is not built")
  validate <- CompreSSoR:::pcodec_native_validate_code_domains
  semantic <- list(z_count = 510L, se_count = 62L, eaf_count = 255L)
  identity <- list(chromosome_lengths = c(10, 10), chromosome_offsets = c(0, 10))
  streams <- c("z", "se", "eaf", "position", "substitution")
  good <- list(z = c(0L, 510L, 511L), se = c(0L, 62L, 63L), eaf = c(0L, 255L),
               position = c(0, 1, 19), substitution = c(1L, 2L, 14L))
  cases <- list(
    good,
    modifyList(good, list(z = c(0L, 512L))),
    modifyList(good, list(z = c(-1L, 0L))),
    modifyList(good, list(se = c(0L, 64L))),
    modifyList(good, list(eaf = c(256L))),
    modifyList(good, list(position = c(0, 2, 1))),
    modifyList(good, list(position = c(0, NaN))),
    modifyList(good, list(position = c(0, 20))),
    modifyList(good, list(position = c(-1, 0))),
    modifyList(good, list(substitution = c(0L))),
    modifyList(good, list(substitution = c(16L))),
    modifyList(good, list(substitution = c(5L, 1L)))
  )
  run <- function(codes, native) {
    old <- options(CompreSSoR.native_validate_codes = native)
    on.exit(options(old))
    tryCatch({
      validate(codes, streams, semantic, identity)
      "ok"
    }, error = conditionMessage)
  }
  for (codes in cases) {
    expect_identical(run(codes, TRUE), run(codes, FALSE))
  }
  expect_identical(run(good, TRUE), "ok")
  # Streams that were not read are not checked.
  expect_silent(validate(list(z = c(9999L)), "se", semantic, identity))
})
