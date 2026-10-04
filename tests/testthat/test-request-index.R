test_that("compressor_identity_code matches the codes stores return", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  input <- make_fixture(3000L)
  store <- compress_sumstats(input, tempfile("identity-code-"), overwrite = TRUE)
  keys <- input$variant_id[c(1L, 17L, 400L, 2999L)]
  got <- read_sumstats_batch(store$path, keys,
                             columns = c("global_position", "substitution"))[[1L]]
  expect_identical(compressor_identity_code(keys),
                   as.numeric(got$global_position) * 16 + as.numeric(got$substitution))
  # Reader-equivalent spellings give the same code; invalid keys give NA.
  expect_identical(compressor_identity_code(c(" chr1:100001:C:A", "1:100001:C:A")),
                   rep(compressor_identity_code("1:100001:C:A"), 2L))
  expect_identical(compressor_identity_code("23:1000:A:G"), compressor_identity_code("X:1000:A:G"))
  invalid <- c("1:100001:C:C", "1:100001:CA:A", "1:100001:c:a", "1:0:C:A",
               "1:300000000:C:A", "MT:5:A:G", "1:100001:C", NA, "")
  expect_true(all(is.na(compressor_identity_code(invalid))))
  expect_identical(compressor_identity_code(character()), numeric())
  expect_false(identical(compressor_identity_code("2:1000:A:G", "GRCh37"),
                         compressor_identity_code("2:1000:A:G", "GRCh38")))
  expect_error(compressor_identity_code(1:3), "character vector")
  expect_true(all(c("key_reads_identity_columns", "request_index", "identity_code") %in%
                    compressor_capabilities()))
})

test_that("read_sumstats_batch(request_index = TRUE) indexes each row's request", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  input <- make_fixture(3000L)
  other <- make_fixture(3000L)
  other$base_pair_location <- other$base_pair_location + 7L
  other$variant_id <- paste0("1:", other$base_pair_location, ":", other$reference_allele, ":",
                             other$alternate_allele)
  stores <- c(
    a = compress_sumstats(input, tempfile("request-index-a-"), overwrite = TRUE)$path,
    b = compress_sumstats(transform(input, beta = beta * 2), tempfile("request-index-b-"),
                          overwrite = TRUE)$path,
    c = compress_sumstats(other, tempfile("request-index-c-"), overwrite = TRUE)$path
  )
  keys <- c(input$variant_id[c(10L, 2L, 500L, 2L)], "1:50:A:C", other$variant_id[3L],
            paste0(" ", input$variant_id[900L]))
  rows <- c(5L, 0L, 5L, 2999L, 4000L, 12L)
  columns <- c("beta", "standard_error", "chromosome", "base_pair_location",
               "reference_allele", "alternate_allele")
  requests <- list(keys, keys, rows)
  previous <- options(CompreSSoR.batch_share_panels = NULL)
  on.exit(options(previous), add = TRUE)
  for (share in c(FALSE, TRUE)) {
    options(CompreSSoR.batch_share_panels = share)
    for (threads in c(1L, 2L)) {
      plain <- read_sumstats_batch(stores, requests, columns = columns, threads = threads)
      indexed <- read_sumstats_batch(stores, requests, columns = columns, threads = threads,
                                     request_index = TRUE)
      expect_named(indexed, names(stores))
      for (i in seq_along(stores)) {
        out <- indexed[[i]]
        expect_identical(names(out), c(columns, "request_index"))
        expected <- plain[[i]]
        expected$request_index <- out$request_index
        expect_identical(out, expected, info = paste(share, threads, i))
        expect_type(out$request_index, "integer")
        expect_false(anyNA(out$request_index))
      }
      # Keys: each row answers the first request with its key.
      for (i in 1:2) {
        out <- indexed[[i]]
        row_keys <- paste(out$chromosome, out$base_pair_location, out$reference_allele,
                          out$alternate_allele, sep = ":")
        expect_identical(trimws(keys[out$request_index]), row_keys)
        expect_identical(out$request_index, match(row_keys, trimws(keys)))
      }
      # Row IDs: each row answers the first request with its row ID.
      expect_identical(rows[indexed[[3L]]$request_index], sort(unique(rows[rows < 3000L])))
      expect_identical(indexed[[3L]]$request_index, c(2L, 1L, 6L, 4L))
    }
  }
  options(previous)
  # Region-only and whole-store reads get NA; identity columns requested by
  # the caller are kept.
  region <- read_sumstats_batch(stores[1:2], NULL, columns = c("global_position", "beta"),
                                region = "chr1:100001-100010", request_index = TRUE)
  expect_identical(names(region[[1L]]), c("global_position", "beta", "request_index"))
  expect_true(all(is.na(region[[1L]]$request_index)))
  expect_identical(nrow(region[[1L]]), 10L)
  keyed <- read_sumstats_batch(stores[[1L]], keys,
                               columns = c("substitution", "beta"), request_index = TRUE)[[1L]]
  expect_identical(names(keyed), c("substitution", "beta", "request_index"))
  # Empty results keep the column.
  none <- read_sumstats_batch(stores[[1L]], "1:50:A:C", columns = "beta",
                              request_index = TRUE)[[1L]]
  expect_identical(none$request_index, integer())
  expect_error(read_sumstats_batch(stores[[1L]], rows, columns = "beta",
                                   region = "chr1:100001-100010", request_index = TRUE),
               "row-ID requests combined with a region")
  expect_error(read_sumstats_batch(stores[[1L]], keys, request_index = NA), "TRUE or FALSE")
  expect_error(read_sumstats_batch(stores[[1L]], keys, columns = c("beta", "request_index"),
                                   request_index = TRUE), "cannot also be requested")
})
