wide_prepared_fixture <- function(n = 20000L, extras = 24L) {
  data <- make_fixture(n)
  for (index in seq_len(extras)) {
    data[[paste0("unused_", index)]] <- paste0("annotation-", index)
  }
  data
}

test_that("native input projection keeps recognized aliases and drops unused columns", {
  input <- wide_prepared_fixture(2000L, extras = 40L)
  projected <- CompreSSoR:::read_sumstats_input(input, project_columns = TRUE)
  metadata <- attr(projected, "input_read_metadata")

  expect_true(all(CompreSSoR:::required_sumstats_columns() %in%
                  names(CompreSSoR:::normalise_sumstats_columns(projected))))
  expect_true(all(c("variant_id", "rsid", "p_value") %in% names(projected)))
  expect_false(any(grepl("^unused_", names(projected))))
  expect_identical(metadata$columns_before, 13L + 40L)
  expect_lt(metadata$columns_read, metadata$columns_before)
  expect_identical(attr(projected, "source_columns"), names(input))
  expect_identical(attr(projected, "source_columns_read"), names(projected))
  expect_true(is.finite(metadata$elapsed_seconds))

  core_projected <- CompreSSoR:::read_sumstats_input(
    input, project_columns = TRUE, core_only = TRUE
  )
  expect_identical(attr(core_projected, "input_read_metadata")$columns_read, 9L)
  expect_false(any(c("variant_id", "rsid", "p_value") %in% names(core_projected)))
})

test_that("projection handles tab, whitespace, and gzip-delimited files", {
  parent <- tempfile("compressor-projection-input-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE, force = TRUE), add = TRUE)

  input <- wide_prepared_fixture(10000L, extras = 12L)
  path <- file.path(parent, "wide.tsv")
  data.table::fwrite(input, path, sep = "\t")
  gzip <- Sys.which("gzip")
  skip_if(!nzchar(gzip), "gzip is required for the compressed-input projection test")
  gz_path <- paste0(path, ".gz")
  expect_identical(system2(gzip, c("-c", shQuote(path)), stdout = gz_path), 0L)

  projected <- CompreSSoR:::read_sumstats_input(
    gz_path, project_columns = TRUE, core_only = TRUE
  )
  metadata <- attr(projected, "input_read_metadata")
  expect_identical(nrow(projected), nrow(input))
  expect_identical(metadata$columns_before, 13L + 12L)
  expect_identical(metadata$columns_read, 9L)
  expect_lt(metadata$columns_read, metadata$columns_before)
  expect_false(any(grepl("^unused_", names(projected))))
  expect_identical(attr(projected, "source_columns"), names(input))
})

test_that("projected native compression preserves provenance and row policy", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- wide_prepared_fixture(32L, extras = 18L)
  input$reference_allele[1L] <- "AT"
  input$alternate_allele[1L] <- "G"
  input$effect_allele[1L] <- "G"
  input$other_allele[1L] <- "AT"

  path <- tempfile("projected-report-")
  store <- compress_sumstats(input, path, row_policy = "report", overwrite = TRUE)
  expect_identical(store$manifest$n_rows, 31L)
  expect_identical(store$manifest$source$columns_before, 13L + 18L)
  expect_identical(store$manifest$source$columns_read, 9L)
  expect_true(store$manifest$source$projected)
  expect_true(all(c("unused_1", "unused_18") %in% store$manifest$source_columns))
  expect_false(any(grepl("^unused_", unlist(store$manifest$source_columns_read))))

  expect_error(
    compress_sumstats(input, tempfile("projected-error-"),
                      row_policy = "error", overwrite = TRUE),
    "explicit REF and ALT"
  )
})

test_that("Parquet keep_extras explicitly retains the full input path", {
  skip_if_not_installed("arrow")
  input <- wide_prepared_fixture(64L, extras = 3L)
  path <- tempfile("projected-parquet-extras-")
  store <- compress_sumstats(input, path, backend = "parquet", profile = "exact",
                             keep_extras = TRUE, overwrite = TRUE)
  got <- decompress_sumstats(store)

  expect_true(all(paste0("unused_", seq_len(3L)) %in% names(got)))
  expect_identical(store$manifest$source$columns_before, 13L + 3L)
  expect_identical(store$manifest$source$columns_read, 13L + 3L)
  expect_false(store$manifest$source$projected)
})

test_that("Parquet default projection also omits unrequested extras", {
  skip_if_not_installed("arrow")
  input <- wide_prepared_fixture(64L, extras = 3L)
  path <- tempfile("projected-parquet-core-")
  store <- compress_sumstats(input, path, backend = "parquet", profile = "exact",
                             keep_extras = FALSE, overwrite = TRUE)

  expect_identical(store$manifest$source$columns_read, 9L)
  expect_true(store$manifest$source$projected)
  expect_false(any(grepl("^unused_", unlist(store$manifest$source_columns_read))))
  expect_false(any(grepl("^unused_", names(decompress_sumstats(store)))))
})

test_that("bounded wide-input projection reports compact timing evidence", {
  input <- wide_prepared_fixture(20000L, extras = 32L)
  started <- unname(proc.time()[["elapsed"]])
  projected <- CompreSSoR:::read_sumstats_input(input, project_columns = TRUE)
  elapsed <- unname(proc.time()[["elapsed"]]) - started
  metadata <- attr(projected, "input_read_metadata")

  expect_identical(nrow(projected), 20000L)
  expect_gt(metadata$columns_before, metadata$columns_read)
  expect_true(is.finite(elapsed) && elapsed >= 0)
  message(sprintf(
    "projection evidence: rows=%d columns_before=%d columns_read=%d elapsed_seconds=%.3f",
    nrow(input), metadata$columns_before, metadata$columns_read, elapsed
  ))
})
