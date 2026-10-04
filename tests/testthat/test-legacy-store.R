test_that("a Z9/SE6 store written by a27b32d decodes identically", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  source_dir <- test_path("fixtures", "legacy-z9se6-a27b32d.cpr")
  path <- file.path(tempfile("legacy-"), "legacy.cpr")
  dir.create(path, recursive = TRUE)
  file.copy(list.files(source_dir, full.names = TRUE), path)
  expect_true(validate_compressor(path)$valid)
  semantic <- open_compressor(path)$manifest$semantic_codec
  # The fixture predates explicit z_count/eaf_count fields, so the reader's
  # legacy fallbacks are what is being exercised here.
  expect_null(semantic$z_count)
  expect_identical(semantic$name, "z9/eaf8/se6")
  expected <- readRDS(test_path("fixtures", "legacy-z9se6-a27b32d-reads.rds"))
  got <- legacy_store_reads(path)
  expect_identical(names(got), names(expected))
  # The reference reads were made by a27b32d on macOS/arm64. Identity
  # columns, row sets and the p-value domains must match exactly; decoded
  # doubles may differ in the last ulp (the decoder's tables use the
  # platform's libm, and since 0.7.2 every reader shares one native decoder
  # compiled without FMA contraction, whereas a27b32d decoded selective reads
  # in R and full reads with clang's default contraction), so they are
  # compared to 1e-13 here.
  for (name in names(expected)) {
    e <- expected[[name]]
    g <- got[[name]]
    if (is.data.frame(e)) {
      expect_identical(names(g), names(e), label = name)
      expect_identical(attributes(g)[setdiff(names(attributes(g)), "row.names")],
                       attributes(e)[setdiff(names(attributes(e)), "row.names")],
                       label = paste(name, "attributes"))
      for (column in names(e)) {
        if (is.double(e[[column]])) {
          expect_equal(g[[column]], e[[column]], tolerance = 1e-13,
                       label = paste(name, column))
          expect_identical(is.na(g[[column]]), is.na(e[[column]]),
                           label = paste(name, column, "NA pattern"))
        } else {
          expect_identical(g[[column]], e[[column]], label = paste(name, column))
        }
      }
    } else {
      expect_identical(g, e, label = name)
    }
  }
})

test_that("the legacy Z9/SE6 profile can still be written and is byte-stable", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native backend not built")
  path <- tempfile("legacy-write-")
  old <- options(CompreSSoR.native_profile = "z9/eaf8/se6")
  on.exit(options(old), add = TRUE)
  compress_sumstats(legacy_store_fixture_data(), path, qc = "none",
                    pvalue_flag = TRUE, pvalue_order = TRUE,
                    pvalue_order_threshold = 0.05, overwrite = TRUE)
  source_dir <- test_path("fixtures", "legacy-z9se6-a27b32d.cpr")
  for (stream in c("position.pco", "substitution.pco", "z.pco", "eaf.pco",
                   "se.pco", "pvalue_flag.pco", "pvalue_order.pco")) {
    expect_identical(
      unname(tools::md5sum(file.path(path, stream))),
      unname(tools::md5sum(file.path(source_dir, stream))),
      label = stream
    )
  }
  # Since 0.7.2 an exception record that exists only for a missing EAF
  # (flags = 4) stores zero in its unread Z and log2(SE) fields; every other
  # record, and every field a reader uses, is unchanged.
  new_exc <- CompreSSoR:::pcodec_native_read_all_exceptions(open_compressor(path))
  old_exc <- CompreSSoR:::pcodec_native_read_all_exceptions(open_compressor(source_dir))
  expect_identical(new_exc[c("row", "eaf", "flags")], old_exc[c("row", "eaf", "flags")])
  eaf_only <- old_exc$flags == 4L
  expect_true(any(eaf_only))
  expect_identical(new_exc[!eaf_only, ], old_exc[!eaf_only, ])
  expect_true(all(new_exc$z[eaf_only] == 0 & new_exc$log2se[eaf_only] == 0))
  # Same streams as the fixture, so the same decoded values on this platform.
  expect_identical(legacy_store_reads(path)$full,
                   legacy_store_reads(source_dir)$full)
})
