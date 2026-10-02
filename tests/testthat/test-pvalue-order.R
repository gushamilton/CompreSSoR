test_that("native stores expose deterministic exact p-value order ranks", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(64L)
  input <- input[c(8:1, 9:64), , drop = FALSE]
  input$p_value <- rep(1e-2, nrow(input))
  input$p_value[match(c(1L, 2L, 3L, 4L), c(8:1, 9:64))] <-
    c(1e-8, 1e-8, 1e-10, 1e-4)
  threshold <- 1e-4
  path <- tempfile("pvalue-order-")
  store <- compress_sumstats(
    input, path, overwrite = TRUE, threads = 1L,
    pvalue_order = TRUE, pvalue_order_threshold = threshold
  )

  domain <- store$manifest$domains$pvalue_order
  expect_identical(domain$format, "aligned_exact_rank_v1")
  expect_identical(domain$dtype, "uint32")
  expect_identical(domain$operator, "<=")
  expect_identical(domain$tie_break, "canonical_variant_key")
  expect_identical(domain$source, "supplied")
  expect_identical(domain$source_column_alias, "p_value_alias")
  expect_equal(domain$threshold, threshold)
  expect_identical(domain$hit_rows, 4L)
  expect_true(file.exists(file.path(path, "pvalue_order.pco")))

  identity <- CompreSSoR:::pcodec_native_identity(input)
  native_order <- order(identity$global_position, identity$substitution,
                        method = "radix")
  native_p <- input$p_value[native_order]
  candidates <- which(native_p <= threshold)
  expected <- candidates[order(native_p[candidates], candidates,
                               method = "radix")] - 1L
  expect_identical(read_pvalue_order(store), as.integer(expected))

  ranks <- read_pvalue_order(store, as = "ranks")
  expect_identical(length(ranks), nrow(input))
  expect_identical(ranks[expected + 1L], seq_along(expected))
  expect_true(all(ranks[-(expected + 1L)] == 0L))
  expect_true(validate_compressor(store, full = TRUE)$valid)
})

test_that("exact p-value order records exact-Z fallback provenance", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(32L)
  input$p_value <- NULL
  store <- compress_sumstats(
    input, tempfile("pvalue-order-derived-"), overwrite = TRUE,
    pvalue_order = TRUE, pvalue_order_threshold = 0.05
  )
  domain <- store$manifest$domains$pvalue_order
  expect_identical(domain$source, "derived_from_exact_prepared_z")
  expect_identical(domain$fallback_statistic,
                   "exact_prepared_z_before_lossy_encoding")
  expect_identical(domain$source_column_present, FALSE)

  supplied <- make_fixture(32L)
  supplied$p_value[c(1L, 2L)] <- c(NA_real_, 1e-8)
  fallback <- compress_sumstats(
    supplied, tempfile("pvalue-order-fallback-"), overwrite = TRUE,
    qc = "none", pvalue_order = TRUE, pvalue_order_threshold = 0.05
  )
  fallback_domain <- fallback$manifest$domains$pvalue_order
  expect_identical(fallback_domain$source,
                   "supplied_with_exact_prepared_z_fallback")
  expect_true(fallback_domain$fallback_rows >= 1L)
  expect_true(fallback_domain$supplied_rows >= 1L)
})

test_that("p-value order fails safely outside its exact threshold", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(48L)
  store <- compress_sumstats(
    input, tempfile("pvalue-order-threshold-"), overwrite = TRUE,
    pvalue_order = TRUE, pvalue_order_threshold = 0.05
  )
  expect_error(
    read_pvalue_order(store, threshold = 0.01),
    "available only at threshold"
  )
  expect_warning(
    approximate <- read_pvalue_order(store, threshold = 0.01,
                                     fallback = "reconstructed"),
    "ordering is approximate"
  )
  expect_type(approximate, "integer")
})

test_that("stores without exact order require explicit reconstructed fallback", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(48L)
  store <- compress_sumstats(
    input, tempfile("pvalue-order-absent-"), overwrite = TRUE,
    pvalue_order = FALSE, pvalue_flag_threshold = 0.05
  )
  expect_null(store$manifest$domains$pvalue_order)
  expect_false(file.exists(file.path(store$path, "pvalue_order.pco")))
  expect_error(read_pvalue_order(store), "no exact p-value ordering domain")

  flag_rows <- read_pvalue_flag(store)
  decoded <- CompreSSoR:::pcodec_native_read_store(
    store, variants = flag_rows, columns = "p_value", threads = 1L
  )
  expected <- as.integer(decoded$row[
    order(decoded$p_value, decoded$row, method = "radix")
  ])
  expect_warning(
    approximate <- read_pvalue_order(
      store, fallback = "reconstructed", threads = 1L
    ),
    "ordering is approximate"
  )
  expect_identical(approximate, expected)
})

test_that("p-value order arguments are validated", {
  expect_error(
    compress_sumstats(make_fixture(8L), tempfile("pvalue-order-invalid-"),
                      pvalue_order = NA),
    "pvalue_order"
  )
  expect_error(
    compress_sumstats(make_fixture(8L), tempfile("pvalue-order-invalid-"),
                      pvalue_order_threshold = 2),
    "pvalue_order_threshold"
  )
  skip_if_not_installed("arrow")
  expect_error(
    compress_sumstats(make_fixture(8L), tempfile("pvalue-order-parquet-"),
                      backend = "parquet", profile = "exact",
                      pvalue_order = TRUE),
    "backend='pcodec' only"
  )
})

test_that("default writes the exact order domain at p <= 0.01 and can be disabled", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  input <- make_fixture(64L)
  store <- compress_sumstats(input, tempfile("order-default-"), overwrite = TRUE)
  domain <- store$manifest$domains$pvalue_order
  expect_false(is.null(domain))
  expect_equal(domain$threshold, 0.01)
  expect_true(isTRUE(domain$standard_default))
  expect_true(validate_compressor(store, full = TRUE)$valid)

  old <- compress_sumstats(input, tempfile("order-off-"), overwrite = TRUE,
                           pvalue_order = FALSE)
  expect_null(old$manifest$domains$pvalue_order)
  expect_null(old$manifest$files$pvalue_order)
  expect_false(file.exists(file.path(old$path, "pvalue_order.pco")))
  expect_error(read_pvalue_order(old), "no exact p-value ordering domain")
  expect_true(is.data.frame(read_candidates(old, 0.01)))
  expect_error(read_candidates(old, 0.01, order = "exact"), "no exact")
  expect_error(read_candidates(old, 0.01, strategy = "exact_order"), "no exact")
  expect_true(validate_compressor(old, full = TRUE)$valid)
})

test_that("exact order and membership match brute force on supplied p edge cases", {
  skip_if_not(CompreSSoR:::pcodec_native_available(),
              "native Pcodec backend is not built")
  n <- 200L
  input <- make_fixture(n)
  set.seed(7)
  p <- stats::runif(n, 0, 0.05)
  p[1:6] <- c(0, 5e-324, 1e-310, 0.01, 0.01, 0.01)       # zero, denormals, at threshold
  p[7:8] <- c(0.01 * (1 + 2.3e-16), 0.01 * (1 - 1.2e-16)) # adjacent doubles around it
  p[9:11] <- 0.003                                         # ties
  p[12:13] <- c(NA_real_, NA_real_)                        # missing supplied p
  input$p_value <- p
  threshold <- 0.01
  store <- suppressWarnings(compress_sumstats(
    input, tempfile("order-edge-"), overwrite = TRUE, qc = "none"
  ))
  identity <- CompreSSoR:::pcodec_native_identity(input)
  o <- order(identity$global_position, identity$substitution, method = "radix")
  np <- p[o]
  # NA p falls back to exact prepared Z; reproduce that for the oracle.
  zz <- suppressWarnings(as.numeric(input$z %||% (input$beta / input$standard_error)))
  eff <- np
  miss <- is.na(eff)
  eff[miss] <- 2 * stats::pnorm(-abs(zz[o][miss]))
  hits <- which(!is.na(eff) & eff <= threshold)
  expected <- hits[order(eff[hits], hits, method = "radix")] - 1L
  expect_identical(read_pvalue_order(store), as.integer(expected))
  expect_identical(store$manifest$domains$pvalue_order$source,
                   "supplied_with_exact_prepared_z_fallback")

  cand <- read_candidates(store, threshold, order = "exact",
                          strategy = "exact_order")
  expect_identical(attr(cand, "candidate_strategy"), "exact_order")
  expect_identical(cand$row, as.integer(expected))
  expect_identical(cand$exact_rank, seq_along(expected))
  # the exactly-at-threshold row is a member under the supplied p
  expect_true((which(np == 0.01)[1L] - 1L) %in% cand$row)

  # default auto semantics are unchanged: reconstructed-p membership
  auto <- read_candidates(store, threshold)
  full <- read_sumstats(store, columns = "p_value")
  expect_identical(sort(auto$row), as.integer(which(!is.na(full$p_value) &
                                                   full$p_value <= threshold) - 1L))
  expect_false(identical(attr(auto, "candidate_strategy"), "exact_order"))
  expect_error(read_candidates(store, 0.001, strategy = "exact_order"),
               "equal to the pvalue_order domain threshold")
})
