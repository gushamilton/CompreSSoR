test_that("vectorised key/row reads match full read and the reference implementation", {
  skip_if_not(CompreSSoR:::pcodec_native_available(), "native Pcodec backend is not built")
  n <- 400000L
  input <- make_fixture(n)
  store <- compress_sumstats(input, tempfile("pcodec-key-read-"), overwrite = TRUE)
  cols <- c("chromosome", "base_pair_location", "reference_allele", "alternate_allele",
            "z", "standard_error", "effect_allele_frequency", "p_value")
  full <- CompreSSoR:::pcodec_native_read_store(store, columns = cols, threads = 1L)
  expect_equal(nrow(full), n)
  bs <- 131072L
  # 0-based rows incl. block boundaries
  rows <- c(0L, bs - 1L, bs, 2L * bs - 1L, 2L * bs, 3L * bs - 1L, 3L * bs, n - 1L, 5000L)
  keys <- input$variant_id[rows + 1L]
  expect_true(nrow(unique(full[rows + 1L, c("chromosome", "base_pair_location")])) == length(rows))

  absent <- c("1:50:A:C", "1:200000000:A:C", paste0("1:", 100001L, ":A:T"))
  sel_keys <- c(keys, keys[1:3], absent, rev(keys))   # duplicates + absent + shuffled
  for (th in c(1L, 2L)) {
    kr <- CompreSSoR:::pcodec_native_read_store(store, variants = sel_keys, columns = cols, threads = th)
    rr <- CompreSSoR:::pcodec_native_read_store(store, variants = rows, columns = cols, threads = th)
    ref_k <- pcodec_native_read_store_reference(store, variants = sel_keys, columns = cols, threads = th)
    ref_r <- pcodec_native_read_store_reference(store, variants = c(rows, rows[1:2]), columns = cols, threads = th)
    new_r <- CompreSSoR:::pcodec_native_read_store(store, variants = c(rows, rows[1:2]), columns = cols, threads = th)
    expect_identical(kr$row, sort(unique(rows)))
    expect_identical(kr[names(kr) != "row"], rr[names(rr) != "row"])
    expect_identical(kr, rr)
    expect_equal(kr, full[match(kr$row, full$row), ] |> `rownames<-`(NULL), tolerance = 1e-12, ignore_attr = TRUE)
    expect_identical(kr, ref_k)
    expect_identical(new_r, ref_r)
  }
  # default columns and identity-only columns
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = sel_keys, threads = 1L),
    pcodec_native_read_store_reference(store, variants = sel_keys, threads = 1L))
  # non-identity row selection (value-only columns), unsorted
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = rev(rows), columns = c("z", "p_value"), threads = 1L),
    pcodec_native_read_store_reference(store, variants = rev(rows), columns = c("z", "p_value"), threads = 1L))
  # larger random key set, spanning all blocks
  set.seed(1)
  big <- input$variant_id[sample.int(n, 3000L)]
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = big, columns = cols, threads = 2L),
    pcodec_native_read_store_reference(store, variants = big, columns = cols, threads = 2L))
  # region reads unchanged
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, region = "chr1:200000-300000", columns = cols, threads = 1L),
    pcodec_native_read_store_reference(store, region = "chr1:200000-300000", columns = cols, threads = 1L))
  # same position, wrong substitution must not match
  wrong <- sub(":([ACGT])$", ":N", keys[1])
  wrong <- paste0("1:", input$base_pair_location[1], ":", input$alternate_allele[1], ":", input$reference_allele[1])
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = wrong, columns = cols, threads = 1L),
    pcodec_native_read_store_reference(store, variants = wrong, columns = cols, threads = 1L))
  # all absent
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = absent[1:2], columns = cols, threads = 1L),
    pcodec_native_read_store_reference(store, variants = absent[1:2], columns = cols, threads = 1L))
  # empty selection
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = character(), columns = cols, threads = 1L),
    pcodec_native_read_store_reference(store, variants = character(), columns = cols, threads = 1L))
  expect_identical(
    CompreSSoR:::pcodec_native_read_store(store, variants = integer(), columns = cols, threads = 1L),
    pcodec_native_read_store_reference(store, variants = integer(), columns = cols, threads = 1L))
})
