test_that("reference descriptors download, cache and enter the manifest", {
  skip_if_not_installed("arrow")
  source_path <- tempfile(fileext = ".tsv")
  writeLines(c(
    "variant_id\tchromosome\tbase_pair_location\teffect_allele\tother_allele",
    "v1\t1\t10\tA\tC",
    "v2\t1\t20\tG\tA"
  ), source_path)
  spec <- list(
    id = "fixture-reference",
    build = "GRCh38",
    source = "test fixture",
    variants = list(
      url = paste0("file://", normalizePath(source_path)),
      filename = "fixture-reference.tsv"
    )
  )
  cache_dir <- tempfile("compressor-reference-cache-")
  resolved <- download_reference(spec, cache_dir = cache_dir)
  expect_true(file.exists(resolved$variants))
  expect_equal(resolved$metadata$id, "fixture-reference")
  expect_true(nzchar(resolved$metadata$sha256))
  resolved_again <- resolve_reference(resolved, cache_dir = cache_dir)
  expect_equal(normalizePath(resolved_again$variants), normalizePath(resolved$variants))

  input <- data.frame(
    variant_id = c("v1", "v2"), chromosome = c("1", "1"),
    base_pair_location = c(10L, 20L), effect_allele = c("A", "G"),
    other_allele = c("C", "A"), beta = c(0.1, -0.2),
    standard_error = c(0.1, 0.2), stringsAsFactors = FALSE
  )
  store <- compress_sumstats(input, tempfile("reference-store-"),
                             reference = resolved, backend = "parquet", overwrite = TRUE)
  expect_equal(store$manifest$reference$id, "fixture-reference")
  expect_equal(store$manifest$reference$rows, 2L)
  expect_true(file.exists(store$manifest$reference$normalized_cache_path))
})

test_that("the built-in reference is an external, pinned descriptor", {
  reference <- grch38_reference()
  expect_equal(reference$build, "GRCh38")
  expect_equal(reference$id, "ebi_ensembl95_grch38_all_v1")
  expect_null(reference$variants)
  expect_error(resolve_reference(reference), "build_ebi_reference")
})
