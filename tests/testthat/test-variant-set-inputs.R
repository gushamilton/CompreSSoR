test_that("gzip TSV variant sets do not require R.utils", {
  path <- tempfile(fileext = ".tsv.gz")
  con <- gzfile(path, open = "wt")
  writeLines(c(
    "variant_id\tchromosome\tbase_pair_location\tother_allele\teffect_allele",
    "1:101:A:G\t1\t101\tA\tG",
    "X:202:C:T\tX\t202\tC\tT"
  ), con)
  close(con)

  got <- CompreSSoR:::read_variant_set(path)
  expect_equal(nrow(got), 2L)
  expect_equal(got$variant_id, c("1:101:A:G", "X:202:C:T"))
  expect_equal(got$chromosome, c("1", "X"))
})

test_that("chromosome panel directories select only requested shards", {
  skip_if_not_installed("data.table")
  shard_dir <- tempfile("panel-shards-")
  dir.create(shard_dir)
  data.table::fwrite(data.frame(variant_id = c("1:100:A:C", "1:200:G:T")),
                     file.path(shard_dir, "chr1.tsv.gz"), sep = "\t", compress = "gzip")
  data.table::fwrite(data.frame(variant_id = "2:100:A:G"),
                     file.path(shard_dir, "chr2.tsv.gz"), sep = "\t", compress = "gzip")

  panel <- CompreSSoR:::read_variant_set(shard_dir, chromosomes = "chr1")
  expect_equal(panel$variant_id, c("1:100:A:C", "1:200:G:T"))
  expect_equal(attr(panel, "variant_set_metadata")$chromosomes, "1")
  expect_true(isTRUE(attr(panel, "variant_set_canonical")))

  multi <- CompreSSoR:::read_variant_set(shard_dir, chromosomes = c("1", "2"))
  expect_equal(multi$variant_id, c("1:100:A:C", "1:200:G:T", "2:100:A:G"))
  expect_equal(attr(multi, "variant_set_metadata")$chromosomes, c("1", "2"))
})
