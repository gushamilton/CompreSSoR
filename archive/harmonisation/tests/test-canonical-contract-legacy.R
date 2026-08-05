test_that("canonical reference is allele-keyed, stable and external", {
  skip_if_not_installed("arrow")
  source <- data.frame(
    chromosome = c("chr2", "1", "1", "1"),
    pos = c(20L, 30L, 10L, 40L),
    ref = c("A", "C", "A", "AT"),
    alt = c("G", "T", "C", "A"),
    rsid = c("rs-shared", "rs-shared", "rs10", "rs-indel"),
    stringsAsFactors = FALSE
  )
  path <- tempfile(fileext = ".parquet")
  descriptor <- build_canonical_reference(source, path, overwrite = TRUE)
  ref <- arrow::read_parquet(path)
  expect_equal(nrow(ref), 4L)
  expect_equal(ref$variant_id, c("1:10:A:C", "1:30:C:T", "1:40:AT:A", "2:20:A:G"))
  expect_equal(ref$variant_index, 0:3)
  expect_equal(ref$rsid[ref$rsid == "rs-shared"], rep("rs-shared", 2))
  expect_equal(descriptor$id, "canonical_grch38_variant_v2")
  expect_equal(descriptor$variants, path)
  expect_true(file.exists(sub("[.]parquet$", ".manifest.json", path)))
  resolved <- resolve_reference(descriptor)
  expect_equal(resolved$metadata$sha256, descriptor$metadata$sha256)
})

test_that("canonical builder reads narrow PVAR and VCF reference sources", {
  skip_if_not_installed("arrow")
  pvar <- tempfile(fileext = ".pvar")
  writeLines(c("#CHROM POS ID REF ALT", "1 10 rs10 A C", "1 20 rs20 G T"), pvar)
  pvar_out <- tempfile(fileext = ".parquet")
  build_canonical_reference(pvar, pvar_out, overwrite = TRUE)
  pvar_ref <- arrow::read_parquet(pvar_out)
  expect_equal(pvar_ref$variant_id, c("1:10:A:C", "1:20:G:T"))

  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.3", "#CHROM POS ID REF ALT QUAL FILTER INFO",
               "1 30 rs30 C G . PASS .", "1 40 rs40 T A . PASS ."), vcf)
  vcf_out <- tempfile(fileext = ".parquet")
  build_canonical_reference(vcf, vcf_out, overwrite = TRUE)
  vcf_ref <- arrow::read_parquet(vcf_out)
  expect_equal(vcf_ref$variant_id, c("1:30:C:G", "1:40:T:A"))
})

test_that("study ingest cannot silently create or choose a GRCh38 reference", {
  input <- make_fixture(5L)
  expect_error(
    compress_sumstats(input, tempfile(), mode = "qc", reference = NULL, overwrite = TRUE),
    "reference is required"
  )
  expect_error(
    harmonise_sumstats(input, reference = data.frame(), input_build = "GRCh37"),
    "liftover"
  )
})

test_that("liftover accepts gzipped UCSC chain files", {
  skip_if_not_installed("rtracklayer")
  chain <- tempfile(fileext = ".chain.gz")
  con <- gzfile(chain, open = "wt")
  writeLines(c(
    "chain 1 chr1 1000 + 0 1000 chr1 1000 + 0 1000 1",
    "1000",
    ""
  ), con)
  close(con)

  input <- data.frame(
    chromosome = "1", base_pair_location = 10L,
    other_allele = "A", effect_allele = "C",
    stringsAsFactors = FALSE
  )
  got <- CompreSSoR:::lift_table_to_grch38(input, "GRCh37", chain)
  expect_equal(got$chromosome, "1")
  expect_equal(got$base_pair_location, 10L)
  expect_equal(got$.compressor_liftover_status, "mapped")
})

test_that("reference alignment is keyed by alleles rather than rsID", {
  ref <- data.frame(chromosome = "1", base_pair_location = 100L,
                    ref = "A", alt = "C", rsid = "rs-any")
  input <- data.frame(chromosome = "1", base_pair_location = 100L,
                      effect_allele = "C", other_allele = "A",
                      beta = .2, standard_error = .1,
                      rsid = "a-completely-different-id")
  got <- harmonise_sumstats(input, ref)
  expect_equal(nrow(got), 1L)
  expect_equal(got$variant_id, "1:100:A:C")
  expect_equal(got$rsid, "rs-any")
})

test_that("canonical keys are build-aware and normalize chromosome aliases", {
  expect_equal(
    compressor_variant_key("chr23", 100, "a", "g", build = "hg19"),
    "X:100:A:G"
  )
  expect_equal(
    compressor_variant_key(c("chr1", "24"), c(1, 2), c("A", "C"),
                           c("G", "T"), build = "GRCh38"),
    c("1:1:A:G", "Y:2:C:T")
  )
  expect_error(
    compressor_variant_key("MT", 1, "A", "G", build = "GRCh38"),
    "chromosome must be one of"
  )
})
