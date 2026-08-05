test_that("canonical builder keeps sequence variants and splits multiallelic VCFs", {
  skip_if_not_installed("arrow")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.3",
    "#CHROM POS ID REF ALT QUAL FILTER INFO",
    "1 10 rs10 A C,G . PASS .",
    "1 20 rs20 AT A . PASS .",
    "1 30 rs30 C TTT . PASS ."
  ), vcf)
  output <- tempfile(fileext = ".parquet")
  descriptor <- build_canonical_reference(vcf, output, overwrite = TRUE)
  reference <- arrow::read_parquet(output)
  expect_equal(reference$variant_id,
               c("1:10:A:C", "1:10:A:G", "1:20:AT:A", "1:30:C:TTT"))
  expect_equal(reference$variant_type, c("SNV", "SNV", "INDEL", "INDEL"))
  expect_equal(descriptor$metadata$rows, 4L)
})

test_that("EBI builder writes one globally indexed partitioned dictionary", {
  skip_if_not_installed("arrow")
  source_dir <- tempfile("ebi-vcf-")
  dir.create(source_dir)
  write_vcf_gz <- function(path, lines) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con)
  }
  write_vcf_gz(file.path(source_dir, "homo_sapiens-chr1.vcf.gz"), c(
    "##fileformat=VCFv4.3",
    "#CHROM POS ID REF ALT QUAL FILTER INFO",
    "1 10 rs10 A C . PASS .",
    "1 20 rs20 AT A . PASS .",
    "1 20 rs20_duplicate AT A . PASS ."
  ))
  write_vcf_gz(file.path(source_dir, "homo_sapiens-chr2.vcf.gz"), c(
    "##fileformat=VCFv4.3",
    "#CHROM POS ID REF ALT QUAL FILTER INFO",
    "2 30 rs30 C T . PASS ."
  ))

  output <- file.path(tempdir(), paste0("compressor-ebi-", as.integer(Sys.time())))
  descriptor <- build_ebi_reference(output, source = source_dir,
                                    chromosomes = c("1", "2"), chunk_rows = 1L)
  expect_true(dir.exists(output))
  expect_true(file.exists(file.path(output, "manifest.json")))
  expect_true(file.exists(file.path(output, "variants.parquet")))
  expect_false(dir.exists(file.path(output, "aliases")))
  expect_equal(descriptor$id, "ebi_ensembl95_grch38_all_v2")
  expect_equal(descriptor$metadata$rows, 3L)

  master <- arrow::read_parquet(file.path(output, "variants.parquet"))
  expect_setequal(names(master), c(
    "chromosome_code", "position_block", "position_delta", "reference_allele",
    "alternate_allele", "rsid_code", "rsid_other"
  ))
  expect_equal(master$position_delta, c(10L, 10L, 30L))
  expect_equal(master$position_block, c(0L, 0L, 1L))
  chr1 <- reference_table(descriptor, chromosome = "1")
  chr2 <- reference_table(descriptor, chromosome = "2")
  expect_equal(chr1$variant_index, 0:1)
  expect_equal(chr2$variant_index, 2L)

  input <- data.frame(
    chromosome = c("1", "2"),
    base_pair_location = c(10L, 30L),
    effect_allele = c("C", "T"),
    other_allele = c("A", "C"),
    beta = c(.1, -.2),
    standard_error = c(.1, .2),
    stringsAsFactors = FALSE
  )
  got <- harmonise_sumstats(input, descriptor)
  expect_equal(got$variant_id, c("1:10:A:C", "2:30:C:T"))

  rsid_input <- data.frame(
    effect_allele = "C", other_allele = "A", rsid = "rs10",
    variant_id = "rs10", beta = .1, standard_error = .1,
    stringsAsFactors = FALSE
  )
  rsid_got <- harmonise_sumstats(rsid_input, descriptor)
  expect_equal(rsid_got$variant_id, "1:10:A:C")
  secondary <- CompreSSoR:::reference_alias_table(descriptor, "rs20_duplicate")
  expect_equal(secondary$variant_id, "1:20:AT:A")

  corrupt <- tempfile("ebi-corrupt-reference-")
  dir.create(corrupt)
  file.copy(file.path(output, "variants.parquet"),
            file.path(corrupt, "variants.parquet"))
  corrupt_manifest <- descriptor$metadata
  corrupt_manifest$local_path <- corrupt
  corrupt_manifest$sha256 <- paste(rep("0", 64L), collapse = "")
  CompreSSoR:::write_manifest(corrupt_manifest, file.path(corrupt, "manifest.json"))
  expect_error(resolve_reference(corrupt), "failed its declared SHA-256")
})

test_that("truncated gzip reference inputs fail integrity validation", {
  source <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(source, open = "wt")
  writeLines(c(
    "##fileformat=VCFv4.3",
    "#CHROM POS ID REF ALT QUAL FILTER INFO",
    sprintf("1 %d rs%d A C . PASS .", seq_len(1000L), seq_len(1000L))
  ), con)
  close(con)
  bytes <- readBin(source, "raw", n = file.info(source)$size)
  truncated <- tempfile(fileext = ".vcf.gz")
  writeBin(bytes[seq_len(length(bytes) - 20L)], truncated)
  expect_error(
    CompreSSoR:::ebi_vcf_stream(truncated, chunk_lines = 100L,
                                callback = function(data) invisible(data)),
    "gzip integrity validation"
  )
})

test_that("rsIDs resolve coordinates through the canonical alias table", {
  skip_if_not_installed("arrow")
  reference <- data.frame(
    chromosome = "1", base_pair_location = 100L,
    reference_allele = "A", alternate_allele = "C",
    rsid = "rs100", stringsAsFactors = FALSE
  )
  input <- data.frame(
    effect_allele = "C", other_allele = "A",
    rsid = "rs100", variant_id = "rs100",
    beta = .2, standard_error = .1,
    stringsAsFactors = FALSE
  )
  got <- harmonise_sumstats(input, reference)
  expect_equal(nrow(got), 1L)
  expect_equal(got$variant_id, "1:100:A:C")
  expect_equal(got$chromosome, "1")
  expect_equal(got$base_pair_location, 100L)
  expect_equal(attr(got, "alignment_stats")$alias_resolution$resolved, 1L)
})
