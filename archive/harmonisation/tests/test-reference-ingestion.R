test_that("PVAR metadata and rows are read across supported transports", {
  lines <- c(
    "##fileformat=VCFv4.3",
    "##contig=<ID=1,length=100>",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">",
    "#CHROM POS ID REF ALT FILTER INFO",
    "chr1 10 rs10 A C PASS AF=0.2",
    "1 20 rs20 G T PASS AF=0.3"
  )
  plain <- tempfile(fileext = ".pvar")
  writeLines(lines, plain)
  paths <- plain

  gzip <- tempfile(fileext = ".pvar.gz")
  con <- gzfile(gzip, open = "wt")
  writeLines(lines, con)
  close(con)
  paths <- c(paths, gzip)

  zstd <- unname(Sys.which("zstd"))
  if (nzchar(zstd)) {
    zstandard <- tempfile(fileext = ".pvar.zst")
    status <- system2(zstd, c("-q", "-f", plain, "-o", zstandard))
    if (identical(as.integer(status), 0L)) paths <- c(paths, zstandard)
  }

  for (path in paths) {
    got <- CompreSSoR:::read_reference_source(path)
    expect_equal(nrow(got), 2L)
    expect_equal(got$chromosome, c("chr1", "1"))
    metadata <- attr(got, "reference_source_metadata")
    expect_equal(metadata$format, CompreSSoR:::reference_source_format(path))
    expect_equal(metadata$columns[1:5], c("CHROM", "POS", "ID", "REF", "ALT"))
    expect_equal(metadata$metadata$contig[[1]]$length, "100")
  }
})

test_that("reference harmonisation keeps provenance and bounded diagnostics", {
  reference <- data.frame(
    chromosome = c("1", "1"), base_pair_location = c(100L, 200L),
    reference_allele = c("A", "A"), alternate_allele = c("C", "G"),
    effect_allele_frequency = c(0.2, 0.3), stringsAsFactors = FALSE
  )
  input <- data.frame(
    chromosome = c("1", "1", "1"), base_pair_location = c(100L, 200L, 300L),
    effect_allele = c("C", "A", "T"), other_allele = c("A", "G", "C"),
    beta = c(0.2, 0.3, 0.1), standard_error = 0.1,
    effect_allele_frequency = c(0.2, 0.9, 0.4), stringsAsFactors = FALSE
  )
  got <- harmonise_sumstats(
    input, reference, drop_unresolved = FALSE,
    qc = list(frequency = "report", example_limit = 1L)
  )
  expect_equal(attr(got, "source_provenance")$kind, "data.frame")
  expect_equal(attr(got, "alignment_stats")$diagnostic_counts$unmatched_rows, 1L)
  expect_lte(length(attr(got, "alignment_stats")$diagnostic_examples$unmatched), 1L)
  expect_equal(attr(got, "audit")$qc$frequency, "report")
  expect_equal(got$frequency_qc_status, c("match", "mismatch", "not_checked"))
})

test_that("strand, palindromic and liftover controls are explicit", {
  reference <- data.frame(
    chromosome = "1", base_pair_location = 100L,
    reference_allele = "A", alternate_allele = "C",
    stringsAsFactors = FALSE
  )
  input <- data.frame(
    chromosome = "1", base_pair_location = 100L,
    effect_allele = "G", other_allele = "T", beta = 0.2,
    standard_error = 0.1, effect_allele_frequency = 0.2,
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(harmonise_sumstats(input, reference, strand = FALSE)), 0L)
  expect_equal(nrow(harmonise_sumstats(input, reference, strand = TRUE)), 1L)

  palindromic_reference <- transform(reference,
    base_pair_location = 200L, alternate_allele = "T",
    effect_allele_frequency = 0.1
  )
  palindromic_input <- transform(input,
    base_pair_location = 200L, effect_allele = "A", other_allele = "T",
    effect_allele_frequency = 0.9
  )
  dropped <- harmonise_sumstats(palindromic_input, palindromic_reference)
  recovered <- harmonise_sumstats(palindromic_input, palindromic_reference,
                                  palindromic = "frequency")
  expect_equal(nrow(dropped), 0L)
  expect_equal(nrow(recovered), 1L)
  expect_true(recovered$harmonisation_flip)
})

test_that("regional and sparse Parquet reference queries retain global indices", {
  skip_if_not_installed("arrow")
  source <- data.frame(
    chromosome = c("1", "1", "1", "2"),
    base_pair_location = c(10L, 20L, 30L, 40L),
    reference_allele = c("A", "G", "C", "T"),
    alternate_allele = c("C", "T", "A", "C"),
    rsid = paste0("rs", 1:4), stringsAsFactors = FALSE
  )
  output <- tempfile(fileext = ".parquet")
  descriptor <- build_canonical_reference(source, output, overwrite = TRUE)
  regional <- CompreSSoR:::reference_table(descriptor, region = "chr1:15-25")
  sparse <- CompreSSoR:::reference_table(descriptor, variant_ids = "1:30:C:A")
  expect_equal(regional$variant_id, "1:20:G:T")
  expect_equal(regional$variant_index, 1L)
  expect_equal(sparse$variant_id, "1:30:C:A")
  expect_equal(sparse$variant_index, 2L)
})
