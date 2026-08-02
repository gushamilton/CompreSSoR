test_that("standard stores round-trip and validate", {
  input <- make_fixture()
  reference <- input[c("chromosome", "base_pair_location", "effect_allele", "other_allele", "variant_id", "effect_allele_frequency")]
  path <- file.path(tempdir(), "standard.cpr")
  store <- compress_sumstats(input, path, reference = reference, overwrite = TRUE)
  expect_s3_class(store, "compressor_store")
  expect_true(validate_compressor(store)$valid)
  expect_equal(store$manifest$benchmark$benchmark_id, "first_pareto_10m_full_read")
  expect_equal(store$manifest$benchmark_comparisons$vcf_tabix$benchmark_id,
               "vcf_tabix_11m")
  expect_equal(store$manifest$benchmark_comparisons$finngen_end_to_end$benchmark_id,
               "finngen_r2_antidepressants_end_to_end")
  expect_equal(store$manifest$benchmark_comparisons$finngen_optimization$benchmark_id,
               "finngen_r2_antidepressants_reference_cache_optimization")
  expect_equal(store$manifest$benchmark_comparisons$modes_edge$benchmark_id,
               "compressor_modes_edge_550k")
  expect_equal(store$manifest$benchmark_comparisons$release_gate$benchmark_id,
               "compressor_release_gate_finngen_16111549")
  expect_equal(store$manifest$benchmark_comparisons$release_gate_followup$benchmark_id,
               "compressor_release_gate_slice_parallel_100k")
  got <- read_sumstats(store, region = "chr1:100100-100250",
                       columns = c("variant_id", "beta", "standard_error", "p_value", "annotation"))
  expect_equal(nrow(got), 151L)
  expect_true(all(c("beta", "standard_error", "p_value", "annotation") %in% names(got)))
  expect_equal(got$variant_id[1], "1_100100")
})

test_that("the shipped benchmark table is available", {
  bench <- benchmark_table()
  expect_equal(nrow(bench), 11L)
  expect_true(all(c("full_read_median_seconds",
                    "compression_ratio_vs_source_gzip",
                    "pareto_frontier") %in% names(bench)))
  expect_equal(bench$full_read_median_seconds[bench$format == "q9 Parquet"], 0.0159)
  vcf <- benchmark_table("vcf")
  expect_equal(vcf$benchmark_id, "vcf_tabix_11m")
  expect_equal(vcf$full_read_median_seconds, 13.28)
  expect_equal(vcf$regional_rows, 40822L)
  finngen <- benchmark_table("finngen")
  expect_equal(finngen$output_rows, 3294606L)
  expect_equal(finngen$validated, TRUE)
  optimization <- benchmark_table("finngen_optimization")
  expect_equal(nrow(optimization), 7L)
  expect_equal(median(optimization$elapsed_seconds[optimization$cache_state == "warm"]), 63.820)
  modes <- benchmark_table("modes")
  expect_equal(nrow(modes), 9L)
  expect_equal(median(modes$elapsed_seconds[modes$scenario == "qc_serial"]), 0.374)
  release <- benchmark_table("release_gate")
  expect_true(nrow(release) >= 60L)
  expect_equal(median(release$elapsed_seconds[release$scenario == "finngen_convert_standard_compress"]), 46.006)
  followup <- benchmark_table("release_gate_followup")
  expect_true(nrow(followup) >= 45L)
  expect_equal(median(followup$elapsed_seconds[followup$scenario == "slice_qc_chrom4_compress"]), 18.985)
})

test_that("exact profile retains numeric values", {
  input <- make_fixture(100L)
  path <- file.path(tempdir(), "exact.cpr")
  store <- compress_sumstats(input, path, reference = NULL, profile = "exact", overwrite = TRUE)
  got <- read_sumstats(store)
  expect_equal(got$beta, input$beta, tolerance = 0)
  expect_equal(got$p_value, input$p_value, tolerance = 0)
})

test_that("reverse allele alignment flips beta and EAF", {
  input <- make_fixture(20L)
  reference <- input[c("chromosome", "base_pair_location", "effect_allele", "other_allele", "variant_id", "effect_allele_frequency")]
  input$effect_allele[1] <- reference$other_allele[1]
  input$other_allele[1] <- reference$effect_allele[1]
  original_beta <- input$beta[1]
  original_eaf <- input$effect_allele_frequency[1]
  path <- file.path(tempdir(), "flip.cpr")
  store <- compress_sumstats(input, path, reference = reference, profile = "exact", overwrite = TRUE)
  got <- read_sumstats(store)
  expect_equal(got$beta[1], -original_beta)
  expect_equal(got$effect_allele_frequency[1], 1 - original_eaf)
})

test_that("strand complements are aligned without changing the effect", {
  input <- data.frame(
    chromosome = "1", base_pair_location = 200L,
    effect_allele = "T", other_allele = "G", beta = 0.2,
    standard_error = 0.1, effect_allele_frequency = 0.2,
    p_value = 0.05, variant_id = "v1", stringsAsFactors = FALSE
  )
  reference <- data.frame(
    chromosome = "1", base_pair_location = 200L,
    effect_allele = "A", other_allele = "C", variant_id = "v1",
    effect_allele_frequency = 0.2, stringsAsFactors = FALSE
  )
  path <- file.path(tempdir(), "complement.cpr")
  store <- compress_sumstats(input, path, reference = reference,
                             profile = "exact", overwrite = TRUE)
  got <- read_sumstats(store)
  expect_equal(got$effect_allele, "A")
  expect_equal(got$other_allele, "C")
  expect_equal(got$beta, 0.2)
})

test_that("harmonisation is available as an explicit step", {
  input <- make_fixture(20L)
  reference <- input[c("chromosome", "base_pair_location", "effect_allele",
                       "other_allele", "variant_id", "effect_allele_frequency")]
  got <- harmonise_sumstats(input, reference)
  expect_equal(nrow(got), nrow(input))
  expect_equal(attr(got, "genome_build"), "GRCh38")
  expect_true(is.character(attr(got, "reference_hash")))
})

test_that("palindromic ambiguity fails closed", {
  input <- data.frame(
    chromosome = "1", base_pair_location = 300L,
    effect_allele = "A", other_allele = "T", beta = 0.2,
    standard_error = 0.1, effect_allele_frequency = 0.5,
    p_value = 0.05, variant_id = "v2", stringsAsFactors = FALSE
  )
  reference <- input[c("chromosome", "base_pair_location", "effect_allele", "other_allele", "variant_id")]
  expect_error(
    compress_sumstats(input, file.path(tempdir(), "ambiguous.cpr"),
                      reference = reference, profile = "exact", strict = TRUE,
                      overwrite = TRUE),
    "ambiguous"
  )
})

test_that("effect-allele frequency is optional", {
  input <- make_fixture(10L)
  input$effect_allele_frequency <- NULL
  path <- file.path(tempdir(), "no-eaf.cpr")
  store <- compress_sumstats(input, path, reference = NULL, profile = "exact", overwrite = TRUE)
  got <- read_sumstats(store)
  expect_true(all(is.na(got$effect_allele_frequency)))
})

test_that("gzipped delimited files and common aliases are accepted", {
  input <- make_fixture(25L)
  aliased <- input[c("chromosome", "base_pair_location", "effect_allele",
                     "other_allele", "beta", "standard_error",
                     "effect_allele_frequency", "p_value", "variant_id")]
  names(aliased) <- c("chr", "pos", "a1", "a2", "b", "se", "eaf", "pval", "SNP")
  input_path <- tempfile(fileext = ".tsv.gz")
  con <- gzfile(input_path, open = "wt")
  write.table(aliased, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  store <- compress_sumstats(input_path, file.path(tempdir(), "aliases.cpr"),
                             reference = NULL, profile = "exact", overwrite = TRUE)
  got <- read_sumstats(store)
  expect_equal(nrow(got), 25L)
  expect_equal(got$variant_id, input$variant_id)
})
