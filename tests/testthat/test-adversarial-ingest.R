test_that("factor-valued numeric columns retain their represented values", {
  input <- data.frame(
    chromosome = factor(c("1", "2")),
    base_pair_location = factor(c("100000", "200000")),
    effect_allele = factor(c("A", "C")), other_allele = factor(c("G", "T")),
    beta = factor(c("0.25", "-0.5")),
    standard_error = factor(c("0.05", "0.1")),
    effect_allele_frequency = factor(c("0.2", "0.8"))
  )
  got <- import_sumstats(input)
  expect_equal(got$base_pair_location, c(100000L, 200000L))
  expect_equal(got$beta, c(0.25, -0.5))
  expect_equal(got$standard_error, c(0.05, 0.1))
  expect_equal(got$effect_allele_frequency, c(0.2, 0.8))
})

test_that("already-numeric columns retain double precision", {
  beta <- sin(seq_len(100L) / 13) / 5
  se <- 0.02 + (seq_len(100L) %% 17) / 1000
  input <- data.frame(
    chromosome = "1",
    base_pair_location = seq_len(100L),
    effect_allele = "A",
    other_allele = "C",
    beta = beta,
    standard_error = se,
    effect_allele_frequency = 0.2
  )
  got <- import_sumstats(input)
  expect_identical(got$beta, beta)
  expect_identical(got$standard_error, se)
  expect_identical(got$z, beta / se)
})

test_that("malformed delimited rows fail instead of truncating input", {
  header <- "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error"
  lines <- c(header, "1\t100\tC\tA\t0.2\t0.1", "1\tBROKEN", "1\t300\tT\tG\t0.3\t0.1")
  plain <- tempfile(fileext = ".tsv")
  writeLines(lines, plain)
  compressed <- tempfile(fileext = ".tsv.gz")
  con <- gzfile(compressed, "wt")
  writeLines(lines, con)
  close(con)
  expect_error(import_sumstats(plain), "malformed delimited input")
  expect_error(import_sumstats(compressed), "malformed delimited input|could not read")
})

test_that("conflicting aliases and effect statistics fail closed", {
  aliases <- data.frame(
    chromosome = "1", position = 100L, effect_allele = "C", other_allele = "A",
    beta = 0.2, BETA = 9, standard_error = 0.1, check.names = FALSE
  )
  expect_error(import_sumstats(aliases), "conflicting aliases for beta")

  inconsistent <- data.frame(
    chromosome = "1", position = 100L, effect_allele = "C", other_allele = "A",
    beta = 0.2, standard_error = 0.1, z = 99
  )
  expect_error(import_sumstats(inconsistent), "inconsistent")
})

test_that("convert mode still validates supplied numeric values", {
  bad <- data.frame(
    chromosome = "1", position = 100L, effect_allele = "C", other_allele = "A",
    beta = Inf, standard_error = -1, effect_allele_frequency = 1.2
  )
  expect_error(harmonise_sumstats(bad, reference = NULL, mode = "convert"),
               "beta must be finite|standard_error must be positive|effect_allele_frequency")
})

test_that("single-trait GWAS-VCF FORMAT values are imported", {
  path <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait",
    "1\t101\trs101\tA\tC\t.\tPASS\t.\tES:SE:LP:AF:SS:INFO\t0.2:0.1:2:0.3:10000:0.98"
  ), path)
  got <- import_sumstats(path)
  expect_equal(got$beta, 0.2)
  expect_equal(got$standard_error, 0.1)
  expect_equal(got$p_value, 0.01)
  expect_equal(got$effect_allele_frequency, 0.3)
  expect_equal(got$sample_size, 10000)
  expect_equal(got$info, 0.98)
})

test_that("strict harmonisation stops on duplicate canonical targets", {
  input <- make_fixture(1L)
  input <- rbind(input, input)
  reference <- input[1L, c("chromosome", "base_pair_location", "effect_allele",
                           "other_allele", "effect_allele_frequency", "rsid")]
  expect_error(harmonise_sumstats(input, reference, strict = TRUE), "duplicate=2")
})

test_that("reverse-chain liftover reverse-complements alleles", {
  skip_if_not_installed("rtracklayer")
  chain <- tempfile(fileext = ".chain")
  writeLines(c(
    "chain 1 chr1 1000 + 0 1000 chr1 1000 - 0 1000 1",
    "1000",
    ""
  ), chain)
  input <- data.frame(
    chromosome = "1", base_pair_location = 10L,
    effect_allele = "C", other_allele = "A",
    beta = 0.2, standard_error = 0.1,
    effect_allele_frequency = 0.3, variant_id = "source-id"
  )
  got <- liftover_sumstats(input, input_build = "GRCh37", chain = chain)
  expect_equal(got$base_pair_location, 991L)
  expect_equal(got$effect_allele, "G")
  expect_equal(got$other_allele, "T")
  expect_equal(got$variant_id, "1_991_T_G")
  expect_equal(got$source_variant_id, "source-id")
})

test_that("liftover rejects unsupported or misspelled source builds", {
  input <- make_fixture(1L)
  expect_error(
    liftover_sumstats(input, input_build = "GRCh36", chain = tempfile()),
    "GRCh37/hg19 or GRCh38/hg38"
  )
  expect_error(
    liftover_sumstats(input, input_build = "GRCh73", chain = tempfile()),
    "GRCh37/hg19 or GRCh38/hg38"
  )
})

test_that("palindromic alleles remain ambiguous despite divergent EAF", {
  input <- data.frame(
    chromosome = "1", base_pair_location = 100L,
    effect_allele = "A", other_allele = "T",
    beta = 0.4, standard_error = 0.1,
    effect_allele_frequency = 0.9
  )
  reference <- data.frame(
    chromosome = "1", base_pair_location = 100L,
    reference_allele = "A", alternate_allele = "T",
    effect_allele_frequency = 0.1
  )
  got <- harmonise_sumstats(input, reference, strict = FALSE)
  expect_equal(nrow(got), 0L)
  expect_equal(attr(got, "alignment_stats")$ambiguous_rows, 1L)
})

test_that("parallel and serial harmonisation account for unresolved coordinates equally", {
  aligned <- make_fixture(1L)
  missing <- aligned
  missing$chromosome <- NA_character_
  missing$base_pair_location <- NA_integer_
  missing$rsid <- "rs-not-in-reference"
  missing$variant_id <- missing$rsid
  input <- rbind(aligned, missing)
  reference <- aligned[c("chromosome", "base_pair_location", "effect_allele",
                         "other_allele", "effect_allele_frequency", "rsid")]
  serial <- harmonise_sumstats(input, reference, chrom_threads = 1L)
  parallel <- harmonise_sumstats(input, reference, chrom_threads = 2L)
  expect_equal(serial, parallel)
  expect_equal(attr(serial, "alignment_stats")$unmatched_rows,
               attr(parallel, "alignment_stats")$unmatched_rows)
})

test_that("GWAS-VCF FORMAT values override site-level INFO annotations", {
  path <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", "trait", sep = "\t"),
    paste("1", "100", "rs1", "A", "C", ".", "PASS",
          "ES=9;SE=3;AF=0.9", "ES:SE:AF", "0.2:0.1:0.3", sep = "\t")
  ), path)
  got <- import_sumstats(path)
  expect_equal(got$beta, 0.2)
  expect_equal(got$standard_error, 0.1)
  expect_equal(got$effect_allele_frequency, 0.3)
})

test_that("parallel strict harmonisation reports aggregate QC counts", {
  skip_on_os("windows")
  input <- data.frame(
    chromosome = c("1", "2", "2"), base_pair_location = c(100L, 200L, 200L),
    effect_allele = c("C", "G", "G"), other_allele = c("A", "A", "A"),
    beta = c(0.1, 0.2, 0.2), standard_error = 0.1,
    effect_allele_frequency = 0.2
  )
  reference <- data.frame(
    chromosome = "2", base_pair_location = 200L,
    reference_allele = "A", alternate_allele = "G"
  )
  expect_error(
    harmonise_sumstats(input, reference, strict = TRUE, chrom_threads = 2L),
    "reference alignment failed: unmatched=1, incompatible=0, ambiguous=0, duplicate=2"
  )
})

test_that("report mode retains bounded malformed-statistic diagnostics", {
  input <- make_fixture(4L)
  input$beta[2] <- "not-a-number"
  input$standard_error[3] <- 0
  input$effect_allele_frequency[4] <- 1.5
  result <- preflight_sumstats(input, max_examples = 2L)
  expect_equal(result$report$rejection_counts[["malformed_beta"]], 1L)
  expect_equal(result$report$rejection_counts[["invalid_standard_error"]], 1L)
  expect_equal(result$report$rejection_counts[["invalid_effect_allele_frequency"]], 1L)
  expect_equal(result$report$dropped_rows, 3L)
  expect_equal(nrow(result$data), 1L)
  expect_error(preflight_sumstats(input, strict = TRUE), "structural QC rejected")
})

test_that("missing required columns are reported in report mode", {
  input <- data.frame(
    chromosome = "1", position = 100L, beta = .2, standard_error = .1,
    stringsAsFactors = FALSE
  )
  result <- preflight_sumstats(input)
  expect_setequal(result$report$missing_columns,
                  c("effect_allele", "other_allele"))
  expect_equal(result$report$rejection_counts[["missing_effect_allele"]], 1L)
  expect_equal(result$report$rejection_counts[["missing_other_allele"]], 1L)
  expect_error(import_sumstats(input, strict = TRUE), "Missing required")
})

test_that("strict import rejects out-of-range coordinates and unsupported contigs", {
  input <- make_fixture(2L)
  input$base_pair_location[1] <- 248956423L
  input$chromosome[2] <- "chrM"
  report <- attr(import_sumstats(input), "structural_qc_report")
  expect_equal(report$rejection_counts[["coordinate_out_of_range"]], 1L)
  expect_equal(report$rejection_counts[["unsupported_contig"]], 1L)
  expect_error(import_sumstats(input, strict = TRUE), "structural QC rejected")
})

test_that("compressed simple VCF input uses the same structural path", {
  skip_if_not(any(nzchar(Sys.which("bzip2")), nzchar(Sys.which("gzip"))))
  path <- tempfile(fileext = ".vcf.bz2")
  con <- bzfile(path, open = "wt")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "24\t101\trs101\tA\tC\t.\tPASS\tES=0.2;SE=0.1;AF=0.3;LP=2"
  ), con)
  close(con)
  got <- import_sumstats(path)
  expect_equal(got$chromosome, "Y")
  expect_equal(got$beta, .2)
  expect_equal(got$p_value, .01)
})
