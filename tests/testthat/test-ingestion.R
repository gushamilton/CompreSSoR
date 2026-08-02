test_that("fruity delimited aliases are canonicalised without losing rsids", {
  input <- data.frame(
    `#CHROM` = c("chr1", "2"),
    POS = c(101L, 202L),
    A1 = c("A", "G"), A2 = c("C", "A"),
    BETA = c("0.2", "-0.3"), SE = c("0.1", "0.2"),
    EAF = c("0.2", "0.7"), P = c("0.05", "0.1"),
    rsids = c("rs101", "rs202"), annotation = c("x", "y"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  got <- harmonise_sumstats(input, reference = NULL, mode = "convert")
  expect_equal(got$chromosome, c("1", "2"))
  expect_equal(got$variant_id, c("rs101", "rs202"))
  expect_equal(got$rsid, c("rs101", "rs202"))
  expect_equal(got$beta, c(0.2, -0.3))
  expect_equal(got$annotation, c("x", "y"))
  expect_equal(attr(got, "genome_build"), "unknown")
})

test_that("delimited files accept comma, whitespace and gzip transport", {
  input <- make_fixture(10L)
  comma_path <- tempfile(fileext = ".csv")
  write.csv(input, comma_path, row.names = FALSE)
  whitespace_path <- tempfile(fileext = ".txt")
  write.table(input, whitespace_path, row.names = FALSE, quote = FALSE)
  gzip_path <- tempfile(fileext = ".tsv.gz")
  con <- gzfile(gzip_path, open = "wt")
  write.table(input, con, sep = "\t", row.names = FALSE, quote = FALSE)
  close(con)
  for (path in c(comma_path, whitespace_path, gzip_path)) {
    store <- compress_sumstats(path, tempfile(), reference = NULL,
                               mode = "convert", profile = "exact", overwrite = TRUE)
    expect_true(validate_compressor(store)$valid)
    expect_equal(nrow(decompress_sumstats(store)), 10L)
  }
})

test_that("VCF and VCF.gz INFO fields are converted safely", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=ES,Number=1,Type=Float,Description=Effect>",
    "##INFO=<ID=SE,Number=1,Type=Float,Description=Standard error>",
    "##INFO=<ID=AF,Number=1,Type=Float,Description=Frequency>",
    "##INFO=<ID=LP,Number=1,Type=Float,Description=-log10 p>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t101\trs101\tC\tA\t.\tPASS\tES=0.2;SE=0.1;AF=0.2;LP=2",
    "2\t202\t.\tA\tG\t.\tPASS\tOR=1.5;SE=0.2;AF=0.7;P=0.1"
  )
  vcf_path <- tempfile(fileext = ".vcf")
  writeLines(lines, vcf_path)
  vcf_gz <- tempfile(fileext = ".vcf.gz")
  con <- gzfile(vcf_gz, open = "wt")
  writeLines(lines, con)
  close(con)
  for (path in c(vcf_path, vcf_gz)) {
    got <- harmonise_sumstats(path, reference = NULL, mode = "convert")
    expect_equal(nrow(got), 2L)
    expect_equal(got$variant_id, c("rs101", "2_202_A_G"))
    expect_equal(got$beta[1], 0.2)
    expect_equal(got$beta[2], log(1.5), tolerance = 1e-12)
    expect_equal(got$p_value, c(0.01, 0.1), tolerance = 1e-12)
    expect_equal(got$effect_allele_frequency, c(0.2, 0.7), tolerance = 1e-12)
  }
})

test_that("multiallelic VCF records fail with an actionable message", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t101\trs101\tC\tA,G\t.\tPASS\tBETA=0.2;SE=0.1"
  )
  path <- tempfile(fileext = ".vcf")
  writeLines(lines, path)
  expect_error(harmonise_sumstats(path, reference = NULL, mode = "convert"),
               "multiallelic")
})

test_that("a BIM reference can anchor sumstats", {
  input <- make_fixture(3L)
  input$other_allele[2] <- "T"
  path <- tempfile(fileext = ".bim")
  writeLines(c(
    "1 1_100001 0 100001 C A",
    "1 1_100002 0 100002 T C",
    "1 1_100003 0 100003 T G"
  ), path)
  got <- harmonise_sumstats(input, reference = path, mode = "qc")
  expect_equal(nrow(got), 3L)
  expect_true(all(got$harmonisation_status == "aligned"))
})
