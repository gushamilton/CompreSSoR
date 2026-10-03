# Reads used to check that stores written by an earlier release decode
# identically. Shared by fixtures/make-legacy-z9se6.R (run against the
# release that wrote the fixture) and test-legacy-store.R.
legacy_store_reads <- function(path) {
  strip <- function(x) {
    attr(x, "source_bytes_read") <- NULL
    x
  }
  all_columns <- c("chromosome", "base_pair_location", "effect_allele",
                   "other_allele", "z", "beta", "standard_error",
                   "effect_allele_frequency", "p_value")
  full <- read_sumstats(path, columns = all_columns)
  keys <- paste(full$chromosome, full$base_pair_location, full$other_allele,
                full$effect_allele, sep = ":")
  list(
    full = strip(full),
    full_threads = strip(read_sumstats(path, columns = all_columns, threads = 2L)),
    region = strip(read_sumstats(path, region = "chr1:100201-101400",
                                 columns = all_columns)),
    variants = strip(read_sumstats(path, variants = keys[c(3, 77, 1500, 2600, 3999)],
                                   columns = all_columns)),
    candidates_1e3 = strip(read_candidates(path, 1e-3)),
    candidates_1e5 = strip(read_candidates(path, 1e-5)),
    candidates_region = strip(read_candidates(path, 0.05, region = "chr2:100001-100800")),
    candidates_flag = strip(read_candidates(path, 5e-8, strategy = "pvalue_flag")),
    flag = read_pvalue_flag(path),
    order = read_pvalue_order(path)
  )
}

legacy_store_fixture_data <- function(n = 4000L, seed = 20261003L) {
  set.seed(seed)
  half <- n / 2L
  pos <- rep(seq.int(100001L, length.out = half), 2L)
  ref <- rep(c("C", "G", "T", "A"), length.out = n)
  alt <- rep(c("A", "C", "G", "T"), length.out = n)
  se <- exp(rnorm(n, log(0.02), 0.6))
  z <- rnorm(n, sd = 1.5)
  hits <- sample.int(n, 80L)
  z[hits] <- sample(c(-1, 1), 80L, TRUE) * runif(80L, 3.4, 30)
  eaf <- runif(n, 0.01, 0.99)
  eaf[sample.int(n, 25L)] <- NA_real_
  se[sample.int(n, 10L)] <- se[sample.int(n, 10L)] * 40
  data.frame(
    chromosome = rep(c("1", "2"), each = half), base_pair_location = pos,
    reference_allele = ref, alternate_allele = alt,
    effect_allele = alt, other_allele = ref,
    beta = z * se, standard_error = se, effect_allele_frequency = eaf,
    p_value = 2 * pnorm(-abs(z)), stringsAsFactors = FALSE
  )
}
