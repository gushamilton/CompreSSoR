# Write-path speedups must not change any result: each faster routine is
# checked against the a27b32d implementation in helper-write-reference.R.

qc_adversarial_frame <- function(n = 400L, seed = 1L) {
  set.seed(seed)
  bases <- c("A", "C", "G", "T")
  chromosome <- sort(sample(c(1:22), n, TRUE))
  position <- integer(n)
  for (chr in unique(chromosome)) {
    rows <- which(chromosome == chr)
    position[rows] <- sort(sample.int(1e6, length(rows)))
  }
  ref <- sample(bases, n, TRUE)
  alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1L), "")
  data <- data.frame(
    chromosome = as.character(chromosome), base_pair_location = position,
    reference_allele = ref, alternate_allele = alt,
    effect_allele = alt, other_allele = ref,
    beta = rnorm(n, 0, 0.1), standard_error = runif(n, 0.01, 0.2),
    effect_allele_frequency = runif(n), p_value = runif(n),
    stringsAsFactors = FALSE
  )
  edit <- function(row, ...) {
    values <- list(...)
    for (name in names(values)) data[row, name] <<- values[[name]]
  }
  row <- 0L
  nxt <- function() { row <<- row + 3L; row }
  edit(nxt(), chromosome = NA)
  edit(nxt(), chromosome = "")
  edit(nxt(), chromosome = "MT")
  edit(nxt(), chromosome = "chrUn_gl000220")
  edit(nxt(), chromosome = "chr7")
  edit(nxt(), chromosome = "23")
  edit(nxt(), base_pair_location = NA)
  edit(nxt(), base_pair_location = 0L)
  edit(nxt(), base_pair_location = -3L)
  edit(nxt(), base_pair_location = 2000000000L)
  edit(nxt(), reference_allele = NA)
  edit(nxt(), alternate_allele = ".")
  edit(nxt(), effect_allele = "")
  edit(nxt(), other_allele = "NA")
  i <- nxt(); edit(i, alternate_allele = "A,G", effect_allele = "A,G")
  i <- nxt(); edit(i, alternate_allele = "<INS>", effect_allele = "<INS>")
  i <- nxt(); edit(i, other_allele = "*", reference_allele = "*")
  i <- nxt(); edit(i, alternate_allele = "ACGT", effect_allele = "ACGT")
  i <- nxt(); edit(i, alternate_allele = "N", effect_allele = "N")
  i <- nxt(); edit(i, alternate_allele = " c ", effect_allele = "c")
  i <- nxt(); edit(i, alternate_allele = data$reference_allele[i],
                   effect_allele = data$reference_allele[i])
  i <- nxt(); edit(i, effect_allele = data$other_allele[i], other_allele = data$effect_allele[i])
  i <- nxt(); edit(i, effect_allele = data$other_allele[i])
  edit(nxt(), beta = Inf)
  edit(nxt(), beta = -Inf)
  edit(nxt(), beta = NA)
  edit(nxt(), beta = NaN)
  edit(nxt(), standard_error = Inf)
  edit(nxt(), standard_error = 0)
  edit(nxt(), standard_error = -1)
  edit(nxt(), standard_error = NA)
  edit(nxt(), effect_allele_frequency = 1.0000001)
  edit(nxt(), effect_allele_frequency = -0.5)
  edit(nxt(), effect_allele_frequency = Inf)
  edit(nxt(), effect_allele_frequency = NaN)
  edit(nxt(), effect_allele_frequency = NA)
  edit(nxt(), p_value = 1.5)
  edit(nxt(), p_value = -0.1)
  edit(nxt(), p_value = -Inf)
  edit(nxt(), p_value = 0)
  edit(nxt(), p_value = NaN)
  # Duplicates: exact copies, a copy with a second reason, an unsorted triple.
  extra <- data[c(200, 200, 250, 260, 270, 270, 280), ]
  extra$beta[3] <- Inf
  data <- rbind(data, extra)
  data <- data[c(seq_len(100), sample(101:nrow(data))), ]
  rownames(data) <- NULL
  data
}

qc_imported <- function(data, construct_variant_id = FALSE, project = TRUE) {
  CompreSSoR:::import_sumstats_impl(
    data, input_build = "GRCh38", row_policy = "report",
    project_columns = project, core_only = project, run_qc = FALSE,
    construct_variant_id = construct_variant_id, include_p_value = TRUE
  )
}

expect_qc_identical <- function(data, ...) {
  for (detail in c("compact", "full")) {
    new <- CompreSSoR:::structural_qc_report(data, detail = detail, ...)
    ref <- reference_structural_qc_report(data, detail = detail, ...)
    for (field in c("input_rows", "rejection_counts", "counts", "rejections",
                    "examples", "rejected_rows", "valid", "missing_columns")) {
      expect_identical(new[[field]], ref[[field]], info = paste(detail, field))
    }
    if (identical(detail, "compact")) {
      expect_identical(new$internal$invalid, ref$internal$invalid)
      expect_identical(new$internal$duplicate, ref$internal$duplicate)
      expect_identical(new$internal$canonical_key, ref$internal$canonical_key)
      duplicate <- which(ref$internal$duplicate)
      expect_identical(!new$internal$invalid_before_duplicate[duplicate],
                       ref$internal$reason_count[duplicate] == 1L)
      expect_identical(new$internal$invalid_rows, which(ref$internal$invalid))
      expect_identical(new$internal$duplicate_rows, which(ref$internal$duplicate))
    } else {
      for (field in c("row_status", "canonical_key", "invalid_rows",
                      "duplicate_rows", "structurally_valid_rows")) {
        expect_identical(new[[field]], ref[[field]], info = paste(detail, field))
      }
    }
  }
}

test_that("structural QC reports match the reference on crafted bad rows", {
  data <- qc_adversarial_frame()
  expect_qc_identical(qc_imported(data))
  expect_qc_identical(qc_imported(data), require_statistics = FALSE)
  expect_qc_identical(qc_imported(data), max_examples = 0L)
  expect_qc_identical(qc_imported(data), input_build = "GRCh37")
  # Fewer than 64 rows takes the non-hashed distinct path.
  expect_qc_identical(qc_imported(data[1:40, ]))
  # A clean, sorted table: every reason empty.
  clean <- qc_adversarial_frame(seed = 2L)[150:300, ]
  clean <- clean[order(as.integer(clean$chromosome), clean$base_pair_location), ]
  expect_qc_identical(qc_imported(clean))
})

test_that("structural QC matches the reference with aliases, parse failures and extra fields", {
  data <- qc_adversarial_frame(seed = 3L)
  n <- nrow(data)
  # Aliases excuse missing chromosome/coordinate rows.
  data$rsid <- ifelse(seq_len(n) %% 5L == 0L, paste0("rs", seq_len(n)), NA_character_)
  data$variant_id <- ifelse(seq_len(n) %% 7L == 0L, "1:100:A:G", NA_character_)
  data$chromosome[c(5, 10, 14)] <- NA
  data$base_pair_location[c(15, 21)] <- NA
  # Malformed text that the report-mode parser records as parse failures.
  data$beta <- as.character(data$beta)
  data$beta[c(8, 30)] <- "not-a-number"
  data$base_pair_location <- as.character(data$base_pair_location)
  data$base_pair_location[c(35, 70)] <- c("12x", "1.5")
  data$minus_log10_p <- -log10(runif(n))
  data$minus_log10_p[c(11, 12)] <- c(-1, Inf)
  data$sample_size <- 1000
  data$sample_size[c(13, 16)] <- c(0, Inf)
  data$info <- runif(n)
  data$info[c(17, 19)] <- c(1.2, -Inf)
  data$odds_ratio <- exp(suppressWarnings(as.numeric(data$beta)))
  data$odds_ratio[c(22, 23)] <- c(-1, Inf)
  data$z <- suppressWarnings(as.numeric(data$beta)) / data$standard_error
  data$z[24] <- Inf
  imported <- suppressWarnings(qc_imported(data, construct_variant_id = TRUE, project = FALSE))
  expect_true(length(attr(imported, "parse_failures")) > 0L)
  expect_qc_identical(imported)
})

test_that("structural QC matches the reference on random tables", {
  for (seed in 11:16) {
    data <- qc_adversarial_frame(n = 300L, seed = seed)
    shuffle <- sample(nrow(data))
    expect_qc_identical(qc_imported(data[shuffle, ]))
  }
})

test_that("identity reuse after QC equals full encoding", {
  data <- qc_adversarial_frame(seed = 4L)
  imported <- qc_imported(data)
  qc <- CompreSSoR:::apply_structural_qc(imported, row_policy = "report",
                                         detail = "compact")
  expect_false(is.null(qc$identity))
  reused <- CompreSSoR:::canonicalize_core_identity(qc$data, build = "GRCh38",
                                                    include_variant_id = FALSE,
                                                    identity = qc$identity)
  fresh <- CompreSSoR:::canonicalize_core_identity(qc$data, build = "GRCh38",
                                                   include_variant_id = FALSE)
  expect_identical(reused, fresh)
  # A kept row without a canonical key (alias-excused missing chromosome)
  # makes QC withhold the codes, so the full encoder runs (and rejects it).
  data$rsid <- paste0("rs", seq_len(nrow(data)))
  data$chromosome[2] <- NA
  imported <- qc_imported(data, construct_variant_id = TRUE, project = FALSE)
  qc <- CompreSSoR:::apply_structural_qc(imported, row_policy = "report",
                                         detail = "compact")
  expect_null(qc$identity)
})

test_that("p-value resolution and position gaps match the reference", {
  set.seed(5)
  n <- 5000L
  base <- data.frame(z = c(rnorm(n - 3L, 0, 3), NA, Inf, NaN))
  cases <- list(
    absent = base,
    supplied = transform(base, p_value = runif(n)),
    mixed = transform(base, p_value = ifelse(runif(n) < 0.2, NA,
                                            ifelse(runif(n) < 0.1, 1.5, runif(n))))
  )
  attr(cases$mixed, "p_value_source_present") <- TRUE
  hidden <- cases$supplied
  attr(hidden, "p_value_source_present") <- FALSE
  cases$hidden <- hidden
  for (name in names(cases)) {
    expect_identical(CompreSSoR:::pcodec_native_pvalue_resolve(cases[[name]]),
                     reference_pvalue_resolve(cases[[name]]), info = name)
  }
  positions <- cumsum(sample.int(1000L, 70000L, TRUE)) - 1
  for (block in c(1L, 1024L, 8192L, 65536L, 100000L)) {
    expect_identical(CompreSSoR:::pcodec_native_position_gaps(positions, block),
                     reference_position_gaps(positions, block))
  }
  expect_identical(CompreSSoR:::pcodec_native_position_gaps(numeric(), 8192L),
                   reference_position_gaps(numeric(), 8192L))
  expect_identical(CompreSSoR:::pcodec_native_position_gaps(7, 8192L),
                   reference_position_gaps(7, 8192L))
  for (bad in list(c(1, 0), c(1, NA), c(1, Inf), c(0.5, 2), c(-1, 2), c(1, 2^32))) {
    expect_error(CompreSSoR:::pcodec_native_position_gaps(bad, 8192L))
    expect_error(reference_position_gaps(bad, 8192L))
  }
})

test_that("parse and normalisation shortcuts return the reference values", {
  x <- c(5L, NA, -3L, 2147483647L)
  expect_identical(CompreSSoR:::parse_integer_column(x, "x"), x)
  d <- c(1.5, NaN, NA, -Inf)
  parsed <- CompreSSoR:::parse_numeric_column(d, "d")
  expect_identical(is.na(parsed), c(FALSE, TRUE, TRUE, FALSE))
  expect_false(any(is.nan(parsed)))
  clean <- c(0.1, 0.2)
  expect_identical(CompreSSoR:::parse_numeric_column(clean, "clean"), clean)
  alleles <- rep(c("A", "C", "G", "T"), 100)
  expect_identical(CompreSSoR:::normalise_allele_vector(alleles), alleles)
  messy <- rep(c(" a", "C", "g "), 100)
  expect_identical(CompreSSoR:::normalise_allele_vector(messy),
                   toupper(trimws(messy)))
  chromosomes <- rep(c("1", "chr2", "23", "X"), 50)
  expect_identical(CompreSSoR:::normalise_chromosome(chromosomes),
                   rep(c("1", "2", "X", "X"), 50))
  canonical <- rep(c("1", "22", "X"), 50)
  expect_identical(CompreSSoR:::normalise_chromosome(canonical), canonical)
  expect_identical(CompreSSoR:::match_unique(rep(c("C", "Z", NA), 30), c("A", "C")),
                   match(rep(c("C", "Z", NA), 30), c("A", "C")))
})

test_that("the column content hash is XXH64 over the documented byte stream", {
  skip_if_not(is.loaded("compressor_column_xxh64", PACKAGE = "CompreSSoR"))
  xxh <- function(columns, threads = 1L) {
    .Call("compressor_column_xxh64", columns, as.integer(threads), PACKAGE = "CompreSSoR")
  }
  for (n in c(0L, 1L, 7L, 8L, 9L, 33L, 1000L)) {
    ints <- seq_len(n) * 7919L - 3L
    expected <- digest::digest(writeBin(ints, raw(), size = 4L, endian = "little"),
                               algo = "xxhash64", serialize = FALSE)
    expect_identical(xxh(list(ints)), expected)
    reals <- ints / 3
    expected <- digest::digest(writeBin(reals, raw(), size = 8L, endian = "little"),
                               algo = "xxhash64", serialize = FALSE)
    expect_identical(xxh(list(reals)), expected)
  }
  strings <- c("A", NA, "", "chr10", strrep("x", 70000))
  bytes <- unlist(lapply(strings, function(s) {
    if (is.na(s)) return(as.raw(c(0xff, 0xff, 0xff, 0xff)))
    c(writeBin(nchar(s, type = "bytes"), raw(), size = 4L, endian = "little"), charToRaw(s))
  }))
  expect_identical(xxh(list(strings)),
                   digest::digest(bytes, algo = "xxhash64", serialize = FALSE))
  # NA and NaN are canonical, whatever their payload bits.
  other_nan <- readBin(as.raw(c(0, 0, 0, 0, 0, 0, 0xf8, 0xff)), "double", size = 8L,
                       endian = "little")
  expect_identical(xxh(list(c(1, NaN))), xxh(list(c(1, other_nan))))
  expect_false(identical(xxh(list(c(1, NaN))), xxh(list(c(1, NA)))))
  columns <- list(1:10, letters, c(TRUE, NA), runif(50))
  expect_identical(xxh(columns, 1L), xxh(columns, 4L))
  expect_identical(xxh(columns, 64L), xxh(columns, 2L))
})

test_that("the input-table provenance hash depends on content, not threads", {
  frame <- data.frame(chromosome = c("1", "2"), base_pair_location = c(10L, 20L),
                      beta = c(0.1, -0.2), stringsAsFactors = FALSE)
  h1 <- CompreSSoR:::columns_content_hash(frame, threads = 1L)
  expect_match(h1, "^[0-9a-f]{64}$")
  expect_identical(CompreSSoR:::columns_content_hash(frame, threads = 8L), h1)
  changed <- frame
  changed$beta[2] <- -0.2000001
  expect_false(identical(CompreSSoR:::columns_content_hash(changed), h1))
  renamed <- frame
  names(renamed)[3] <- "BETA"
  expect_false(identical(CompreSSoR:::columns_content_hash(renamed), h1))
  as_factor <- frame
  as_factor$chromosome <- factor(as_factor$chromosome)
  expect_match(CompreSSoR:::columns_content_hash(as_factor), "^[0-9a-f]{64}$")
  listed <- frame
  listed$extra <- I(list(1, 2))
  expect_match(CompreSSoR:::columns_content_hash(listed), "^[0-9a-f]{64}$")
})

test_that("manifest normalisation serializes like a JSON round trip", {
  rt <- function(x) {
    path <- tempfile(fileext = ".json")
    on.exit(unlink(path))
    CompreSSoR:::write_manifest(x, path)
    CompreSSoR:::read_manifest(path)
  }
  pretty <- function(x) {
    as.character(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null",
                                  digits = 17))
  }
  manifest <- list(
    format = "CompreSSoR", rows = 10L, ratio = 1 / 3, whole = 65536,
    big = 4294967295, tiny = 1e-300, flag = TRUE,
    vector = c(-3.5, 3.5), named_vector = c(a = 1L, b = 2L),
    strings = c("z", "standard_error"), empty = character(), empty_list = list(),
    named_empty = structure(list(), names = character()),
    with_na = c(1, NA), logical_na = NA, character_na = NA_character_,
    numeric_na = NA_real_, inf = c(1, Inf), factor = factor(c("x", "y")),
    frame = data.frame(a = 1:2, b = c("u", NA), stringsAsFactors = FALSE),
    nested = list(blocks = lapply(1:3, function(i) list(row_start = i - 1L, offset = i * 1.5)),
                  missing = NULL, inner = list(values = 1:4, label = "x")),
    matrix = matrix(1:4, 2)
  )
  normalised <- CompreSSoR:::manifest_json_normalise(manifest)
  expect_identical(pretty(normalised), pretty(rt(manifest)))
  expect_identical(pretty(normalised), pretty(rt(rt(manifest))))
  expect_identical(CompreSSoR:::pcodec_canonical_manifest_sha256(normalised),
                   CompreSSoR:::pcodec_canonical_manifest_sha256(rt(manifest)))
})

test_that("allele_columns writes the payload of the explicit REF/ALT table", {
  skip_if_not(CompreSSoR:::pcodec_native_available())
  set.seed(7)
  n <- 3000L
  bases <- c("A", "C", "G", "T")
  ref <- sample(bases, n, TRUE)
  alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1L), "")
  ssf <- data.frame(
    chromosome = "3", base_pair_location = sort(sample.int(1e7, n)),
    effect_allele = alt, other_allele = ref, beta = rnorm(n, 0, 0.05),
    standard_error = runif(n, 0.01, 0.1), effect_allele_frequency = runif(n),
    p_value = runif(n), stringsAsFactors = FALSE
  )
  prepared <- ssf
  prepared$reference_allele <- prepared$other_allele
  prepared$alternate_allele <- prepared$effect_allele
  payload <- function(store) {
    files <- setdiff(list.files(store$path), c("manifest.json", "manifest.sha256"))
    stats::setNames(unname(tools::md5sum(file.path(store$path, files))), files)
  }
  for (qc in c("compact", "none")) {
    a <- compress_sumstats(prepared, tempfile(fileext = ".cpr"), qc = qc)
    b <- compress_sumstats(ssf, tempfile(fileext = ".cpr"), qc = qc,
                           allele_columns = c(ref = "other_allele", alt = "effect_allele"))
    expect_identical(payload(b), payload(a))
    expect_identical(b$manifest$integrity$payload_sha256, a$manifest$integrity$payload_sha256)
    expect_identical(b$manifest$source$allele_columns,
                     list(reference_allele = "other_allele", alternate_allele = "effect_allele"))
    expect_null(a$manifest$source$allele_columns)
  }
  path <- tempfile(fileext = ".tsv")
  utils::write.table(ssf, path, sep = "\t", quote = FALSE, row.names = FALSE)
  from_file <- compress_sumstats(path, tempfile(fileext = ".cpr"),
                                 allele_columns = c(alt = "effect_allele", ref = "other_allele"))
  expect_identical(from_file$manifest$source$allele_columns$reference_allele, "other_allele")
  expect_error(compress_sumstats(ssf, tempfile(fileext = ".cpr")), "explicit REF and ALT")
  expect_error(compress_sumstats(ssf, tempfile(fileext = ".cpr"),
                                 allele_columns = c(ref = "a2", alt = "effect_allele")),
               "not in the input")
  expect_error(compress_sumstats(prepared, tempfile(fileext = ".cpr"),
                                 allele_columns = c(ref = "other_allele", alt = "effect_allele")),
               "explicit REF/ALT")
  for (bad in list("other_allele", c(ref = "x", alt = "x"), c(a = "x", b = "y"),
                   c(ref = NA, alt = "y"), c(ref = "x", alt = ""))) {
    expect_error(compress_sumstats(ssf, tempfile(fileext = ".cpr"), allele_columns = bad),
                 "allele_columns")
  }
})

test_that("data-frame provenance records the column hash; files keep sha256", {
  skip_if_not(CompreSSoR:::pcodec_native_available())
  example <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
  skip_if_not(nzchar(example))
  frame <- utils::read.delim(example, stringsAsFactors = FALSE)
  from_frame <- compress_sumstats(frame, tempfile(fileext = ".cpr"))
  source <- from_frame$manifest$source
  expect_identical(source$hash_algorithm, "xxh64_per_column_sha256_combined_v1")
  expect_match(source$hash, "^[0-9a-f]{64}$")
  expect_null(source$sha256)
  again <- compress_sumstats(frame, tempfile(fileext = ".cpr"), threads = 2L)
  expect_identical(again$manifest$source$hash, source$hash)
  from_file <- compress_sumstats(example, tempfile(fileext = ".cpr"))
  expect_identical(from_file$manifest$source$sha256,
                   digest::digest(example, algo = "sha256", file = TRUE))
  expect_null(from_file$manifest$source$hash_algorithm)
})

test_that("a data-frame store written by a27b32d still opens, validates and reads", {
  skip_if_not(CompreSSoR:::pcodec_native_available())
  # Written by a27b32d from a data frame, so it carries the earlier
  # source$sha256 (serialized-object SHA-256) and no source$hash_algorithm.
  fixture <- test_path("fixtures", "legacy-z9se6-a27b32d.cpr")
  skip_if_not(dir.exists(fixture))
  store_dir <- file.path(tempdir(), paste0("a27b32d-", Sys.getpid()))
  unlink(store_dir, recursive = TRUE)
  dir.create(store_dir)
  file.copy(list.files(fixture, full.names = TRUE), store_dir)
  store <- open_compressor(store_dir)
  expect_identical(store$manifest$source$kind, "data.frame")
  expect_match(store$manifest$source$sha256, "^[0-9a-f]{64}$")
  expect_null(store$manifest$source$hash_algorithm)
  expect_no_error(CompreSSoR:::pcodec_validate_store(store, full = TRUE))
  values <- read_sumstats(store)
  expect_identical(nrow(values), as.integer(store$manifest$n_rows))
})
