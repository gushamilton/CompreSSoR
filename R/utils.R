`%||%` <- function(x, y) if (is.null(x)) y else x

phase_clock <- function() unname(proc.time()[["elapsed"]])

phase_seconds <- function(start) {
  as.numeric(phase_clock() - start)
}

require_parquet_backend <- function(feature = "this operation", dplyr = FALSE) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(feature, " requires the optional R package 'arrow'", call. = FALSE)
  }
  if (isTRUE(dplyr) && !requireNamespace("dplyr", quietly = TRUE)) {
    stop(feature, " requires the optional R package 'dplyr'", call. = FALSE)
  }
  invisible(TRUE)
}

required_sumstats_columns <- function() {
  c("chromosome", "base_pair_location", "reference_allele",
    "alternate_allele", "effect_allele", "other_allele", "beta",
    "standard_error")
}

sumstats_column_alias_map <- function() {
  list(
    chromosome = c("chromosome", "chr", "CHR", "#chrom", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "bp", "BP", "GENPOS"),
    reference_allele = c("reference_allele", "reference", "ref", "REF"),
    alternate_allele = c("alternate_allele", "alternate", "alt", "ALT"),
    effect_allele = c("effect_allele", "ea", "EA", "a1", "A1",
                      "ALLELE1", "alleleB", "ALLELEB"),
    other_allele = c("other_allele", "oa", "NEA", "nea", "a2", "A2",
                     "ALLELE0", "alleleA", "ALLELEA"),
    # A single-letter B is used by some sources for expected non-reference
    # counts. It must not compete with an exact beta column (or be silently
    # interpreted as an effect size when its meaning is schema-dependent).
    beta = c("beta", "BETA", "effect", "effect_size", "estimate",
             "ES", "LOGOR", "LOG_OR", "log_odds", "BETA_LINREG", "BETA_LMM"),
    z = c("z", "Z", "zscore", "Z_SCORE", "z_stat", "ZSTAT", "zstatistic",
          "Z_STATISTIC", "Z_BOLT_LMM"),
    odds_ratio = c("odds_ratio", "OR", "ODDSRATIO", "oddsratio"),
    standard_error = c("standard_error", "SE", "se", "sebeta", "SEBETA", "stderr", "std_err"),
    effect_allele_frequency = c("effect_allele_frequency", "eaf", "EAF", "af", "AF",
                                "effect_af", "A1FREQ", "ALT_FREQ", "ALT_AF", "EAF_ALT"),
    p_value = c("p_value", "P", "p", "pval", "pvalue", "P_VALUE", "P_BOLT_LMM"),
    minus_log10_p = c("minus_log10_p", "LP", "lp", "neglog10p", "-LOG10P", "LOG10P"),
    variant_id = c("variant_id", "SNPID", "SNP_ID", "SNP", "snp", "variant", "ID"),
    sample_size = c("sample_size", "N", "n", "SS", "N_TOTAL", "TOTALSAMPLESIZE", "N_SUMMARY"),
    info = c("info", "INFO_SCORE", "IMPINFO", "INFO", "R2", "RSQ", "INFO_SCORE_IMP"),
    expected_count = c("expected_count", "B", "expected_nonref_count", "nonref_count")
  )
}

sumstats_resolution_matrix <- function() {
  data.frame(
    field = c("beta", "standard_error", "z", "odds_ratio", "p_value"),
    primary_route = c(
      "explicit beta alias",
      "explicit standard_error alias",
      "explicit z alias checked against beta/standard_error, otherwise derived",
      "positive OR alias resolved to log(OR) when beta is absent",
      "input aid only; never stored"
    ),
    strict_fallback = c(
      "positive OR -> log(OR) only at the boundary when beta is absent",
      "none",
      "beta / standard_error",
      "positive OR -> log(OR) when beta is absent; SE remains required",
      "none; p-to-SE inference is disabled"
    ),
    stringsAsFactors = FALSE
  )
}

sumstats_resolution_contract <- function(allow_p_to_se = FALSE,
                                         p_to_se_rows = integer()) {
  list(
    id = "beta_se_primary_v1",
    beta = "explicit_beta_alias_or_unambiguous_positive_or_boundary_conversion",
    standard_error = "explicit_standard_error_alias_or_explicit_opt_in_p_to_se",
    z = "explicit_z_checked_against_beta_over_standard_error_or_derived",
    odds_ratio = "positive_OR_to_log_OR_only_when_beta_is_absent",
    p_value = "input_aid_only;_never_stored",
    p_to_se = list(
      enabled = isTRUE(allow_p_to_se),
      method = if (isTRUE(allow_p_to_se))
        "explicit_opt_in_conflict_checked" else "disabled_in_strict_core",
      rows = as.integer(length(p_to_se_rows)),
      example_rows = as.integer(utils::head(p_to_se_rows, 5L))
    )
  )
}

sumstats_projection_aliases <- function(core_only = FALSE,
                                        allow_p_to_se = FALSE) {
  alias_map <- sumstats_column_alias_map()
  if (isTRUE(core_only)) {
    alias_map <- alias_map[c(
      "chromosome", "base_pair_location", "reference_allele",
      "alternate_allele", "effect_allele", "other_allele", "beta",
      "standard_error", "effect_allele_frequency", "z", "odds_ratio"
    )]
  }
  if (isTRUE(allow_p_to_se) && isTRUE(core_only)) {
    alias_map$p_value <- sumstats_column_alias_map()$p_value
  }
  unique(c(
    if (isTRUE(core_only)) character() else
      c("rsid", "rsids", "rsID", "RSID", "rs_id", "RS_ID", "dbsnp"),
    unlist(alias_map, use.names = FALSE)
  ))
}

sumstats_projection_indices <- function(source_columns, core_only = FALSE,
                                        allow_p_to_se = FALSE) {
  source_columns <- as.character(source_columns %||% character())
  indices <- which(alias_key(source_columns) %in%
                     alias_key(sumstats_projection_aliases(
                       core_only = core_only, allow_p_to_se = allow_p_to_se
                     )))
  # Retain one source field when no recognized names exist. This preserves the
  # row count so the normal schema validator can issue its precise missing-
  # column diagnostic without materialising the rest of a wide input.
  if (!length(indices) && length(source_columns)) indices <- 1L
  indices
}

sumstats_primary_chromosomes <- function() {
  c(as.character(seq_len(22L)), "X", "Y")
}

sumstats_chromosome_lengths <- function(input_build = "GRCh38") {
  build <- toupper(normalise_build_name(input_build))
  lengths <- switch(
    build,
    GRCH37 = c(
      `1` = 249250621, `2` = 243199373, `3` = 198022430,
      `4` = 191154276, `5` = 180915260, `6` = 171115067,
      `7` = 159138663, `8` = 146364022, `9` = 141213431,
      `10` = 135534747, `11` = 135006516, `12` = 133851895,
      `13` = 115169878, `14` = 107349540, `15` = 102531392,
      `16` = 90354753, `17` = 81195210, `18` = 78077248,
      `19` = 59128983, `20` = 63025520, `21` = 48129895,
      `22` = 51304566, X = 155270560, Y = 59373566
    ),
    GRCH38 = c(
      `1` = 248956422, `2` = 242193529, `3` = 198295559,
      `4` = 190214555, `5` = 181538259, `6` = 170805979,
      `7` = 159345973, `8` = 145138636, `9` = 138394717,
      `10` = 133797422, `11` = 135086622, `12` = 133275309,
      `13` = 114364328, `14` = 107043718, `15` = 101991189,
      `16` = 90338345, `17` = 83257441, `18` = 80373285,
      `19` = 58617616, `20` = 64444167, `21` = 46709983,
      `22` = 50818468, X = 156040895, Y = 57227415
    ),
    stop("input_build must be GRCh37/hg19 or GRCh38/hg38", call. = FALSE)
  )
  lengths
}

normalise_build_name <- function(input_build) {
  if (length(input_build) != 1L || !is.character(input_build) ||
      is.na(input_build) || !nzchar(trimws(input_build))) {
    stop("input_build must be GRCh37/hg19 or GRCh38/hg38", call. = FALSE)
  }
  key <- toupper(gsub("[. _-]", "", trimws(input_build)))
  if (key %in% c("GRCH37", "HG19", "37")) return("GRCh37")
  if (key %in% c("GRCH38", "HG38", "38")) return("GRCh38")
  stop("input_build must be GRCh37/hg19 or GRCh38/hg38", call. = FALSE)
}

clean_input_names <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  names(data) <- sub("^\\ufeff", "", trimws(names(data)))
  names(data)[is.na(names(data)) | !nzchar(names(data))] <- paste0(
    "unnamed_", which(is.na(names(data)) | !nzchar(names(data)))
  )
  data
}

strict_fread <- function(...) {
  fatal <- NULL
  result <- tryCatch(
    withCallingHandlers(
      data.table::fread(...),
      warning = function(w) {
        message <- conditionMessage(w)
        if (grepl("stopped early|expected [0-9]+ fields|discarded single-line footer|line [0-9]+.*fields",
                  message, ignore.case = TRUE)) {
          fatal <<- message
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {
      stop("malformed delimited input: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.null(fatal)) stop("malformed delimited input: ", fatal, call. = FALSE)
  result
}

input_compression <- function(path) {
  name <- tolower(basename(as.character(path)))
  if (grepl("[.]bgz$", name)) return("bgz")
  if (grepl("[.]gz$", name)) return("gz")
  if (grepl("[.]bz2$", name)) return("bz2")
  if (grepl("[.]xz$", name)) return("xz")
  if (grepl("[.]zst$", name)) return("zst")
  "none"
}

compressed_read_command <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  compression <- input_compression(path)
  if (identical(compression, "none")) return(shQuote(path))
  override <- Sys.getenv("COMPRESSOR_DECOMPRESSOR", unset = "")
  decompressor <- if (nzchar(override)) {
    override
  } else {
    switch(compression,
      gz = if (nzchar(Sys.which("pigz"))) "pigz" else "gzip",
      bgz = if (nzchar(Sys.which("pigz"))) "pigz" else "gzip",
      bz2 = "bzip2",
      xz = "xz",
      zst = "zstd",
      "gzip"
    )
  }
  flags <- if (identical(compression, "zst")) "-dc --quiet" else "-dc"
  paste(shQuote(decompressor), flags, shQuote(path))
}

detect_delimited_separator <- function(lines) {
  lines <- lines[!grepl("^\\s*##", lines)]
  if (!length(lines)) return("")
  header <- lines[[1L]]
  candidates <- c("\t", ",", ";", "|")
  counts <- vapply(candidates, function(pattern) {
    hit <- gregexpr(pattern, header, fixed = TRUE)[[1L]]
    sum(hit > 0L)
  }, numeric(1))
  if (!max(counts)) return("")
  candidates[which.max(counts)]
}

numeric_missing <- function(x) {
  text <- trimws(as.character(x))
  is.na(x) | is.na(text) | !nzchar(text) |
    tolower(text) %in% c(".", "na", "nan", "null")
}

parse_numeric_column <- function(x, name, invalid = c("error", "report")) {
  invalid <- match.arg(invalid)
  # Preserve numeric inputs bit for bit. A double -> character -> double
  # round-trip can alter least-significant bits. Factors deliberately take the
  # text path below so their labels, rather than level codes, are parsed.
  if (is.numeric(x) && !is.factor(x)) {
    out <- as.numeric(x)
    out[is.na(out)] <- NA_real_
    return(out)
  }
  text <- trimws(as.character(x))
  missing <- numeric_missing(x)
  out <- suppressWarnings(as.numeric(text))
  failed <- !missing & is.na(out)
  if (any(failed)) {
    rows <- which(failed)
    if (identical(invalid, "error")) {
      preview <- paste(utils::head(rows, 5L), collapse = ", ")
      stop("cannot parse ", name, " as numeric at row(s): ", preview,
           if (length(rows) > 5L) ", ..." else "", call. = FALSE)
    }
  }
  out[missing] <- NA_real_
  if (any(failed)) attr(out, "parse_failures") <- which(failed)
  out
}

eaf_coverage_metadata <- function(values) {
  if (is.null(values)) values <- numeric()
  values <- if (is.numeric(values) && !is.factor(values)) {
    as.numeric(values)
  } else {
    suppressWarnings(as.numeric(as.character(values)))
  }
  present <- is.finite(values) & values >= 0 & values <= 1
  missing <- is.na(values) | !is.finite(values)
  invalid <- is.finite(values) & (values < 0 | values > 1)
  rows <- length(values)
  list(
    rows = as.integer(rows),
    present = as.integer(sum(present)),
    missing = as.integer(sum(missing)),
    invalid = as.integer(sum(invalid)),
    coverage = if (rows) as.numeric(sum(present) / rows) else NA_real_,
    definition = "present means finite EAF in [0, 1]; missing is non-finite EAF"
  )
}

parse_integer_column <- function(x, name, invalid = c("error", "report")) {
  invalid <- match.arg(invalid)
  value <- parse_numeric_column(x, name, invalid = invalid)
  parse_failures <- attr(value, "parse_failures") %||% integer()
  bad <- !is.na(value) & (!is.finite(value) | value != trunc(value) |
                            value < -.Machine$integer.max | value > .Machine$integer.max)
  if (any(bad)) {
    rows <- which(bad)
    if (identical(invalid, "error")) {
      stop(name, " must contain whole 32-bit integers; invalid row(s): ",
           paste(utils::head(rows, 5L), collapse = ", "), call. = FALSE)
    }
    parse_failures <- sort(unique(c(parse_failures, rows)))
    value[bad] <- NA_real_
  }
  out <- as.integer(value)
  if (length(parse_failures)) attr(out, "parse_failures") <- parse_failures
  out
}

alias_values_equal <- function(x, y) {
  x_text <- trimws(as.character(x))
  y_text <- trimws(as.character(y))
  x_missing <- numeric_missing(x)
  y_missing <- numeric_missing(y)
  same <- (x_missing & y_missing) | (!x_missing & !y_missing & x_text == y_text)
  x_num <- suppressWarnings(as.numeric(x_text))
  y_num <- suppressWarnings(as.numeric(y_text))
  numeric_same <- !x_missing & !y_missing & is.finite(x_num) & is.finite(y_num) &
    abs(x_num - y_num) <= 1e-12 * pmax(1, abs(x_num), abs(y_num))
  all(same | numeric_same)
}

alias_key <- function(x) gsub("[^a-z0-9#]", "", tolower(as.character(x)), perl = TRUE)

normalise_chromosome <- function(x) {
  out <- toupper(trimws(sub("^chr", "", as.character(x), ignore.case = TRUE)))
  out[out == "23"] <- "X"
  out[out == "24"] <- "Y"
  out[out %in% c("M", "25", "26")] <- "MT"
  out[out %in% c("", ".", "NA")] <- NA_character_
  out
}

first_alias_index <- function(data, aliases) {
  keys <- alias_key(names(data))
  hits <- match(alias_key(aliases), keys)
  hits <- hits[!is.na(hits)]
  if (!length(hits)) NA_integer_ else hits[1L]
}

parse_vcf_info <- function(info) {
  fields <- strsplit(ifelse(is.na(info), "", as.character(info)), ";", fixed = TRUE)
  out <- lapply(fields, function(parts) {
    result <- list()
    for (part in parts) {
      bits <- strsplit(part, "=", fixed = TRUE)[[1L]]
      if (length(bits) >= 2L) result[[alias_key(bits[1L])]] <- bits[2L]
    }
    result
  })
  keys <- unique(unlist(lapply(out, names), use.names = FALSE))
  empty_columns <- replicate(length(keys), rep(NA_character_, length(info)), simplify = FALSE)
  names(empty_columns) <- keys
  result <- as.data.frame(empty_columns, stringsAsFactors = FALSE, check.names = FALSE)
  for (key in keys) {
    result[[key]] <- vapply(out, function(row) row[[key]] %||% NA_character_, character(1))
  }
  result
}

parse_vcf_format <- function(format, sample) {
  format <- as.character(format)
  sample <- as.character(sample)
  keys_by_row <- strsplit(ifelse(is.na(format), "", format), ":", fixed = TRUE)
  values_by_row <- strsplit(ifelse(is.na(sample), "", sample), ":", fixed = TRUE)
  keys <- unique(unlist(keys_by_row, use.names = FALSE))
  keys <- alias_key(keys[nzchar(keys)])
  result <- as.data.frame(stats::setNames(
    replicate(length(keys), rep(NA_character_, length(format)), simplify = FALSE),
    keys
  ), stringsAsFactors = FALSE, check.names = FALSE)
  for (i in seq_along(format)) {
    row_keys <- alias_key(keys_by_row[[i]])
    row_values <- values_by_row[[i]]
    n <- min(length(row_keys), length(row_values))
    if (!n) next
    for (j in seq_len(n)) {
      if (nzchar(row_keys[j]) && row_keys[j] %in% names(result)) {
        result[[row_keys[j]]][i] <- row_values[j]
      }
    }
  }
  result
}

input_probe_lines <- function(input, n = 32L) {
  compression <- input_compression(input)
  connection <- switch(compression,
    gz = gzfile(input, open = "rt"),
    bgz = gzfile(input, open = "rt"),
    bz2 = bzfile(input, open = "rt"),
    xz = xzfile(input, open = "rt"),
    file(input, open = "rt")
  )
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  readLines(connection, n = n, warn = FALSE)
}

read_delimited_header <- function(input, skip = NULL) {
  compression <- input_compression(input)
  probe <- if (!identical(compression, "zst")) {
    tryCatch(input_probe_lines(input), error = function(e) character())
  } else {
    character()
  }
  separator <- detect_delimited_separator(probe)
  arguments <- list(data.table = FALSE, showProgress = FALSE, check.names = FALSE,
                    fill = FALSE, nrows = 0L)
  if (nzchar(separator)) arguments$sep <- separator
  if (!is.null(skip)) arguments$skip <- skip
  if (identical(compression, "none")) {
    arguments$file <- input
  } else {
    arguments$cmd <- compressed_read_command(input)
  }
  clean_input_names(do.call(strict_fread, arguments))
}

read_delimited_file <- function(input, skip = NULL, select = NULL) {
  compression <- input_compression(input)
  probe <- if (!identical(compression, "zst")) {
    tryCatch(input_probe_lines(input), error = function(e) character())
  } else {
    character()
  }
  separator <- detect_delimited_separator(probe)
  arguments <- list(data.table = FALSE, showProgress = FALSE, check.names = FALSE,
                    fill = FALSE)
  if (nzchar(separator)) arguments$sep <- separator
  if (!is.null(skip)) arguments$skip <- skip
  if (!is.null(select)) arguments$select <- select
  if (identical(compression, "none")) {
    arguments$file <- input
  } else {
    # Use the explicit decompressor command for every compressed transport.
    # data.table's direct .gz path requires the optional R.utils package,
    # whereas the command path keeps gzip support available on a minimal R
    # installation and remains stream-oriented for large inputs.
    arguments$cmd <- compressed_read_command(input)
  }
  do.call(strict_fread, arguments)
}

read_vcf_input <- function(input, parse_policy = c("error", "report")) {
  parse_policy <- match.arg(parse_policy)
  data <- read_delimited_file(input, skip = "#CHROM")
  if (!ncol(data)) stop("VCF contains no header or records", call. = FALSE)
  data <- clean_input_names(data)
  chrom_i <- first_alias_index(data, c("#CHROM", "CHROM", "chromosome"))
  pos_i <- first_alias_index(data, c("POS", "position", "base_pair_location"))
  ref_i <- first_alias_index(data, c("REF", "other_allele"))
  alt_i <- first_alias_index(data, c("ALT", "effect_allele"))
  id_i <- first_alias_index(data, c("ID", "variant_id", "SNP"))
  info_i <- first_alias_index(data, c("INFO", "info"))
  format_i <- first_alias_index(data, c("FORMAT", "format"))
  required <- c(chrom_i, pos_i, ref_i, alt_i)
  if (anyNA(required)) stop("VCF requires #CHROM, POS, REF and ALT columns", call. = FALSE)
  alt <- as.character(data[[alt_i]])
  if (any(grepl(",", alt, fixed = TRUE), na.rm = TRUE)) {
    stop("VCF contains multiallelic ALT values; split multiallelic records before conversion", call. = FALSE)
  }
  info <- if (!is.na(info_i)) parse_vcf_info(data[[info_i]]) else data.frame()
  format_values <- data.frame()
  sample_i <- integer()
  if (!is.na(format_i)) {
    sample_i <- seq.int(format_i + 1L, ncol(data))
    sample_i <- sample_i[sample_i <= ncol(data)]
    if (length(sample_i) != 1L) {
      stop("GWAS-VCF FORMAT input must contain exactly one trait/sample column", call. = FALSE)
    }
    format_values <- parse_vcf_format(data[[format_i]], data[[sample_i]])
  }
  field <- function(aliases, direct = TRUE) {
    direct_i <- if (isTRUE(direct)) first_alias_index(data, aliases) else NA_integer_
    if (!is.na(direct_i)) return(as.character(data[[direct_i]]))
    # A GWAS-VCF's single trait/sample FORMAT values are trait-specific and
    # therefore take precedence over site-level INFO annotations.
    hits <- intersect(alias_key(aliases), names(format_values))
    if (length(hits)) return(format_values[[hits[1L]]])
    hits <- intersect(alias_key(aliases), names(info))
    if (length(hits)) return(info[[hits[1L]]])
    rep(NA_character_, nrow(data))
  }
  beta <- field(c("BETA", "ES", "EFFECT", "LOGOR"))
  or <- field(c("OR", "ODDSRATIO"))
  invalid_or <- is.na(suppressWarnings(as.numeric(beta))) &
    !is.na(suppressWarnings(as.numeric(or))) &
    (!is.finite(suppressWarnings(as.numeric(or))) | suppressWarnings(as.numeric(or)) <= 0)
  if (any(invalid_or) && identical(parse_policy, "error")) {
    stop("VCF odds ratio must be finite and positive", call. = FALSE)
  }
  beta_missing <- is.na(suppressWarnings(as.numeric(beta))) &
    is.finite(suppressWarnings(as.numeric(or))) & suppressWarnings(as.numeric(or)) > 0
  beta[beta_missing] <- as.character(log(as.numeric(or[beta_missing])))
  lp <- field(c("LP", "-LOG10P"))
  p_value <- field(c("P", "PVAL", "PVALUE"))
  lp_missing <- is.na(suppressWarnings(as.numeric(p_value))) & is.finite(suppressWarnings(as.numeric(lp)))
  p_value[lp_missing] <- as.character(10^(-as.numeric(lp[lp_missing])))
  out <- data.frame(
    chromosome = data[[chrom_i]],
    base_pair_location = data[[pos_i]],
    reference_allele = data[[ref_i]],
    alternate_allele = alt,
    effect_allele = alt,
    other_allele = data[[ref_i]],
    beta = beta,
    standard_error = field(c("SE", "STDERR", "SEBETA")),
    effect_allele_frequency = field(c("EAF", "AF", "A1FREQ")),
    p_value = p_value,
    variant_id = if (!is.na(id_i)) {
      ids <- as.character(data[[id_i]])
      ids[ids == "."] <- NA_character_
      ids
    } else NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  z <- field(c("Z", "ZSCORE", "Z_SCORE", "ZSTAT", "Z_STATISTIC"))
  out$z <- z
  out$sample_size <- field(c("SS", "N", "N_TOTAL", "TOTALSAMPLESIZE"))
  out$info <- field(c("INFO_SCORE", "IMPINFO", "R2", "RSQ", "INFO"), direct = FALSE)
  if (all(is.na(suppressWarnings(as.numeric(out$beta)))) &&
      any(is.finite(suppressWarnings(as.numeric(or))))) {
    out$beta <- as.character(log(as.numeric(or)))
  }
  if (all(is.na(suppressWarnings(as.numeric(out$standard_error)))) &&
      any(is.finite(suppressWarnings(as.numeric(out$z))))) {
    # SE still has to be supplied or derivable from p; this branch only keeps
    # the z column available for the common z/se input shape.
    out$standard_error <- field(c("SE", "STDERR", "SEBETA"))
  }
  exclude <- names(data)[required]
  if (!is.na(id_i)) exclude <- c(exclude, names(data)[id_i])
  if (!is.na(info_i)) exclude <- c(exclude, names(data)[info_i])
  if (!is.na(format_i)) exclude <- c(exclude, names(data)[format_i], names(data)[sample_i])
  extra_names <- setdiff(names(data), exclude)
  if (length(extra_names)) out[extra_names] <- data[extra_names]
  attr(out, "source_columns") <- names(data)
  out
}

looks_like_vcf <- function(input) {
  if (grepl("[.]vcf(?:[.]bgz|[.]gz|[.]bz2|[.]xz|[.]zst)?$", input,
            ignore.case = TRUE)) return(TRUE)
  probe <- tryCatch(input_probe_lines(input, n = 16L), error = function(e) character())
  any(grepl("^##fileformat[=]VCF", probe, ignore.case = TRUE)) ||
    any(grepl("^#CHROM(?:\\t|\\s)", probe, ignore.case = TRUE))
}

read_sumstats_input <- function(input, parse_policy = c("error", "report"),
                                project_columns = FALSE, core_only = FALSE,
                                allow_p_to_se = FALSE) {
  parse_policy <- match.arg(parse_policy)
  if (length(project_columns) != 1L || !is.logical(project_columns) || is.na(project_columns) ||
      length(core_only) != 1L || !is.logical(core_only) || is.na(core_only)) {
    stop("project_columns and core_only must be TRUE or FALSE", call. = FALSE)
  }
  if (length(allow_p_to_se) != 1L || !is.logical(allow_p_to_se) || is.na(allow_p_to_se)) {
    stop("allow_p_to_se must be TRUE or FALSE", call. = FALSE)
  }
  started <- unname(proc.time()[["elapsed"]])
  if (is.data.frame(input)) {
    data <- clean_input_names(input)
    source_columns <- names(data)
    selected <- if (isTRUE(project_columns)) {
      sumstats_projection_indices(source_columns, core_only = core_only,
                                  allow_p_to_se = allow_p_to_se)
    } else {
      seq_along(source_columns)
    }
    data <- data[, selected, drop = FALSE]
    attr(data, "source_columns") <- source_columns
    attr(data, "source_columns_read") <- names(data)
    attr(data, "input_read_metadata") <- list(
      projected = isTRUE(project_columns),
      columns_before = as.integer(length(source_columns)),
      columns_read = as.integer(length(selected)),
      elapsed_seconds = unname(proc.time()[["elapsed"]]) - started
    )
    return(data)
  }
  if (length(input) != 1L || !is.character(input) || !file.exists(input)) {
    stop("input must be a data.frame or an existing delimited file path", call. = FALSE)
  }
  if (looks_like_vcf(input)) {
    data <- read_vcf_input(input, parse_policy = parse_policy)
    source_columns <- attr(data, "source_columns") %||% names(data)
    attr(data, "source_columns_read") <- source_columns
    attr(data, "input_read_metadata") <- list(
      projected = FALSE,
      columns_before = as.integer(length(source_columns)),
      columns_read = as.integer(length(source_columns)),
      elapsed_seconds = unname(proc.time()[["elapsed"]]) - started
    )
    return(data)
  }
  header <- read_delimited_header(input)
  source_columns <- names(header)
  selected <- if (isTRUE(project_columns)) {
    sumstats_projection_indices(source_columns, core_only = core_only,
                                allow_p_to_se = allow_p_to_se)
  } else {
    seq_along(source_columns)
  }
  data <- clean_input_names(read_delimited_file(input, select = selected))
  attr(data, "source_columns") <- source_columns
  attr(data, "source_columns_read") <- names(data)
  attr(data, "input_read_metadata") <- list(
    projected = isTRUE(project_columns),
    columns_before = as.integer(length(source_columns)),
    columns_read = as.integer(length(selected)),
    elapsed_seconds = unname(proc.time()[["elapsed"]]) - started
  )
  data
}

rename_first_alias <- function(data, target, aliases) {
  keys <- alias_key(names(data))
  aliases <- unique(c(target, aliases))
  alias_keys <- alias_key(aliases)
  hits <- which(keys %in% alias_keys)
  if (!length(hits)) return(data)
  priority <- match(keys[hits], alias_keys)
  hit <- hits[order(priority, hits)][1L]
  redundant <- setdiff(hits, hit)
  if (length(redundant)) {
    conflicts <- redundant[!vapply(redundant, function(i) {
      alias_values_equal(data[[hit]], data[[i]])
    }, logical(1))]
    if (length(conflicts)) {
      stop("conflicting aliases for ", target, ": ",
           paste(names(data)[c(hit, conflicts)], collapse = ", "), call. = FALSE)
    }
    names(data)[hit] <- target
    data <- data[-redundant]
    return(data)
  }
  names(data)[hit] <- target
  data
}

normalise_sumstats_columns <- function(data, parse_policy = c("error", "report"),
                                       allow_p_to_se = FALSE,
                                       construct_variant_id = TRUE) {
  parse_policy <- match.arg(parse_policy)
  if (length(allow_p_to_se) != 1L || !is.logical(allow_p_to_se) ||
      is.na(allow_p_to_se)) {
    stop("allow_p_to_se must be TRUE or FALSE", call. = FALSE)
  }
  data <- clean_input_names(data)
  parse_failures <- list()
  assign_parsed <- function(target, value) {
    failures <- attr(value, "parse_failures") %||% integer()
    if (length(failures)) parse_failures[[target]] <<- as.integer(failures)
    data[[target]] <<- value
  }
  data <- rename_first_alias(data, "rsid",
                             c("rsid", "rsids", "rsID", "RSID", "rs_id", "RS_ID", "dbsnp"))
  alias_map <- sumstats_column_alias_map()
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  beta_source_present <- "beta" %in% names(data)
  standard_error_source_present <- "standard_error" %in% names(data)
  z_source_present <- "z" %in% names(data)
  odds_ratio_source_present <- "odds_ratio" %in% names(data)
  # Keep source-level missingness for provenance, but delay the parse-policy
  # error until the unambiguous positive-OR boundary conversion has had a
  # chance to provide beta. The strict writer still checks source columns and
  # therefore never treats OR as a strict-core beta field.
  missing <- setdiff(required_sumstats_columns(), names(data))
  for (column in missing) data[[column]] <- NA_character_
  if (!"chromosome" %in% names(data)) data$chromosome <- NA_character_
  if (!"base_pair_location" %in% names(data)) data$base_pair_location <- NA_integer_
  data$chromosome <- normalise_chromosome(data$chromosome)
  assign_parsed("base_pair_location",
                parse_integer_column(data$base_pair_location, "base_pair_location",
                                     invalid = parse_policy))
  data$reference_allele <- toupper(trimws(as.character(data$reference_allele)))
  data$alternate_allele <- toupper(trimws(as.character(data$alternate_allele)))
  data$effect_allele <- toupper(trimws(as.character(data$effect_allele)))
  data$other_allele <- toupper(trimws(as.character(data$other_allele)))
  data$reference_allele[data$reference_allele %in% c("", ".", "NA", "N/A")] <- NA_character_
  data$alternate_allele[data$alternate_allele %in% c("", ".", "NA", "N/A")] <- NA_character_
  data$effect_allele[data$effect_allele %in% c("", ".", "NA", "N/A")] <- NA_character_
  data$other_allele[data$other_allele %in% c("", ".", "NA", "N/A")] <- NA_character_
  if (!"beta" %in% names(data)) data$beta <- NA_real_
  if (!"standard_error" %in% names(data)) data$standard_error <- NA_real_
  if (!"z" %in% names(data)) data$z <- NA_real_
  if (!"odds_ratio" %in% names(data)) data$odds_ratio <- NA_real_
  assign_parsed("beta", parse_numeric_column(data$beta, "beta", invalid = parse_policy))
  assign_parsed("standard_error", parse_numeric_column(data$standard_error, "standard_error",
                                                         invalid = parse_policy))
  assign_parsed("z", parse_numeric_column(data$z, "z", invalid = parse_policy))
  assign_parsed("odds_ratio", parse_numeric_column(data$odds_ratio, "odds_ratio",
                                                     invalid = parse_policy))
  bad_or <- !is.na(data$odds_ratio) & (!is.finite(data$odds_ratio) | data$odds_ratio <= 0)
  if (any(bad_or) && identical(parse_policy, "error")) {
    stop("odds_ratio must be finite and positive when supplied", call. = FALSE)
  }
  if (beta_source_present && odds_ratio_source_present) {
    comparable <- is.finite(data$beta) & is.finite(data$odds_ratio) &
      data$odds_ratio > 0
    if (any(comparable)) {
      beta_from_or <- log(data$odds_ratio[comparable])
      tolerance <- 1e-6 + 0.01 * pmax(abs(data$beta[comparable]),
                                       abs(beta_from_or))
      conflict <- abs(data$beta[comparable] - beta_from_or) > tolerance
      if (any(conflict)) {
        rows <- which(comparable)[conflict]
        stop("beta and odds_ratio are inconsistent at row(s): ",
             paste(utils::head(rows, 5L), collapse = ", "), call. = FALSE)
      }
    }
  }
  beta_from_or <- integer()
  if (odds_ratio_source_present) {
    beta_from_or <- which(is.na(data$beta) & is.finite(data$odds_ratio) &
                            data$odds_ratio > 0)
    if (length(beta_from_or)) data$beta[beta_from_or] <- log(data$odds_ratio[beta_from_or])
  }
  if (!"effect_allele_frequency" %in% names(data)) data$effect_allele_frequency <- NA_real_
  assign_parsed("effect_allele_frequency",
                parse_numeric_column(data$effect_allele_frequency,
                                     "effect_allele_frequency", invalid = parse_policy))
  if (!"p_value" %in% names(data)) data$p_value <- NA_real_
  assign_parsed("p_value", parse_numeric_column(data$p_value, "p_value", invalid = parse_policy))
  if ("minus_log10_p" %in% names(data)) {
    lp <- parse_numeric_column(data$minus_log10_p, "minus_log10_p", invalid = parse_policy)
    lp_failures <- attr(lp, "parse_failures") %||% integer()
    if (length(lp_failures)) parse_failures$minus_log10_p <- as.integer(lp_failures)
    fill <- is.na(data$p_value) & is.finite(lp)
    data$p_value[fill] <- 10^(-lp[fill])
  }
  # The strict route is explicit beta + standard error. An optional supplied Z
  # is checked against beta / standard_error, or derived when absent. P-value
  # to SE inference is deliberately opt-in and is conflict-checked below.
  z_missing <- is.na(data$z) & is.finite(data$beta) & is.finite(data$standard_error) & data$standard_error > 0
  data$z[z_missing] <- data$beta[z_missing] / data$standard_error[z_missing]
  p_to_se_rows <- integer()
  if (isTRUE(allow_p_to_se)) {
    p_to_se_candidate <- is.finite(data$beta) & is.finite(data$p_value) &
      data$p_value > 0 & data$p_value < 1 & data$beta != 0
    derived_se <- rep(NA_real_, nrow(data))
    derived_se[p_to_se_candidate] <- abs(data$beta[p_to_se_candidate]) /
      abs(stats::qnorm(data$p_value[p_to_se_candidate] / 2))
    supplied_se <- is.finite(data$standard_error) & data$standard_error > 0
    conflict <- p_to_se_candidate & supplied_se &
      abs(data$standard_error - derived_se) >
        (1e-6 + 0.01 * pmax(abs(data$standard_error), abs(derived_se)))
    if (any(conflict)) {
      stop("standard_error conflicts with the explicitly enabled p-value-derived value at row(s): ",
           paste(utils::head(which(conflict), 5L), collapse = ", "), call. = FALSE)
    }
    p_to_se_rows <- which(p_to_se_candidate & !supplied_se)
    if (length(p_to_se_rows)) data$standard_error[p_to_se_rows] <- derived_se[p_to_se_rows]
  }
  if (identical(parse_policy, "error")) {
    resolved_missing <- missing
    if (length(beta_from_or)) resolved_missing <- setdiff(resolved_missing, "beta")
    if (isTRUE(allow_p_to_se) && length(p_to_se_rows)) {
      resolved_missing <- setdiff(resolved_missing, "standard_error")
    }
    if (length(resolved_missing)) {
      stop("Missing required summary-statistics columns: ",
           paste(resolved_missing, collapse = ", "), call. = FALSE)
    }
  }
  z_missing <- is.na(data$z) & is.finite(data$beta) & is.finite(data$standard_error) & data$standard_error > 0
  data$z[z_missing] <- data$beta[z_missing] / data$standard_error[z_missing]
  consistent <- z_source_present & is.finite(data$beta) & is.finite(data$z) &
    is.finite(data$standard_error) & data$standard_error > 0
  if (any(consistent)) {
    beta_from_z <- data$z[consistent] * data$standard_error[consistent]
    supplied_beta <- data$beta[consistent]
    tolerance <- 1e-6 + 0.01 * pmax(abs(supplied_beta), abs(beta_from_z))
    conflict <- abs(supplied_beta - beta_from_z) > tolerance
    if (any(conflict)) {
      rows <- which(consistent)[conflict]
      stop("beta, z and standard_error are inconsistent at row(s): ",
           paste(utils::head(rows, 5L), collapse = ", "), call. = FALSE)
    }
  }
  if (isTRUE(construct_variant_id)) {
    if (!"variant_id" %in% names(data) && "rsid" %in% names(data)) {
      data$variant_id <- as.character(data$rsid)
    }
    if (!"variant_id" %in% names(data)) {
      data$variant_id <- paste(data$chromosome, data$base_pair_location,
                               data$other_allele, data$effect_allele, sep = "_")
    } else {
      data$variant_id <- as.character(data$variant_id)
      missing_id <- is.na(data$variant_id) | !nzchar(data$variant_id) | data$variant_id == "."
      if (any(missing_id)) {
        data$variant_id[missing_id] <- paste(data$chromosome[missing_id],
                                             data$base_pair_location[missing_id],
                                             data$other_allele[missing_id],
                                             data$effect_allele[missing_id], sep = "_")
      }
    }
  }
  if (isTRUE(construct_variant_id) && !"rsid" %in% names(data)) {
    data$rsid <- NA_character_
  }
  if ("sample_size" %in% names(data)) {
    assign_parsed("sample_size", parse_numeric_column(data$sample_size, "sample_size",
                                                       invalid = parse_policy))
  }
  if ("info" %in% names(data)) {
    assign_parsed("info", parse_numeric_column(data$info, "info", invalid = parse_policy))
  }
  if (isTRUE(construct_variant_id)) {
    alias_id <- !is.na(data$variant_id) & grepl("^rs", data$variant_id, ignore.case = TRUE)
    fill_rsid <- (is.na(data$rsid) | !nzchar(data$rsid)) & alias_id
    data$rsid[fill_rsid] <- data$variant_id[fill_rsid]
  }
  attr(data, "parse_failures") <- parse_failures
  attr(data, "missing_columns") <- missing
  attr(data, "resolution_provenance") <- sumstats_resolution_contract(
    allow_p_to_se = allow_p_to_se, p_to_se_rows = p_to_se_rows
  )
  attr(data, "resolution_provenance")$beta_route <- if (length(beta_from_or)) {
    "positive_OR_to_log_OR_boundary_conversion"
  } else if (beta_source_present) {
    "explicit_beta_alias"
  } else {
    "missing"
  }
  attr(data, "resolution_provenance")$standard_error_route <- if (standard_error_source_present) {
    "explicit_standard_error_alias"
  } else if (length(p_to_se_rows)) {
    "explicit_opt_in_p_value_to_standard_error_conversion"
  } else {
    "missing"
  }
  attr(data, "resolution_provenance")$z_source <- if (z_source_present) {
    "supplied_and_checked_when_beta_and_standard_error_are_finite"
  } else {
    "derived_from_beta_over_standard_error"
  }
  data
}

# Fast boundary for inputs that have already been prepared into the exact
# compressor schema.  This deliberately does not resolve aliases, construct a
# temporary variant_id, derive OR, or attach row-level QC state.  Identity and
# codec code still perform their own strict safety checks later in the write.
normalise_prepared_core_columns <- function(data, input_build = "GRCh38") {
  required <- c("chromosome", "base_pair_location", "reference_allele",
                "alternate_allele", "effect_allele", "other_allele", "beta",
                "standard_error")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("qc='none' requires already-canonical columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- data
  out$chromosome <- normalise_chromosome(out$chromosome)
  out$base_pair_location <- suppressWarnings(as.numeric(as.character(out$base_pair_location)))
  for (field in c("reference_allele", "alternate_allele", "effect_allele", "other_allele")) {
    out[[field]] <- toupper(trimws(as.character(out[[field]])))
  }
  numeric_fast <- function(value) {
    if (is.numeric(value) && !is.factor(value)) return(as.numeric(value))
    suppressWarnings(as.numeric(as.character(value)))
  }
  for (field in c("beta", "standard_error", "effect_allele_frequency", "z")) {
    if (field %in% names(out)) out[[field]] <- numeric_fast(out[[field]])
  }
  if (!"effect_allele_frequency" %in% names(out)) {
    out$effect_allele_frequency <- rep(NA_real_, nrow(out))
  }
  z_supplied <- "z" %in% names(out)
  if (!"z" %in% names(out)) out$z <- rep(NA_real_, nrow(out))
  z_supplied <- z_supplied & !is.na(out$z)
  derive_z <- is.na(out$z) & is.finite(out$beta) & is.finite(out$standard_error) &
    out$standard_error > 0
  out$z[derive_z] <- out$beta[derive_z] / out$standard_error[derive_z]
  # P-values are derived from Z at read time and are not part of the native
  # payload. Preserve a supplied canonical p_value for callers that explicitly
  # provided it, but do not allocate a full-length placeholder on the fast
  # prepared-input path.
  if ("p_value" %in% names(out)) out$p_value <- numeric_fast(out$p_value)
  attr(out, "missing_columns") <- character()
  attr(out, "parse_failures") <- list()
  attr(out, "z_supplied") <- z_supplied
  attr(out, "resolution_provenance") <- sumstats_resolution_contract(allow_p_to_se = FALSE)
  attr(out, "resolution_provenance")$beta_route <- "explicit_beta_canonical"
  attr(out, "resolution_provenance")$standard_error_route <- "explicit_standard_error_canonical"
  attr(out, "resolution_provenance")$z_source <- if (any(derive_z)) {
    "derived_from_beta_over_standard_error"
  } else {
    "supplied_canonical"
  }
  attr(out, "input_build") <- compressor_normalize_build(input_build)
  out
}

new_structural_qc_report <- function(n, input_build, max_examples = 5L,
                                     detail = c("full", "compact")) {
  detail <- match.arg(detail)
  report <- list(
    input_rows = as.integer(n),
    accepted_rows = NA_integer_,
    rejected_rows = NA_integer_,
    dropped_rows = NA_integer_,
    input_build = normalise_build_name(input_build),
    valid = FALSE,
    rejection_counts = integer(),
    examples = list(),
    max_examples = as.integer(max_examples),
    detail = detail
  )
  if (identical(detail, "full")) {
    report$row_status <- data.frame(
      row = seq_len(n), structurally_valid = logical(n), accepted = logical(n),
      reasons = rep("", n), stringsAsFactors = FALSE
    )
  }
  report
}

structural_qc_report <- function(data, input_build = "GRCh38",
                                 require_statistics = TRUE, max_examples = 5L,
                                 detail = c("full", "compact")) {
  detail <- match.arg(detail)
  data <- clean_input_names(data)
  n <- nrow(data)
  report <- new_structural_qc_report(
    n, input_build, max_examples = max_examples, detail = detail
  )
  report$missing_columns <- attr(data, "missing_columns") %||% character()
  if (!n) {
    report$accepted_rows <- 0L
    report$rejected_rows <- 0L
    report$dropped_rows <- 0L
    report$valid <- TRUE
    if (identical(detail, "full")) {
      report$row_status$structurally_valid <- logical()
      report$row_status$accepted <- logical()
    }
    return(report)
  }

  reasons <- if (identical(detail, "full")) list() else NULL
  rejection_counts <- integer()
  examples <- list()
  invalid <- rep(FALSE, n)
  reason_count <- integer(n)
  add_reason <- function(name, mask) {
    mask <- as.logical(mask)
    mask[is.na(mask)] <- FALSE
    if (length(mask) != n) stop("structural QC mask has the wrong number of rows",
                                call. = FALSE)
    if (identical(detail, "full")) reasons[[name]] <<- mask
    count <- sum(mask)
    rejection_counts[name] <<- as.integer(count)
    examples[[name]] <<- utils::head(which(mask), max(0L, as.integer(max_examples)))
    invalid <<- invalid | mask
    reason_count <<- reason_count + as.integer(mask)
  }
  values <- function(name, default = NA_real_) {
    if (name %in% names(data)) return(data[[name]])
    rep(default, n)
  }
  text_values <- function(name) {
    value <- values(name, NA_character_)
    out <- trimws(as.character(value))
    out[is.na(value) | out %in% c("", ".", "NA", "N/A")] <- NA_character_
    out
  }

  chromosome <- normalise_chromosome(values("chromosome", NA_character_))
  position <- suppressWarnings(as.numeric(values("base_pair_location", NA_real_)))
  reference <- toupper(text_values("reference_allele"))
  alternate <- toupper(text_values("alternate_allele"))
  effect <- toupper(text_values("effect_allele"))
  other <- toupper(text_values("other_allele"))
  rsid <- text_values("rsid")
  variant_id <- text_values("variant_id")
  alias_available <- (!is.na(rsid) & nzchar(rsid)) |
    grepl("^rs", variant_id, ignore.case = TRUE) |
    grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$",
          variant_id, ignore.case = TRUE)
  lengths <- sumstats_chromosome_lengths(input_build)
  primary <- sumstats_primary_chromosomes()

  add_reason("missing_chromosome",
             (is.na(chromosome) | !nzchar(chromosome)) & !alias_available)
  add_reason("unsupported_contig", !is.na(chromosome) & nzchar(chromosome) &
               !(chromosome %in% primary))
  add_reason("missing_coordinate", is.na(position) & !alias_available)
  parse_failures <- attr(data, "parse_failures") %||% list()
  add_reason("malformed_coordinate",
             position %in% c(NA_real_, NaN) &
               seq_len(n) %in% (parse_failures$base_pair_location %||% integer()) &
               !alias_available)
  add_reason("nonpositive_coordinate", !is.na(position) & is.finite(position) & position < 1)
  add_reason("noninteger_coordinate", !is.na(position) & is.finite(position) &
               position != trunc(position))
  known_chromosome <- !is.na(chromosome) & chromosome %in% names(lengths)
  add_reason("coordinate_out_of_range", known_chromosome & is.finite(position) &
               position >= 1 & position > unname(lengths[chromosome]))

  add_reason("missing_reference_allele", is.na(reference) | !nzchar(reference))
  add_reason("missing_alternate_allele", is.na(alternate) | !nzchar(alternate))
  add_reason("missing_effect_allele", is.na(effect) | !nzchar(effect))
  add_reason("missing_other_allele", is.na(other) | !nzchar(other))
  allele_values <- list(reference, alternate, effect, other)
  any_allele <- function(predicate) {
    Reduce(`|`, Map(function(value) {
      present <- !is.na(value)
      present & predicate(value)
    }, allele_values), init = rep(FALSE, n))
  }
  add_reason("multiallelic", any_allele(function(value) grepl(",", value, fixed = TRUE)))
  add_reason("symbolic_allele", any_allele(function(value) grepl("^<|^\\*", value)))
  add_reason("indel", any_allele(function(value) nchar(value) != 1L))
  add_reason("invalid_allele", any_allele(function(value) !grepl("^[ACGT]$", value)))
  add_reason("same_alleles", (
    !is.na(reference) & !is.na(alternate) & reference == alternate
  ) | (
    !is.na(effect) & !is.na(other) & effect == other
  ))
  orientation_defined <- !is.na(reference) & !is.na(alternate) &
    !is.na(effect) & !is.na(other)
  add_reason("orientation_mismatch", orientation_defined &
               (effect != alternate | other != reference))

  numeric_fields <- c("beta", "standard_error", "z", "effect_allele_frequency",
                      "p_value", "odds_ratio", "sample_size", "info")
  for (field in numeric_fields) {
    value <- suppressWarnings(as.numeric(values(field, NA_real_)))
    malformed <- seq_len(n) %in% (parse_failures[[field]] %||% integer())
    add_reason(paste0("malformed_", field), malformed)
    add_reason(paste0("nonfinite_", field), !is.na(value) & !is.finite(value))
  }
  beta <- suppressWarnings(as.numeric(values("beta", NA_real_)))
  se <- suppressWarnings(as.numeric(values("standard_error", NA_real_)))
  z <- suppressWarnings(as.numeric(values("z", NA_real_)))
  eaf <- suppressWarnings(as.numeric(values("effect_allele_frequency", NA_real_)))
  p_value <- suppressWarnings(as.numeric(values("p_value", NA_real_)))
  odds_ratio <- suppressWarnings(as.numeric(values("odds_ratio", NA_real_)))
  sample_size <- suppressWarnings(as.numeric(values("sample_size", NA_real_)))
  info <- suppressWarnings(as.numeric(values("info", NA_real_)))
  add_reason("invalid_standard_error", !is.na(se) & is.finite(se) & se <= 0)
  add_reason("invalid_effect_allele_frequency", !is.na(eaf) & is.finite(eaf) &
               (eaf < 0 | eaf > 1))
  add_reason("invalid_p_value", !is.na(p_value) & is.finite(p_value) &
               (p_value < 0 | p_value > 1))
  add_reason("invalid_odds_ratio", !is.na(odds_ratio) & is.finite(odds_ratio) & odds_ratio <= 0)
  add_reason("invalid_sample_size", !is.na(sample_size) & is.finite(sample_size) & sample_size <= 0)
  add_reason("invalid_info", !is.na(info) & is.finite(info) & (info < 0 | info > 1))
  if ("minus_log10_p" %in% names(data)) {
    lp <- suppressWarnings(as.numeric(data$minus_log10_p))
    add_reason("invalid_minus_log10_p", !is.na(lp) & is.finite(lp) & lp < 0)
    add_reason("nonfinite_minus_log10_p", !is.na(lp) & !is.finite(lp))
  }
  if (isTRUE(require_statistics)) {
    add_reason("missing_statistics", !is.finite(beta) | !is.finite(z) |
                 !is.finite(se) | se <= 0)
  }

  valid_key <- known_chromosome & is.finite(position) & position >= 1 &
    position == floor(position) & position <= unname(lengths[chromosome]) &
    !is.na(other) & !is.na(effect) & other %in% c("A", "C", "G", "T") &
    effect %in% c("A", "C", "G", "T") & other != effect
  key <- rep(NA_real_, n)
  if (any(valid_key)) {
    identity <- compressor_encode_variant_identity(
      chromosome[valid_key], position[valid_key], other[valid_key], effect[valid_key],
      build = input_build
    )
    key[valid_key] <- compressor_identity_code(identity$global_position,
                                               identity$substitution)
  }
  duplicate <- !is.na(key) & (duplicated(key) | duplicated(key, fromLast = TRUE))
  add_reason("duplicate_variant", duplicate)
  report$rejection_counts <- sort(rejection_counts, decreasing = TRUE)
  report$counts <- report$rejection_counts
  report$rejections <- report$rejection_counts[report$rejection_counts > 0L]
  report$examples <- examples
  report$rejected_rows <- as.integer(sum(invalid))
  report$valid <- !any(invalid)
  if (identical(detail, "full")) {
    reason_names <- names(reasons)
    reason_text <- rep("", n)
    for (name in reason_names) {
      rows <- reasons[[name]]
      reason_text[rows] <- ifelse(
        nzchar(reason_text[rows]), paste0(reason_text[rows], ",", name), name
      )
    }
    report$row_status$structurally_valid <- !invalid
    report$row_status$reasons <- reason_text
    report$canonical_key <- key
    report$invalid_rows <- which(invalid)
    report$duplicate_rows <- which(duplicate)
    report$structurally_valid_rows <- which(!invalid)
  } else {
    # Private, short-lived state consumed by apply_structural_qc(). It is
    # removed before the compact report is attached to a store or returned.
    report$internal <- list(
      invalid = invalid,
      duplicate = duplicate,
      canonical_key = key,
      reason_count = reason_count
    )
  }
  report
}

format_structural_qc_failure <- function(report) {
  counts <- report$rejection_counts
  counts <- counts[counts > 0L]
  summary <- if (length(counts)) {
    paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", ")
  } else {
    "none"
  }
  examples <- report$examples
  example_text <- vapply(names(examples), function(name) {
    rows <- examples[[name]]
    if (!length(rows)) return("")
    paste0(name, " rows ", paste(rows, collapse = ","))
  }, character(1))
  example_text <- example_text[nzchar(example_text)]
  if (length(example_text)) summary <- paste(summary, paste(example_text, collapse = "; "), sep = "; ")
  rejected <- report$rejected_rows %||% length(report$invalid_rows %||% integer())
  paste0("structural QC rejected ", as.integer(rejected), " row(s): ", summary)
}

apply_structural_qc <- function(data, input_build = "GRCh38", strict = FALSE,
                                row_policy = c("report", "error"),
                                require_statistics = TRUE, max_examples = 5L,
                                drop_duplicates = TRUE, check_duplicates = TRUE,
                                detail = c("full", "compact")) {
  row_policy <- match.arg(row_policy)
  detail <- match.arg(detail)
  if (isTRUE(strict)) row_policy <- "error"
  report <- structural_qc_report(data, input_build = input_build,
                                 require_statistics = require_statistics,
                                 max_examples = max_examples, detail = detail)
  if (identical(detail, "compact")) {
    internal <- report$internal %||% list(
      invalid = logical(nrow(data)), duplicate = logical(nrow(data)),
      canonical_key = character(nrow(data)), reason_count = integer(nrow(data))
    )
  } else {
    internal <- list(
      invalid = rep(FALSE, nrow(data)),
      duplicate = rep(FALSE, nrow(data)),
      canonical_key = report$canonical_key %||% rep(NA_character_, nrow(data)),
      reason_count = integer(nrow(data))
    )
    internal$invalid[report$invalid_rows %||% integer()] <- TRUE
    internal$duplicate[report$duplicate_rows %||% integer()] <- TRUE
  }
  invalid_mask <- internal$invalid
  duplicate_mask <- internal$duplicate
  invalid <- which(invalid_mask)
  if (!isTRUE(check_duplicates)) {
    invalid <- setdiff(invalid, which(duplicate_mask))
  }
  keep <- rep(TRUE, nrow(data))
  if (length(invalid)) {
    keep[invalid] <- FALSE
    # In report mode keep the first copy of an otherwise valid duplicate and
    # reject later copies. Strict/error mode still rejects every duplicate row.
    duplicate <- which(duplicate_mask)
    if (isTRUE(drop_duplicates) && identical(row_policy, "report") && length(duplicate)) {
      duplicate_key <- internal$canonical_key[duplicate]
      first <- !duplicated(duplicate_key)
      duplicate_first <- duplicate[first]
      only_duplicate_reason <- if (identical(detail, "compact")) {
        internal$reason_count[duplicate_first] == 1L
      } else {
        report$row_status$reasons[duplicate_first] == "duplicate_variant"
      }
      keep[duplicate_first[only_duplicate_reason]] <- TRUE
    }
  }
  if (identical(row_policy, "error") && length(invalid)) {
    failure <- format_structural_qc_failure(report)
    orientation_rejections <- report$rejection_counts[["orientation_mismatch"]] %||% 0L
    identity_rejections <- sum(as.numeric(report$rejection_counts[c(
      "missing_reference_allele", "missing_alternate_allele",
      "invalid_allele", "same_alleles"
    )]), na.rm = TRUE)
    if (orientation_rejections > 0L) {
      failure <- paste(
        "effect_allele and other_allele are inconsistent with explicit REF/ALT;",
        failure
      )
    } else if (identity_rejections > 0L) {
      failure <- paste(
        "explicit REF and ALT must be distinct single A/C/G/T alleles;",
        failure
      )
    }
    stop(failure, call. = FALSE)
  }
  if (identical(detail, "full")) report$row_status$accepted <- keep
  report$accepted_rows <- as.integer(sum(keep))
  report$rejected_rows <- as.integer(length(invalid))
  report$dropped_rows <- as.integer(sum(!keep))
  report$valid <- !length(invalid)
  report$internal <- NULL
  out <- data[keep, , drop = FALSE]
  attr(out, "structural_qc_report") <- report
  list(data = out, report = report)
}

compact_structural_qc_report <- function(report) {
  if (is.null(report)) return(NULL)
  compact <- report
  # Named atomic vectors are serialized by jsonlite as unnamed arrays. Keep
  # rejection reasons addressable in manifests by storing these aggregates as
  # named lists instead.
  compact$rejection_counts <- as.list(report$rejection_counts)
  compact$counts <- compact$rejection_counts
  compact$rejections <- as.list(report$rejections)
  # The public workflows retain aggregate counts and bounded examples in
  # manifests. Full row status/key vectors remain available from
  # preflight_sumstats() for callers that explicitly request an audit table,
  # but must not be carried through a large compression job.
  compact$row_status <- NULL
  compact$canonical_key <- NULL
  compact$invalid_rows <- NULL
  compact$duplicate_rows <- NULL
  compact$structurally_valid_rows <- NULL
  compact$internal <- NULL
  compact
}

validate_sumstats_values <- function(data, require_identity = TRUE) {
  if (isTRUE(require_identity) &&
      (anyNA(data$chromosome) || anyNA(data$base_pair_location) ||
       any(data$base_pair_location < 1L, na.rm = TRUE))) {
    stop("chromosome and positive base-pair locations are required", call. = FALSE)
  }
  if (anyNA(data$effect_allele) || anyNA(data$other_allele) ||
      any(data$effect_allele == "" | data$other_allele == "", na.rm = TRUE)) {
    stop("effect_allele and other_allele must be non-empty", call. = FALSE)
  }
  bad_eaf <- !is.na(data$effect_allele_frequency) & (data$effect_allele_frequency < 0 | data$effect_allele_frequency > 1)
  if (any(bad_eaf)) stop("effect_allele_frequency must be between 0 and 1", call. = FALSE)
  bad_beta <- !is.na(data$beta) & !is.finite(data$beta)
  if (any(bad_beta)) stop("beta must be finite when supplied", call. = FALSE)
  bad_se <- !is.na(data$standard_error) & (!is.finite(data$standard_error) | data$standard_error <= 0)
  if (any(bad_se)) stop("standard_error must be positive when supplied", call. = FALSE)
  bad_p <- !is.na(data$p_value) & (!is.finite(data$p_value) | data$p_value < 0 | data$p_value > 1)
  if (any(bad_p)) stop("p_value must be between 0 and 1 when supplied", call. = FALSE)
  bad_z <- !is.na(data$z) & !is.finite(data$z)
  if (any(bad_z)) stop("z must be finite when supplied", call. = FALSE)
  if ("sample_size" %in% names(data)) {
    bad_n <- !is.na(data$sample_size) & (!is.finite(data$sample_size) | data$sample_size <= 0)
    if (any(bad_n)) stop("sample_size must be positive when supplied", call. = FALSE)
  }
  if ("info" %in% names(data)) {
    bad_info <- !is.na(data$info) & !is.finite(data$info)
    if (any(bad_info)) stop("info must be finite when supplied", call. = FALSE)
  }
  invisible(data)
}

store_path <- function(store) {
  if (inherits(store, "compressor_store")) return(store$path)
  if (length(store) == 1L && is.character(store)) return(normalizePath(store, mustWork = FALSE))
  stop("store must be a compressor_store or a store directory path", call. = FALSE)
}

manifest_path <- function(store) file.path(store_path(store), "manifest.json")

write_manifest <- function(manifest, path) {
  jsonlite::write_json(
    manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17
  )
  invisible(path)
}

read_manifest <- function(path) {
  if (!file.exists(path)) stop("Missing CompreSSoR manifest: ", path, call. = FALSE)
  jsonlite::read_json(path, simplifyVector = FALSE)
}

now_utc <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

safe_file_size <- function(path) if (file.exists(path)) file.info(path)$size else NA_real_

stage_store_output <- function(output, overwrite = FALSE) {
  if (length(output) != 1L || !is.character(output) || is.na(output) || !nzchar(output)) {
    stop("output must be one directory path", call. = FALSE)
  }
  parent <- dirname(path.expand(output))
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(parent, mustWork = TRUE)
  target <- file.path(parent, basename(output))
  exists <- file.exists(target) || dir.exists(target)
  if (exists && !isTRUE(overwrite)) {
    stop("output already exists; use overwrite = TRUE", call. = FALSE)
  }
  staging <- tempfile(paste0(".", basename(output), "-staging-"), tmpdir = parent)
  if (!dir.create(staging, recursive = FALSE, showWarnings = FALSE)) {
    stop("could not create staging directory beside output", call. = FALSE)
  }
  list(target = target, staging = staging, overwrite = isTRUE(overwrite))
}

commit_store_output <- function(transaction) {
  target <- transaction$target
  staging <- transaction$staging
  backup <- NULL
  target_exists <- file.exists(target) || dir.exists(target)
  if (target_exists) {
    backup <- tempfile(paste0(".", basename(target), "-backup-"), tmpdir = dirname(target))
    if (!file.rename(target, backup)) {
      stop("could not move the existing output aside for atomic replacement", call. = FALSE)
    }
  }
  if (!file.rename(staging, target)) {
    if (!is.null(backup)) file.rename(backup, target)
    stop("could not move the completed store into place", call. = FALSE)
  }
  if (!is.null(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  invisible(target)
}
