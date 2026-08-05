`%||%` <- function(x, y) if (is.null(x)) y else x

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
  # Coordinates can be recovered from the canonical rsID alias table. Alleles
  # remain required because they determine orientation and effect flipping.
  c("effect_allele", "other_allele")
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

read_delimited_file <- function(input, skip = NULL) {
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
  attr(out, "explicit_ref_alt") <- TRUE
  out
}

looks_like_vcf <- function(input) {
  if (grepl("[.]vcf(?:[.]bgz|[.]gz|[.]bz2|[.]xz|[.]zst)?$", input,
            ignore.case = TRUE)) return(TRUE)
  probe <- tryCatch(input_probe_lines(input, n = 16L), error = function(e) character())
  any(grepl("^##fileformat[=]VCF", probe, ignore.case = TRUE)) ||
    any(grepl("^#CHROM(?:\\t|\\s)", probe, ignore.case = TRUE))
}

read_sumstats_input <- function(input, parse_policy = c("error", "report")) {
  parse_policy <- match.arg(parse_policy)
  if (is.data.frame(input)) return(clean_input_names(input))
  if (length(input) != 1L || !is.character(input) || !file.exists(input)) {
    stop("input must be a data.frame or an existing delimited file path", call. = FALSE)
  }
  if (looks_like_vcf(input)) {
    return(read_vcf_input(input, parse_policy = parse_policy))
  }
  clean_input_names(read_delimited_file(input))
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

normalise_sumstats_columns <- function(data, parse_policy = c("error", "report")) {
  parse_policy <- match.arg(parse_policy)
  data <- clean_input_names(data)
  parse_failures <- list()
  assign_parsed <- function(target, value) {
    failures <- attr(value, "parse_failures") %||% integer()
    if (length(failures)) parse_failures[[target]] <<- as.integer(failures)
    data[[target]] <<- value
  }
  data <- rename_first_alias(data, "rsid",
                             c("rsid", "rsids", "rsID", "RSID", "rs_id", "RS_ID", "dbsnp"))
  alias_map <- list(
    chromosome = c("chromosome", "chr", "CHR", "#chrom", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "bp", "BP", "GENPOS"),
    effect_allele = c("effect_allele", "ea", "EA", "a1", "A1", "alt", "ALT", "ALLELE1"),
    other_allele = c("other_allele", "oa", "NEA", "nea", "a2", "A2", "ref", "REF", "ALLELE0"),
    beta = c("beta", "BETA", "b", "effect", "effect_size", "estimate",
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
    info = c("info", "INFO_SCORE", "IMPINFO", "INFO", "R2", "RSQ", "INFO_SCORE_IMP")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  missing <- setdiff(required_sumstats_columns(), names(data))
  if (length(missing)) {
    if (identical(parse_policy, "error")) {
      stop("Missing required summary-statistics columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    for (column in missing) data[[column]] <- NA_character_
  }
  if (!"chromosome" %in% names(data)) data$chromosome <- NA_character_
  if (!"base_pair_location" %in% names(data)) data$base_pair_location <- NA_integer_
  data$chromosome <- normalise_chromosome(data$chromosome)
  assign_parsed("base_pair_location",
                parse_integer_column(data$base_pair_location, "base_pair_location",
                                     invalid = parse_policy))
  data$effect_allele <- toupper(trimws(as.character(data$effect_allele)))
  data$other_allele <- toupper(trimws(as.character(data$other_allele)))
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
  beta_missing <- is.na(data$beta) & is.finite(data$odds_ratio) & data$odds_ratio > 0
  data$beta[beta_missing] <- log(data$odds_ratio[beta_missing])
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
  # Accept beta/SE, explicit Z/SE, or beta plus p (with the usual Wald
  # approximation). P is an input aid only; it is never written to a store.
  z_missing <- is.na(data$z) & is.finite(data$beta) & is.finite(data$standard_error) & data$standard_error > 0
  data$z[z_missing] <- data$beta[z_missing] / data$standard_error[z_missing]
  beta_missing_from_z <- is.na(data$beta) & is.finite(data$z) & is.finite(data$standard_error) & data$standard_error > 0
  data$beta[beta_missing_from_z] <- data$z[beta_missing_from_z] * data$standard_error[beta_missing_from_z]
  se_missing <- is.na(data$standard_error) & is.finite(data$beta) & is.finite(data$p_value) &
    data$p_value > 0 & data$p_value < 1 & data$beta != 0
  if (any(se_missing)) {
    data$standard_error[se_missing] <- abs(data$beta[se_missing]) /
      abs(stats::qnorm(data$p_value[se_missing] / 2))
  }
  z_missing <- is.na(data$z) & is.finite(data$beta) & is.finite(data$standard_error) & data$standard_error > 0
  data$z[z_missing] <- data$beta[z_missing] / data$standard_error[z_missing]
  consistent <- is.finite(data$beta) & is.finite(data$z) &
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
  if (!"rsid" %in% names(data)) data$rsid <- NA_character_
  if ("sample_size" %in% names(data)) {
    assign_parsed("sample_size", parse_numeric_column(data$sample_size, "sample_size",
                                                       invalid = parse_policy))
  }
  if ("info" %in% names(data)) {
    assign_parsed("info", parse_numeric_column(data$info, "info", invalid = parse_policy))
  }
  alias_id <- !is.na(data$variant_id) & grepl("^rs", data$variant_id, ignore.case = TRUE)
  fill_rsid <- (is.na(data$rsid) | !nzchar(data$rsid)) & alias_id
  data$rsid[fill_rsid] <- data$variant_id[fill_rsid]
  attr(data, "parse_failures") <- parse_failures
  attr(data, "missing_columns") <- missing
  data
}

new_structural_qc_report <- function(n, input_build, max_examples = 5L) {
  list(
    input_rows = as.integer(n),
    accepted_rows = NA_integer_,
    rejected_rows = NA_integer_,
    dropped_rows = NA_integer_,
    input_build = normalise_build_name(input_build),
    valid = FALSE,
    row_status = data.frame(
      row = seq_len(n), structurally_valid = logical(n), accepted = logical(n),
      reasons = rep("", n), stringsAsFactors = FALSE
    ),
    rejection_counts = integer(),
    examples = list(),
    max_examples = as.integer(max_examples)
  )
}

structural_qc_report <- function(data, input_build = "GRCh38",
                                 require_statistics = TRUE, max_examples = 5L) {
  data <- clean_input_names(data)
  n <- nrow(data)
  report <- new_structural_qc_report(n, input_build, max_examples = max_examples)
  report$missing_columns <- attr(data, "missing_columns") %||% character()
  if (!n) {
    report$accepted_rows <- 0L
    report$rejected_rows <- 0L
    report$dropped_rows <- 0L
    report$valid <- TRUE
    report$row_status$structurally_valid <- logical()
    report$row_status$accepted <- logical()
    return(report)
  }

  reasons <- list()
  add_reason <- function(name, mask) {
    mask <- as.logical(mask)
    mask[is.na(mask)] <- FALSE
    if (length(mask) != n) stop("structural QC mask has the wrong number of rows",
                                call. = FALSE)
    reasons[[name]] <<- mask
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

  add_reason("missing_effect_allele", is.na(effect) | !nzchar(effect))
  add_reason("missing_other_allele", is.na(other) | !nzchar(other))
  add_reason("multiallelic", (!is.na(effect) & grepl(",", effect, fixed = TRUE)) |
               (!is.na(other) & grepl(",", other, fixed = TRUE)))
  add_reason("symbolic_allele", (!is.na(effect) & grepl("^(?:<.*>|\\*)$", effect)) |
               (!is.na(other) & grepl("^(?:<.*>|\\*)$", other)))
  add_reason("indel", (!is.na(effect) & nchar(effect) != 1L) |
               (!is.na(other) & nchar(other) != 1L))
  add_reason("invalid_allele", (!is.na(effect) & !grepl("^[ACGT]$", effect)) |
               (!is.na(other) & !grepl("^[ACGT]$", other)))
  add_reason("same_alleles", !is.na(effect) & !is.na(other) & effect == other)

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
    add_reason("missing_statistics", !is.finite(z) | !is.finite(se) | se <= 0)
  }

  key <- ifelse(!is.na(chromosome) & is.finite(position) &
                  !is.na(effect) & !is.na(other),
                paste(chromosome, as.integer(position), other, effect, sep = ":"), NA_character_)
  duplicate <- !is.na(key) & (duplicated(key) | duplicated(key, fromLast = TRUE))
  add_reason("duplicate_variant", duplicate)
  report$canonical_key <- key

  reason_names <- names(reasons)
  reason_matrix <- if (length(reason_names)) {
    do.call(cbind, lapply(reasons, as.integer))
  } else {
    matrix(0L, nrow = n, ncol = 0L)
  }
  if (length(reason_names)) colnames(reason_matrix) <- reason_names
  invalid <- if (ncol(reason_matrix)) rowSums(reason_matrix) > 0L else rep(FALSE, n)
  reason_text <- if (length(reason_names)) {
    apply(reason_matrix, 1L, function(row) {
      paste(reason_names[which(row != 0L)], collapse = ",")
    })
  } else {
    rep("", n)
  }
  report$row_status$structurally_valid <- !invalid
  report$row_status$reasons <- reason_text
  report$rejection_counts <- sort(
    vapply(reasons, sum, integer(1)), decreasing = TRUE
  )
  report$counts <- report$rejection_counts
  report$rejections <- report$rejection_counts[report$rejection_counts > 0L]
  report$examples <- lapply(reasons, function(mask) {
    utils::head(which(mask), max(0L, as.integer(max_examples)))
  })
  report$invalid_rows <- which(invalid)
  report$duplicate_rows <- which(duplicate)
  report$structurally_valid_rows <- which(!invalid)
  report$valid <- !any(invalid)
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
  paste0("structural QC rejected ", length(report$invalid_rows), " row(s): ", summary)
}

apply_structural_qc <- function(data, input_build = "GRCh38", strict = FALSE,
                                row_policy = c("report", "error"),
                                require_statistics = TRUE, max_examples = 5L,
                                drop_duplicates = TRUE, check_duplicates = TRUE) {
  row_policy <- match.arg(row_policy)
  if (isTRUE(strict)) row_policy <- "error"
  report <- structural_qc_report(data, input_build = input_build,
                                  require_statistics = require_statistics,
                                  max_examples = max_examples)
  invalid <- report$invalid_rows
  if (!isTRUE(check_duplicates)) {
    invalid <- setdiff(invalid, report$duplicate_rows)
  }
  keep <- rep(TRUE, nrow(data))
  if (length(invalid)) {
    keep[invalid] <- FALSE
    # In report mode keep the first copy of an otherwise valid duplicate and
    # reject later copies. Strict/error mode still rejects every duplicate row.
    duplicate <- report$duplicate_rows
    if (isTRUE(drop_duplicates) && identical(row_policy, "report") && length(duplicate)) {
      duplicate_key <- report$canonical_key[duplicate]
      first <- !duplicated(duplicate_key)
      duplicate_first <- duplicate[first]
      only_duplicate_reason <- report$row_status$reasons[duplicate_first] == "duplicate_variant"
      keep[duplicate_first[only_duplicate_reason]] <- TRUE
    }
  }
  if (identical(row_policy, "error") && length(invalid)) {
    stop(format_structural_qc_failure(report), call. = FALSE)
  }
  report$row_status$accepted <- keep
  report$accepted_rows <- as.integer(sum(keep))
  report$rejected_rows <- as.integer(length(invalid))
  report$dropped_rows <- as.integer(sum(!keep))
  out <- data[keep, , drop = FALSE]
  attr(out, "structural_qc_report") <- report
  list(data = out, report = report)
}

compact_structural_qc_report <- function(report) {
  if (is.null(report)) return(NULL)
  compact <- report
  # The public workflows retain aggregate counts and bounded examples in
  # manifests. Full row status/key vectors remain available from
  # preflight_sumstats() for callers that explicitly request an audit table,
  # but must not be carried through a large compression job.
  compact$row_status <- NULL
  compact$canonical_key <- NULL
  compact$invalid_rows <- NULL
  compact$duplicate_rows <- NULL
  compact$structurally_valid_rows <- NULL
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
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}

read_manifest <- function(path) {
  if (!file.exists(path)) stop("Missing CompreSSoR manifest: ", path, call. = FALSE)
  jsonlite::read_json(path, simplifyVector = FALSE)
}

now_utc <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

safe_file_size <- function(path) if (file.exists(path)) file.info(path)$size else NA_real_

chain_manifest_metadata <- function(chain) {
  if (is.null(chain)) return(NULL)
  path <- normalizePath(chain, mustWork = TRUE)
  list(path = path, sha256 = digest::digest(path, algo = "sha256", file = TRUE))
}

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
