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

clean_input_names <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  names(data) <- sub("^\\ufeff", "", trimws(names(data)))
  data
}

strict_fread <- function(...) {
  fatal <- NULL
  result <- withCallingHandlers(
    data.table::fread(...),
    warning = function(w) {
      message <- conditionMessage(w)
      if (grepl("stopped early|expected [0-9]+ fields|discarded single-line footer",
                message, ignore.case = TRUE)) {
        fatal <<- message
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!is.null(fatal)) stop("malformed delimited input: ", fatal, call. = FALSE)
  result
}

compressed_read_command <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  decompressor <- Sys.getenv("COMPRESSOR_DECOMPRESSOR", unset = "")
  if (!nzchar(decompressor)) decompressor <- unname(Sys.which("pigz"))
  if (length(decompressor) != 1L || !nzchar(decompressor)) decompressor <- "gzip"
  paste(shQuote(decompressor), "-dc", shQuote(path))
}

numeric_missing <- function(x) {
  text <- trimws(as.character(x))
  is.na(x) | is.na(text) | !nzchar(text) |
    tolower(text) %in% c(".", "na", "nan", "null")
}

parse_numeric_column <- function(x, name) {
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
    preview <- paste(utils::head(rows, 5L), collapse = ", ")
    stop("cannot parse ", name, " as numeric at row(s): ", preview,
         if (length(rows) > 5L) ", ..." else "", call. = FALSE)
  }
  out[missing] <- NA_real_
  out
}

parse_integer_column <- function(x, name) {
  value <- parse_numeric_column(x, name)
  bad <- !is.na(value) & (!is.finite(value) | value != trunc(value) |
                            value < -.Machine$integer.max | value > .Machine$integer.max)
  if (any(bad)) {
    rows <- which(bad)
    stop(name, " must contain whole 32-bit integers; invalid row(s): ",
         paste(utils::head(rows, 5L), collapse = ", "), call. = FALSE)
  }
  as.integer(value)
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

read_vcf_input <- function(input) {
  if (grepl("[.](?:gz|bgz)$", input, ignore.case = TRUE) && .Platform$OS.type != "windows") {
    command <- compressed_read_command(input)
    data <- strict_fread(cmd = command, skip = "#CHROM", data.table = FALSE,
                         showProgress = FALSE, check.names = FALSE)
  } else if (grepl("[.](?:gz|bgz)$", input, ignore.case = TRUE)) {
    temporary <- tempfile(fileext = ".vcf")
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    source <- gzfile(input, open = "rb")
    destination <- file(temporary, open = "wb")
    on.exit(try(close(source), silent = TRUE), add = TRUE)
    on.exit(try(close(destination), silent = TRUE), add = TRUE)
    repeat {
      bytes <- readBin(source, "raw", n = 1024^2)
      if (!length(bytes)) break
      writeBin(bytes, destination)
    }
    close(source)
    close(destination)
    data <- strict_fread(temporary, skip = "#CHROM", data.table = FALSE,
                         showProgress = FALSE, check.names = FALSE)
  } else {
    data <- strict_fread(input, skip = "#CHROM", data.table = FALSE,
                         showProgress = FALSE, check.names = FALSE)
  }
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
  if (any(invalid_or)) stop("VCF odds ratio must be finite and positive", call. = FALSE)
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

read_sumstats_input <- function(input) {
  if (is.data.frame(input)) return(clean_input_names(input))
  if (length(input) != 1L || !is.character(input) || !file.exists(input)) {
    stop("input must be a data.frame or an existing delimited file path", call. = FALSE)
  }
  if (grepl("[.]vcf(?:[.]bgz|[.]gz)?$", input, ignore.case = TRUE)) return(read_vcf_input(input))
  if (grepl("[.]gz$", input, ignore.case = TRUE)) {
    native <- tryCatch(
      strict_fread(input, data.table = FALSE, showProgress = FALSE),
      error = function(e) NULL
    )
    if (!is.null(native)) return(native)
    if (.Platform$OS.type == "windows") {
      stop("could not read compressed delimited input with data.table::fread", call. = FALSE)
    }
    command <- compressed_read_command(input)
    return(clean_input_names(strict_fread(cmd = command, data.table = FALSE, showProgress = FALSE)))
  }
  clean_input_names(strict_fread(input, data.table = FALSE, showProgress = FALSE))
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

normalise_sumstats_columns <- function(data) {
  data <- clean_input_names(data)
  data <- rename_first_alias(data, "rsid",
                             c("rsid", "rsids", "rsID", "RSID", "rs_id", "RS_ID", "dbsnp"))
  alias_map <- list(
    chromosome = c("chromosome", "chr", "CHR", "#chrom", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "bp"),
    effect_allele = c("effect_allele", "ea", "EA", "a1", "A1", "alt", "ALT"),
    other_allele = c("other_allele", "oa", "NEA", "nea", "a2", "A2", "ref", "REF"),
    beta = c("beta", "BETA", "b", "effect", "effect_size", "estimate",
             "ES", "LOGOR", "LOG_OR", "log_odds"),
    z = c("z", "Z", "zscore", "Z_SCORE", "z_stat", "ZSTAT", "zstatistic"),
    odds_ratio = c("odds_ratio", "OR", "ODDSRATIO", "oddsratio"),
    standard_error = c("standard_error", "SE", "se", "sebeta", "SEBETA", "stderr", "std_err"),
    effect_allele_frequency = c("effect_allele_frequency", "eaf", "EAF", "af", "AF",
                                "effect_af", "A1FREQ", "ALT_FREQ", "ALT_AF"),
    p_value = c("p_value", "P", "p", "pval", "pvalue"),
    minus_log10_p = c("minus_log10_p", "LP", "lp", "neglog10p", "-LOG10P"),
    variant_id = c("variant_id", "SNPID", "SNP_ID", "SNP", "snp", "variant", "ID"),
    sample_size = c("sample_size", "N", "n", "SS", "N_TOTAL", "TOTALSAMPLESIZE"),
    info = c("info", "INFO_SCORE", "IMPINFO", "INFO", "R2", "RSQ")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  missing <- setdiff(required_sumstats_columns(), names(data))
  if (length(missing)) {
    stop("Missing required summary-statistics columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!"chromosome" %in% names(data)) data$chromosome <- NA_character_
  if (!"base_pair_location" %in% names(data)) data$base_pair_location <- NA_integer_
  data$chromosome <- normalise_chromosome(data$chromosome)
  data$base_pair_location <- parse_integer_column(data$base_pair_location, "base_pair_location")
  data$effect_allele <- toupper(trimws(as.character(data$effect_allele)))
  data$other_allele <- toupper(trimws(as.character(data$other_allele)))
  if (!"beta" %in% names(data)) data$beta <- NA_real_
  if (!"standard_error" %in% names(data)) data$standard_error <- NA_real_
  if (!"z" %in% names(data)) data$z <- NA_real_
  if (!"odds_ratio" %in% names(data)) data$odds_ratio <- NA_real_
  data$beta <- parse_numeric_column(data$beta, "beta")
  data$standard_error <- parse_numeric_column(data$standard_error, "standard_error")
  data$z <- parse_numeric_column(data$z, "z")
  data$odds_ratio <- parse_numeric_column(data$odds_ratio, "odds_ratio")
  bad_or <- !is.na(data$odds_ratio) & (!is.finite(data$odds_ratio) | data$odds_ratio <= 0)
  if (any(bad_or)) stop("odds_ratio must be finite and positive when supplied", call. = FALSE)
  beta_missing <- is.na(data$beta) & is.finite(data$odds_ratio) & data$odds_ratio > 0
  data$beta[beta_missing] <- log(data$odds_ratio[beta_missing])
  if (!"effect_allele_frequency" %in% names(data)) data$effect_allele_frequency <- NA_real_
  data$effect_allele_frequency <- parse_numeric_column(data$effect_allele_frequency,
                                                        "effect_allele_frequency")
  if (!"p_value" %in% names(data)) data$p_value <- NA_real_
  data$p_value <- parse_numeric_column(data$p_value, "p_value")
  if ("minus_log10_p" %in% names(data)) {
    lp <- parse_numeric_column(data$minus_log10_p, "minus_log10_p")
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
    data$sample_size <- parse_numeric_column(data$sample_size, "sample_size")
  }
  if ("info" %in% names(data)) data$info <- parse_numeric_column(data$info, "info")
  alias_id <- !is.na(data$variant_id) & grepl("^rs", data$variant_id, ignore.case = TRUE)
  fill_rsid <- (is.na(data$rsid) | !nzchar(data$rsid)) & alias_id
  data$rsid[fill_rsid] <- data$variant_id[fill_rsid]
  data
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
