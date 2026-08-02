`%||%` <- function(x, y) if (is.null(x)) y else x

required_sumstats_columns <- function() {
  c("chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error")
}

clean_input_names <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  names(data) <- sub("^\\ufeff", "", trimws(names(data)))
  data
}

alias_key <- function(x) gsub("[^a-z0-9#]", "", tolower(as.character(x)), perl = TRUE)

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

read_vcf_input <- function(input) {
  con <- if (grepl("[.]gz$", input, ignore.case = TRUE)) gzfile(input, open = "rt") else file(input, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!grepl("^##", lines)]
  if (!length(lines)) stop("VCF contains no header or records", call. = FALSE)
  data <- data.table::fread(text = paste(lines, collapse = "\n"), data.table = FALSE,
                            showProgress = FALSE, check.names = FALSE)
  data <- clean_input_names(data)
  chrom_i <- first_alias_index(data, c("#CHROM", "CHROM", "chromosome"))
  pos_i <- first_alias_index(data, c("POS", "position", "base_pair_location"))
  ref_i <- first_alias_index(data, c("REF", "other_allele"))
  alt_i <- first_alias_index(data, c("ALT", "effect_allele"))
  id_i <- first_alias_index(data, c("ID", "variant_id", "SNP"))
  info_i <- first_alias_index(data, c("INFO", "info"))
  required <- c(chrom_i, pos_i, ref_i, alt_i)
  if (anyNA(required)) stop("VCF requires #CHROM, POS, REF and ALT columns", call. = FALSE)
  alt <- as.character(data[[alt_i]])
  if (any(grepl(",", alt, fixed = TRUE), na.rm = TRUE)) {
    stop("VCF contains multiallelic ALT values; split multiallelic records before conversion", call. = FALSE)
  }
  info <- if (!is.na(info_i)) parse_vcf_info(data[[info_i]]) else data.frame()
  field <- function(aliases) {
    direct_i <- first_alias_index(data, aliases)
    if (!is.na(direct_i)) return(as.character(data[[direct_i]]))
    hits <- intersect(alias_key(aliases), names(info))
    if (length(hits)) return(info[[hits[1L]]])
    rep(NA_character_, nrow(data))
  }
  beta <- field(c("BETA", "ES", "EFFECT", "LOGOR"))
  or <- field(c("OR", "ODDSRATIO"))
  beta_missing <- is.na(suppressWarnings(as.numeric(beta))) & is.finite(suppressWarnings(as.numeric(or)))
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
    standard_error = field(c("SE", "STDERR", "STDERR")),
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
  exclude <- names(data)[required]
  if (!is.na(id_i)) exclude <- c(exclude, names(data)[id_i])
  if (!is.na(info_i)) exclude <- c(exclude, names(data)[info_i])
  extra_names <- setdiff(names(data), exclude)
  if (length(extra_names)) out[extra_names] <- data[extra_names]
  out
}

read_sumstats_input <- function(input) {
  if (is.data.frame(input)) return(clean_input_names(input))
  if (length(input) != 1L || !is.character(input) || !file.exists(input)) {
    stop("input must be a data.frame or an existing delimited file path", call. = FALSE)
  }
  if (grepl("[.]vcf(?:[.]bgz|[.]gz)?$", input, ignore.case = TRUE)) return(read_vcf_input(input))
  if (grepl("[.]gz$", input, ignore.case = TRUE)) {
    if (.Platform$OS.type == "windows") {
      con <- gzfile(input, open = "rt")
      on.exit(close(con), add = TRUE)
      return(utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE))
    }
    native <- tryCatch(
      data.table::fread(input, data.table = FALSE, showProgress = FALSE),
      error = function(e) NULL
    )
    if (!is.null(native)) return(native)
    command <- paste("gzip -dc", shQuote(normalizePath(input)))
    return(clean_input_names(data.table::fread(cmd = command, data.table = FALSE, showProgress = FALSE)))
  }
  clean_input_names(data.table::fread(input, data.table = FALSE, showProgress = FALSE))
}

rename_first_alias <- function(data, target, aliases) {
  if (target %in% names(data)) return(data)
  hit <- first_alias_index(data, aliases)
  if (!is.na(hit)) names(data)[hit] <- target
  data
}

normalise_sumstats_columns <- function(data) {
  data <- clean_input_names(data)
  if (!"rsid" %in% names(data)) {
    rsid_aliases <- c("rsid", "rsids", "rsID", "RSID", "rs_id", "RS_ID")
    rsid_hit <- first_alias_index(data, rsid_aliases)
    if (!is.na(rsid_hit)) data$rsid <- as.character(data[[rsid_hit]])
  }
  alias_map <- list(
    chromosome = c("chromosome", "chr", "CHR", "#chrom", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "bp"),
    effect_allele = c("effect_allele", "ea", "EA", "a1", "A1", "alt", "ALT"),
    other_allele = c("other_allele", "oa", "NEA", "nea", "a2", "A2", "ref", "REF"),
    beta = c("beta", "BETA", "b", "effect", "effect_size", "estimate"),
    standard_error = c("standard_error", "SE", "se", "sebeta", "SEBETA", "stderr", "std_err"),
    effect_allele_frequency = c("effect_allele_frequency", "eaf", "EAF", "af", "AF", "effect_af"),
    p_value = c("p_value", "P", "p", "pval", "pvalue"),
    variant_id = c("variant_id", "SNPID", "SNP_ID", "SNP", "snp", "variant", "ID")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  missing <- setdiff(required_sumstats_columns(), names(data))
  if (length(missing)) {
    stop("Missing required summary-statistics columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data$chromosome <- sub("^chr", "", as.character(data$chromosome), ignore.case = TRUE)
  data$base_pair_location <- as.integer(data$base_pair_location)
  data$effect_allele <- toupper(as.character(data$effect_allele))
  data$other_allele <- toupper(as.character(data$other_allele))
  data$beta <- as.numeric(data$beta)
  data$standard_error <- as.numeric(data$standard_error)
  if (!"effect_allele_frequency" %in% names(data)) data$effect_allele_frequency <- NA_real_
  data$effect_allele_frequency <- as.numeric(data$effect_allele_frequency)
  if (!"p_value" %in% names(data)) data$p_value <- 2 * stats::pnorm(-abs(data$beta / data$standard_error))
  data$p_value <- as.numeric(data$p_value)
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
  data
}

validate_sumstats_values <- function(data) {
  if (anyNA(data$chromosome) || anyNA(data$base_pair_location) || any(data$base_pair_location < 1L, na.rm = TRUE)) {
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
