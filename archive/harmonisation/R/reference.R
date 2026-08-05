parse_canonical_variant_keys <- function(x) {
  x <- as.character(x)
  pieces <- strsplit(x, ":", fixed = TRUE)
  out <- data.frame(
    chromosome = rep(NA_character_, length(x)),
    base_pair_location = rep(NA_integer_, length(x)),
    reference_allele = rep(NA_character_, length(x)),
    alternate_allele = rep(NA_character_, length(x)),
    stringsAsFactors = FALSE
  )
  valid <- lengths(pieces) == 4L
  if (any(valid)) {
    valid_pieces <- pieces[valid]
    out$chromosome[valid] <- vapply(valid_pieces, function(z) z[[1L]], character(1))
    out$base_pair_location[valid] <- suppressWarnings(as.integer(
      vapply(valid_pieces, function(z) z[[2L]], character(1))
    ))
    out$reference_allele[valid] <- vapply(valid_pieces, function(z) z[[3L]], character(1))
    out$alternate_allele[valid] <- vapply(valid_pieces, function(z) z[[4L]], character(1))
  }
  out
}

compact_chromosome_labels <- function(x) {
  code <- suppressWarnings(as.integer(x))
  labels <- c(as.character(1:22), "X", "Y", "MT")
  out <- as.character(x)
  valid <- !is.na(code) & code >= 1L & code <= length(labels)
  out[valid] <- labels[code[valid]]
  out
}

compact_variant_type_labels <- function(x) {
  code <- suppressWarnings(as.integer(x))
  labels <- c("SNV", "MNV", "INDEL")
  out <- rep(NA_character_, length(code))
  valid <- !is.na(code) & code >= 1L & code <= length(labels)
  out[valid] <- labels[code[valid]]
  out
}

compact_position_labels <- function(data) {
  if (!"position_delta" %in% names(data) || !nrow(data)) return(data)
  if ("position_block" %in% names(data)) {
    block <- suppressWarnings(as.integer(data$position_block))
    if (anyNA(block)) stop("compact reference has missing position_block values", call. = FALSE)
    data <- data[order(block), , drop = FALSE]
    block <- block[order(block)]
  } else {
    block <- NULL
  }
  codes <- if ("chromosome_code" %in% names(data)) {
    suppressWarnings(as.integer(data$chromosome_code))
  } else {
    rep(NA_integer_, nrow(data))
  }
  delta <- suppressWarnings(as.integer(data$position_delta))
  position <- integer(length(delta))
  if (!is.null(block)) {
    starts <- c(1L, which(diff(block) != 0L) + 1L)
  } else {
    starts <- c(1L, which(!is.na(codes[-1L]) & !is.na(codes[-length(codes)]) &
                           codes[-1L] != codes[-length(codes)]) + 1L)
  }
  ends <- c(starts[-1L] - 1L, length(delta))
  for (i in seq_along(starts)) {
    span <- starts[i]:ends[i]
    position[span] <- cumsum(delta[span])
  }
  data$base_pair_location <- position
  data
}

compact_rsid_labels <- function(data) {
  if (!"rsid_code" %in% names(data) || "rsid" %in% names(data)) return(data)
  code <- data$rsid_code
  code_text <- as.character(code)
  rsid <- as.character(data$rsid_other %||% rep(NA_character_, length(code)))
  valid <- !is.na(code) & !is.na(code_text) & nzchar(code_text) & code_text != "NA"
  rsid[valid] <- paste0("rs", code_text[valid])
  data$rsid <- rsid
  data
}

normalise_reference_columns <- function(data) {
  source_metadata <- attr(data, "reference_source_metadata")
  source_format <- attr(data, "reference_source_format")
  pvar_metadata <- attr(data, "pvar_metadata")
  data <- clean_input_names(data)
  if (!"chromosome" %in% names(data) && "chromosome_code" %in% names(data)) {
    data$chromosome <- compact_chromosome_labels(data$chromosome_code)
  }
  if (!"variant_type" %in% names(data) && "variant_type_code" %in% names(data)) {
    data$variant_type <- compact_variant_type_labels(data$variant_type_code)
  }
  data <- compact_position_labels(data)
  data <- compact_rsid_labels(data)
  alias_map <- list(
    chromosome = c("chromosome", "chr", "CHR", "#chrom", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "bp", "BP"),
    reference_allele = c("reference_allele", "ref", "REF", "other_allele", "NEA", "nea", "a2", "A2"),
    alternate_allele = c("alternate_allele", "alt", "ALT", "effect_allele", "EA", "ea", "a1", "A1"),
    variant_index = c("variant_index", "reference_index", "vid", "VID", "index"),
    rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID", "dbsnp", "ID"),
    effect_allele_frequency = c("effect_allele_frequency", "eaf", "EAF", "af", "AF", "effect_af"),
    variant_id = c("variant_id", "variant_key", "canonical_variant")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])

  if (!"variant_id" %in% names(data)) {
    data$variant_id <- rep(NA_character_, nrow(data))
  }
  parsed <- parse_canonical_variant_keys(data$variant_id)
  for (field in names(parsed)) {
    if (!field %in% names(data)) {
      data[[field]] <- parsed[[field]]
    } else {
      missing_field <- is.na(data[[field]]) | !nzchar(as.character(data[[field]]))
      data[[field]][missing_field] <- parsed[[field]][missing_field]
    }
  }
  required <- c("chromosome", "base_pair_location", "reference_allele", "alternate_allele")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("canonical reference is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data$chromosome <- normalise_chromosome(data$chromosome)
  data$base_pair_location <- suppressWarnings(as.integer(data$base_pair_location))
  data$reference_allele <- toupper(trimws(as.character(data$reference_allele)))
  data$alternate_allele <- toupper(trimws(as.character(data$alternate_allele)))
  data$reference_allele[data$reference_allele %in% c("", ".", "NA")] <- NA_character_
  data$alternate_allele[data$alternate_allele %in% c("", ".", "NA")] <- NA_character_

  if (!"rsid" %in% names(data)) data$rsid <- rep(NA_character_, nrow(data))
  data$rsid <- as.character(data$rsid)
  data$rsid[data$rsid %in% c("", ".", "NA")] <- NA_character_
  if (!"effect_allele_frequency" %in% names(data)) {
    data$effect_allele_frequency <- rep(NA_real_, nrow(data))
  }
  data$effect_allele_frequency <- suppressWarnings(as.numeric(data$effect_allele_frequency))

  if (!"variant_type" %in% names(data)) {
    data$variant_type <- reference_variant_type(data$reference_allele,
                                                data$alternate_allele)
  }

  bad_eaf <- !is.na(data$effect_allele_frequency) &
    (data$effect_allele_frequency < 0 | data$effect_allele_frequency > 1)
  if (any(bad_eaf)) stop("reference effect_allele_frequency must be between 0 and 1", call. = FALSE)

  data$variant_key <- paste(data$chromosome, data$base_pair_location,
                            data$reference_allele, data$alternate_allele, sep = ":")
  bad_key <- is.na(data$chromosome) | is.na(data$base_pair_location) |
    is.na(data$reference_allele) | is.na(data$alternate_allele)
  data$variant_key[bad_key] <- NA_character_
  # variant_id is deliberately the canonical allele key. rsid is a secondary
  # annotation and is allowed to be duplicated or missing.
  data$variant_id <- if (length(data$variant_key)) data$variant_key else character()
  if (!is.null(source_metadata)) attr(data, "reference_source_metadata") <- source_metadata
  if (!is.null(source_format)) attr(data, "reference_source_format") <- source_format
  if (!is.null(pvar_metadata)) attr(data, "pvar_metadata") <- pvar_metadata
  data
}

complement_allele <- function(x) {
  x <- toupper(as.character(x))
  valid <- !is.na(x) & nchar(x) == 1L & x %in% c("A", "T", "C", "G")
  out <- rep(NA_character_, length(x))
  out[valid] <- chartr("ATCG", "TAGC", x[valid])
  out
}

reverse_complement_allele <- function(x) {
  x <- toupper(as.character(x))
  out <- rep(NA_character_, length(x))
  valid <- !is.na(x) & nzchar(x) & grepl("^[ACGTN]+$", x)
  out[valid] <- vapply(strsplit(x[valid], "", fixed = TRUE), function(chars) {
    paste(rev(chartr("ATCGN", "TAGCN", chars)), collapse = "")
  }, character(1))
  out
}

is_snp_allele <- function(x) !is.na(x) & nchar(x) == 1L & x %in% c("A", "C", "G", "T")

is_reference_allele <- function(x) {
  !is.na(x) & nzchar(x) & grepl("^[ACGTN]+$", x)
}

reference_variant_type <- function(reference_allele, alternate_allele) {
  snv <- is_snp_allele(reference_allele) & is_snp_allele(alternate_allele)
  mnv <- !snv & nchar(reference_allele) == nchar(alternate_allele)
  out <- ifelse(snv, "SNV", ifelse(mnv, "MNV", "INDEL"))
  out[is.na(reference_allele) | is.na(alternate_allele)] <- NA_character_
  out
}

split_reference_alternates <- function(data) {
  if (!nrow(data) || !"alternate_allele" %in% names(data)) return(data)
  alternate <- strsplit(as.character(data$alternate_allele), ",", fixed = TRUE)
  counts <- lengths(alternate)
  if (!any(counts > 1L)) return(data)
  rows <- rep(seq_len(nrow(data)), counts)
  out <- data[rows, , drop = FALSE]
  out$alternate_allele <- unlist(alternate, use.names = FALSE)
  row.names(out) <- NULL
  out
}

chromosome_sort_key <- function(x) {
  x <- toupper(sub("^chr", "", as.character(x), ignore.case = TRUE))
  numeric <- suppressWarnings(as.integer(x))
  numeric[is.na(numeric) & x == "X"] <- 23L
  numeric[is.na(numeric) & x == "Y"] <- 24L
  numeric[is.na(numeric) & x %in% c("M", "MT")] <- 25L
  numeric[is.na(numeric)] <- 1000L
  numeric
}

reference_normalized_cache_path <- function(resolved) {
  if (is.data.frame(resolved$variants)) return(NULL)
  path <- resolved$variants
  if (!is.character(path) || length(path) != 1L || grepl("[.]parquet$", path, ignore.case = TRUE)) {
    return(NULL)
  }
  sha256 <- resolved$metadata$sha256 %||% NULL
  if (is.null(sha256) || length(sha256) != 1L || !nzchar(sha256)) return(NULL)
  file.path(reference_cache_dir(resolved$cache_dir %||% NULL),
            paste0("normalized-", tolower(sha256), ".parquet"))
}

write_normalized_reference_cache <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("compressor-normalized-reference-", tmpdir = dirname(path), fileext = ".parquet")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  tryCatch({
    arrow::write_parquet(data, temporary, compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, use_dictionary = TRUE,
                         chunk_size = 65536L)
    if (file.exists(path)) unlink(path, force = TRUE)
    file.rename(temporary, path)
  }, error = function(e) NULL)
  invisible(path)
}

read_reference_vcf_source <- function(source) {
  # A canonical imputation VCF can contain tens of millions of rows and very
  # wide INFO/sample columns. Read only the five fixed VCF columns needed for
  # the reference, while leaving the full sum-statistics VCF reader unchanged.
  source <- normalizePath(source, mustWork = TRUE)
  command <- if (grepl("[.]gz$", source, ignore.case = TRUE)) {
    paste(compressed_read_command(source), "| grep -v '^##'")
  } else {
    paste("grep -v '^##'", shQuote(source))
  }
  out <- data.table::fread(cmd = command, data.table = FALSE, select = 1:5,
                           showProgress = FALSE, check.names = FALSE)
  if (ncol(out) < 5L) stop("VCF reference needs CHROM, POS, ID, REF and ALT columns", call. = FALSE)
  names(out)[seq_len(5L)] <- c("chromosome", "base_pair_location", "rsid",
                                "reference_allele", "alternate_allele")
  split_reference_alternates(out)
}

reference_source_format <- function(source) {
  if (!is.character(source) || length(source) != 1L) return("data.frame")
  if (grepl("[.]pvar$", source, ignore.case = TRUE)) return("pvar")
  if (grepl("[.]pvar[.]gz$", source, ignore.case = TRUE)) return("pvar.gz")
  if (grepl("[.]pvar[.]zst$", source, ignore.case = TRUE)) return("pvar.zst")
  if (grepl("[.]vcf(?:[.]gz|[.]bgz)?$", source, ignore.case = TRUE)) return("vcf")
  if (grepl("[.]bim(?:[.]gz)?$", source, ignore.case = TRUE)) return("bim")
  if (grepl("[.]parquet$", source, ignore.case = TRUE)) return("parquet")
  "delimited"
}

reference_source_command <- function(source) {
  source <- normalizePath(source, mustWork = TRUE)
  if (grepl("[.]gz$", source, ignore.case = TRUE)) {
    return(compressed_read_command(source))
  }
  if (grepl("[.]zst$", source, ignore.case = TRUE)) {
    zstd <- unname(Sys.which("zstd"))
    if (!nzchar(zstd)) {
      stop("reading a Zstandard reference requires the 'zstd' executable", call. = FALSE)
    }
    return(paste(shQuote(zstd), "-dc --quiet", shQuote(source)))
  }
  NULL
}

parse_reference_meta_value <- function(value) {
  value <- sub("^<", "", sub(">$", "", trimws(value)))
  pieces <- strsplit(value, ",(?=[A-Za-z_][A-Za-z0-9_]*=)", perl = TRUE)[[1L]]
  fields <- list()
  for (piece in pieces) {
    bits <- strsplit(piece, "=", fixed = TRUE)[[1L]]
    if (length(bits) < 2L) next
    key <- trimws(bits[1L])
    fields[[key]] <- paste(bits[-1L], collapse = "=")
  }
  fields
}

parse_reference_metadata_lines <- function(lines, format = "pvar") {
  lines <- sub("\\r$", "", as.character(lines))
  metadata_lines <- lines[startsWith(lines, "##")]
  header <- lines[startsWith(lines, "#CHROM") | startsWith(lines, "#chrom")]
  entries <- list()
  for (line in metadata_lines) {
    body <- sub("^##", "", line)
    key <- sub("=.*$", "", body)
    value <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else ""
    key <- trimws(key)
    if (!nzchar(key)) next
    parsed <- if (startsWith(trimws(value), "<")) {
      parse_reference_meta_value(value)
    } else {
      value
    }
    if (is.null(entries[[key]])) entries[[key]] <- list(parsed) else
      entries[[key]] <- c(entries[[key]], list(parsed))
  }
  columns <- if (length(header)) {
    strsplit(sub("^#", "", header[1L]), "[[:space:]]+", perl = TRUE)[[1L]]
  } else character()
  list(
    format = format,
    header = header[1L] %||% NA_character_,
    columns = columns,
    metadata = entries,
    raw = metadata_lines
  )
}

read_reference_metadata <- function(source, format = reference_source_format(source)) {
  if (!is.character(source) || length(source) != 1L || !file.exists(source)) return(NULL)
  if (!grepl("^pvar", format, ignore.case = TRUE)) return(NULL)
  command <- reference_source_command(source)
  con <- if (is.null(command)) file(normalizePath(source, mustWork = TRUE), open = "rt") else
    pipe(command, open = "rt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  lines <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    lines <- c(lines, line)
    if (startsWith(line, "#CHROM") || startsWith(line, "#chrom")) break
    if (!startsWith(line, "#")) break
  }
  parse_reference_metadata_lines(lines, format = format)
}

read_reference_pvar_source <- function(source) {
  format <- reference_source_format(source)
  command <- reference_source_command(source)
  args <- list(data.table = FALSE, skip = "#CHROM", showProgress = FALSE,
               check.names = FALSE)
  if (is.null(command)) args$input <- normalizePath(source, mustWork = TRUE) else
    args$cmd <- command
  out <- do.call(data.table::fread, args)
  if (ncol(out) < 5L) {
    stop("PVAR reference needs #CHROM, POS, ID, REF and ALT columns", call. = FALSE)
  }
  names(out)[seq_len(5L)] <- c("chromosome", "base_pair_location", "rsid",
                                "reference_allele", "alternate_allele")
  metadata <- read_reference_metadata(source, format)
  attr(out, "reference_source_metadata") <- metadata
  attr(out, "pvar_metadata") <- metadata
  attr(out, "reference_source_format") <- format
  split_reference_alternates(out)
}

read_reference_source <- function(source) {
  if (is.data.frame(source)) return(as.data.frame(source, stringsAsFactors = FALSE))
  if (length(source) != 1L || !is.character(source) || !file.exists(source)) {
    stop("reference source must be a data.frame or an existing reference file", call. = FALSE)
  }
  if (grepl("[.]parquet$", source, ignore.case = TRUE)) return(as.data.frame(arrow::read_parquet(source)))
  if (grepl("[.]bim(?:[.]gz)?$", source, ignore.case = TRUE)) {
    args <- list(data.table = FALSE, header = FALSE,
                 col.names = c("chromosome", "rsid", "genetic_distance",
                               "base_pair_location", "reference_allele", "alternate_allele"),
                 showProgress = FALSE)
    if (grepl("[.]gz$", source, ignore.case = TRUE)) {
      args$cmd <- compressed_read_command(source)
    } else {
      args$input <- source
    }
    return(do.call(data.table::fread, args))
  }
  if (grepl("[.]vcf(?:[.]gz|[.]bgz)?$", source, ignore.case = TRUE)) {
    return(read_reference_vcf_source(source))
  }
  if (grepl("[.]pvar(?:[.]gz|[.]zst)?$", source, ignore.case = TRUE)) {
    return(read_reference_pvar_source(source))
  }
  if (grepl("[.]gz$", source, ignore.case = TRUE)) {
    return(data.table::fread(cmd = compressed_read_command(source),
                             data.table = FALSE, showProgress = FALSE))
  }
  data.table::fread(source, data.table = FALSE, showProgress = FALSE)
}

lift_table_to_grch38 <- function(data, input_build, chain = NULL) {
  build <- toupper(gsub("[. -]", "", as.character(input_build %||% "GRCh38")))
  if (build %in% c("GRCH38", "HG38", "38")) {
    data$.compressor_liftover_status <- rep("not_needed", nrow(data))
    return(data)
  }
  if (!build %in% c("GRCH37", "HG19", "37")) {
    stop("input_build must be GRCh37/hg19 or GRCh38/hg38", call. = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("liftover requires the Bioconductor package 'rtracklayer'", call. = FALSE)
  }
  if (is.null(chain)) {
    stop("liftover from non-GRCh38 input requires a local GRCh37-to-GRCh38 chain file via chain", call. = FALSE)
  }
  if (!file.exists(chain)) stop("liftover chain file does not exist: ", chain, call. = FALSE)
  out <- data
  if ("variant_id" %in% names(out) && !".compressor_source_variant_id" %in% names(out)) {
    out$.compressor_source_variant_id <- as.character(out$variant_id)
  }
  out$.compressor_liftover_status <- "unmapped"
  out$.compressor_liftover_row <- seq_len(nrow(out))
  valid <- !is.na(data$chromosome) & nzchar(data$chromosome) &
    !is.na(data$base_pair_location) & data$base_pair_location >= 1L &
    !is.na(data$other_allele) & nzchar(data$other_allele)
  if (!any(valid)) return(out)

  source_rows <- which(valid)
  width <- pmax(1L, nchar(as.character(data$other_allele[valid])))
  width[is.na(width)] <- 1L
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", data$chromosome[valid]),
    ranges = IRanges::IRanges(start = data$base_pair_location[valid], width = width),
    # Source alleles are expressed on the forward genomic strand. Supplying
    # '+' here lets rtracklayer expose reverse-chain mappings as '-' so their
    # sequence alleles can be reverse-complemented correctly.
    strand = "+"
  )
  # rtracklayer versions in the supported Bioconductor range do not all
  # transparently decompress .chain.gz paths. Materialise a temporary
  # uncompressed chain only for the import step; the caller still supplies
  # and owns the original local chain file.
  chain_for_import <- normalizePath(chain, mustWork = TRUE)
  chain_tmp <- NULL
  if (grepl("[.]gz$", chain_for_import, ignore.case = TRUE)) {
    chain_tmp <- tempfile(fileext = ".chain")
    on.exit(unlink(chain_tmp), add = TRUE)
    input <- gzfile(chain_for_import, open = "rb")
    output <- file(chain_tmp, open = "wb")
    on.exit(try(close(input), silent = TRUE), add = TRUE)
    on.exit(try(close(output), silent = TRUE), add = TRUE)
    repeat {
      bytes <- readBin(input, what = "raw", n = 1024^2)
      if (!length(bytes)) break
      writeBin(bytes, output)
    }
    close(input)
    close(output)
    chain_for_import <- chain_tmp
  }
  chain_object <- rtracklayer::import.chain(chain_for_import)
  lifted <- rtracklayer::liftOver(gr, chain_object)
  lengths <- lengths(lifted)
  one <- lengths == 1L
  out$.compressor_liftover_status[source_rows] <- ifelse(
    lengths == 0L, "unmapped", ifelse(one, "mapped", "multi_mapped")
  )
  not_one <- !one
  out$chromosome[source_rows[not_one]] <- NA_character_
  out$base_pair_location[source_rows[not_one]] <- NA_integer_
  if (any(one)) {
    mapped <- unlist(lifted[one], use.names = FALSE)
    mapped_rows_all <- source_rows[one]
    mapped_chromosome <- sub(
      "^chr", "", as.character(GenomicRanges::seqnames(mapped)), ignore.case = TRUE
    )
    mapped_position <- GenomicRanges::start(mapped)
    chromosome_limit <- unname(
      compressor_grch38_chromosome_lengths[mapped_chromosome]
    )
    primary <- !is.na(chromosome_limit) & mapped_position >= 1L &
      mapped_position <= chromosome_limit
    if (any(!primary)) {
      rejected <- mapped_rows_all[!primary]
      out$.compressor_liftover_status[rejected] <- "non_primary_target"
      out$chromosome[rejected] <- NA_character_
      out$base_pair_location[rejected] <- NA_integer_
    }
    mapped_rows <- mapped_rows_all[primary]
    mapped <- mapped[primary]
    out$chromosome[mapped_rows] <- mapped_chromosome[primary]
    out$base_pair_location[mapped_rows] <- mapped_position[primary]
    reverse <- as.character(GenomicRanges::strand(mapped)) == "-"
    if (any(reverse)) {
      rows <- mapped_rows[reverse]
      old_effect <- out$effect_allele[rows]
      old_other <- out$other_allele[rows]
      out$effect_allele[rows] <- reverse_complement_allele(old_effect)
      out$other_allele[rows] <- reverse_complement_allele(old_other)
    }
    if ("variant_id" %in% names(out)) {
      out$variant_id[mapped_rows] <- paste(
        out$chromosome[mapped_rows], out$base_pair_location[mapped_rows],
        out$other_allele[mapped_rows], out$effect_allele[mapped_rows], sep = "_"
      )
    }
  }
  out
}

canonicalize_reference_data <- function(data, source_build = "GRCh38", chain = NULL,
                                        snps_only = FALSE) {
  source_metadata <- attr(data, "reference_source_metadata")
  source_format <- attr(data, "reference_source_format")
  pvar_metadata <- attr(data, "pvar_metadata")
  data <- split_reference_alternates(data)
  data <- normalise_reference_columns(data)
  # The liftover helper operates on summary-statistic allele names.
  data$effect_allele <- data$alternate_allele
  data$other_allele <- data$reference_allele
  data <- lift_table_to_grch38(data, source_build, chain = chain)
  data$reference_allele <- toupper(as.character(data$other_allele))
  data$alternate_allele <- toupper(as.character(data$effect_allele))
  valid <- is_reference_allele(data$reference_allele) &
    is_reference_allele(data$alternate_allele)
  if (isTRUE(snps_only)) {
    valid <- valid & is_snp_allele(data$reference_allele) &
      is_snp_allele(data$alternate_allele)
  }
  bad <- is.na(data$chromosome) | is.na(data$base_pair_location) |
    data$base_pair_location < 1L | !valid |
    data$reference_allele == data$alternate_allele
  data <- data[!bad, , drop = FALSE]
  data$chromosome <- sub("^chr", "", as.character(data$chromosome), ignore.case = TRUE)
  data$variant_key <- paste(data$chromosome, data$base_pair_location,
                            data$reference_allele, data$alternate_allele, sep = ":")

  # A VCF may carry more than one synonym for the same canonical allele. Keep
  # one dictionary row and preserve all non-empty IDs as a semicolon-delimited
  # secondary annotation.
  if (anyDuplicated(data$variant_key)) {
    groups <- split(seq_len(nrow(data)), data$variant_key)
    first <- vapply(groups, function(x) x[1L], integer(1))
    merged <- data[first, , drop = FALSE]
    merged$rsid <- vapply(groups, function(idx) {
      ids <- unique(unlist(strsplit(data$rsid[idx], ";", fixed = TRUE),
                           use.names = FALSE))
      ids <- ids[!is.na(ids) & nzchar(ids) & ids != "."]
      if (length(ids)) paste(ids, collapse = ";") else NA_character_
    }, character(1))
    data <- merged
  }
  data <- data[order(chromosome_sort_key(data$chromosome), data$base_pair_location,
                     data$reference_allele, data$alternate_allele), , drop = FALSE]
  row.names(data) <- NULL
  data$variant_index <- seq_len(nrow(data)) - 1L
  data$variant_id <- data$variant_key
  data$variant_type <- reference_variant_type(data$reference_allele,
                                               data$alternate_allele)
  if (!is.null(source_metadata)) attr(data, "reference_source_metadata") <- source_metadata
  if (!is.null(source_format)) attr(data, "reference_source_format") <- source_format
  if (!is.null(pvar_metadata)) attr(data, "pvar_metadata") <- pvar_metadata
  data
}

#' Build the immutable canonical GRCh38 variant reference
#'
#' The input is read once, normalised, sorted and written as one Parquet file.
#' The resulting variant_id is always chrom:pos:ref:alt; rsids are stored only
#' as secondary annotations. A study is never allowed to create this reference
#' implicitly.
#'
#' @param source A reference data.frame or VCF, PLINK BIM/PVAR, TSV or Parquet
#'   file. Compressed .gz and .zst PVAR sources are supported.
#' @param output Path of the canonical Parquet file to write.
#' @param source_build Build of source; non-GRCh38 inputs use the default
#'   local GRCh37-to-GRCh38 chain via chain.
#' @param chain Optional GRCh37-to-GRCh38 chain file used by rtracklayer.
#' @param snps_only Keep only single-base A/C/G/T variants. The default keeps
#'   SNVs, MNVs and indels with explicit sequence alleles.
#' @param overwrite Whether an existing output may be replaced.
#' @return The canonical reference descriptor, invisibly.
#' @export
build_canonical_reference <- function(source, output, source_build = "GRCh38", chain = NULL,
                                      snps_only = FALSE, overwrite = FALSE) {
  require_parquet_backend("build_canonical_reference()")
  if (length(output) != 1L || !is.character(output)) stop("output must be one file path", call. = FALSE)
  if (file.exists(output) && !isTRUE(overwrite)) stop("output already exists; use overwrite = TRUE", call. = FALSE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  source_data <- read_reference_source(source)
  source_metadata <- attr(source_data, "reference_source_metadata")
  source_format <- attr(source_data, "reference_source_format") %||%
    reference_source_format(source)
  data <- canonicalize_reference_data(source_data, source_build, chain,
                                      snps_only = snps_only)
  keep <- c("variant_index", "variant_id", "variant_key", "chromosome",
            "base_pair_location", "reference_allele", "alternate_allele", "rsid",
            "effect_allele_frequency", "variant_type")
  data <- data[intersect(keep, names(data))]
  # Keep a small, stable schema: the full reference is external to study stores.
  arrow::write_parquet(data, output, compression = "zstd", compression_level = 9,
                       write_statistics = TRUE, use_dictionary = TRUE, chunk_size = 65536L)
  sha256 <- reference_file_digest(output, "sha256")
  id <- if (isTRUE(snps_only)) "canonical_grch38_snp_v1" else "canonical_grch38_variant_v2"
  metadata <- list(
    id = id,
    build = "GRCh38",
    source = if (is.character(source)) normalizePath(source, mustWork = FALSE) else "in_memory",
    source_format = source_format,
    source_metadata = source_metadata,
    source_build = as.character(source_build),
    local_path = normalizePath(output, mustWork = FALSE),
    sha256 = sha256,
    rows = nrow(data),
    key = "chrom:pos:ref:alt",
    rsid = "secondary_annotation",
    snps_only = isTRUE(snps_only),
    variant_types = if (isTRUE(snps_only)) "SNV" else c("SNV", "MNV", "INDEL"),
    created_utc = now_utc()
  )
  write_manifest(metadata, sub("[.]parquet$", ".manifest.json", output, ignore.case = TRUE))
  descriptor <- list(id = metadata$id, build = metadata$build, source = metadata$source,
                     variants = output, cache_dir = dirname(output), metadata = metadata)
  structure(descriptor, class = "compressor_reference_descriptor")
}

empty_reference_table <- function() {
  normalise_reference_columns(data.frame(
    chromosome = character(), base_pair_location = integer(),
    reference_allele = character(), alternate_allele = character(),
    rsid = character(), variant_index = integer(),
    effect_allele_frequency = numeric(),
    stringsAsFactors = FALSE
  ))
}

reference_partition_paths <- function(path, chromosome = NULL) {
  master <- file.path(path, "variants.parquet")
  if (file.exists(master)) return(master)
  files <- list.files(path, pattern = "^chr[^/]+[.]parquet$",
                      full.names = TRUE, ignore.case = TRUE)
  fragment_order <- function(x) {
    basename_x <- basename(x)
    order(grepl("-part-[0-9]+[.]parquet$", basename_x),
          basename_x, na.last = TRUE)
  }
  if (is.null(chromosome)) return(files[fragment_order(files)])
  chromosome <- sub("^chr", "", as.character(chromosome), ignore.case = TRUE)
  pattern <- paste0("^chr", chromosome,
                    "(-part-[0-9]+)?[.]parquet$")
  hits <- files[grepl(pattern, basename(files), ignore.case = TRUE)]
  hits[fragment_order(hits)]
}

reference_partition_path <- function(path, chromosome) {
  paths <- reference_partition_paths(path, chromosome)
  if (length(paths)) paths[1L] else NULL
}

reference_partition_offset <- function(metadata, chromosome) {
  parts <- metadata$partitions %||% NULL
  if (is.null(parts)) return(0L)
  if (is.data.frame(parts) && all(c("chromosome", "first_index") %in% names(parts))) {
    hit <- as.character(parts$chromosome) == as.character(chromosome)
    if (any(hit) && !is.na(parts$first_index[which(hit)[1L]])) {
      return(as.integer(parts$first_index[which(hit)[1L]]))
    }
  }
  if (is.list(parts) && !is.data.frame(parts)) {
    hit <- vapply(parts, function(x) as.character(x$chromosome %||% "") ==
                    as.character(chromosome), logical(1L))
    if (any(hit)) {
      value <- parts[[which(hit)[1L]]]$first_index %||% NA_integer_
      if (!is.na(value)) return(as.integer(value))
    }
  }
  0L
}

reference_query_region <- function(chromosome = NULL, start = NULL, end = NULL,
                                   region = NULL) {
  if (!is.null(region)) {
    if (is.list(region)) {
      chromosome <- region$chromosome %||% region$chr %||% chromosome
      start <- region$start %||% region$begin %||% start
      end <- region$end %||% region$stop %||% end
    } else if (is.character(region) && length(region) == 1L) {
      match <- regexec("^(?:chr)?([^:]+):([0-9]+)(?:-([0-9]+))?$", region,
                       ignore.case = TRUE)
      pieces <- regmatches(region, match)[[1L]]
      if (!length(pieces)) stop("region must look like chr1:100-200", call. = FALSE)
      chromosome <- pieces[2L]
      start <- pieces[3L]
      end <- if (length(pieces) >= 4L && nzchar(pieces[4L])) pieces[4L] else pieces[3L]
    } else {
      stop("region must be a chr:start-end string or a list", call. = FALSE)
    }
  }
  if (!is.null(chromosome)) chromosome <- normalise_chromosome(chromosome)[1L]
  if (!is.null(start)) start <- parse_integer_column(start, "reference region start")
  if (!is.null(end)) end <- parse_integer_column(end, "reference region end")
  if (!is.null(start) && is.na(start)) stop("reference region start is missing", call. = FALSE)
  if (!is.null(end) && is.na(end)) stop("reference region end is missing", call. = FALSE)
  if (!is.null(start) && start < 1L) stop("reference region start must be positive", call. = FALSE)
  if (!is.null(end) && end < 1L) stop("reference region end must be positive", call. = FALSE)
  if (!is.null(start) && !is.null(end) && start > end) {
    stop("reference region start must not exceed end", call. = FALSE)
  }
  list(chromosome = chromosome, start = start, end = end)
}

reference_query_filter <- function(data, query, variant_ids = NULL, rsids = NULL,
                                   variant_indices = NULL) {
  keep <- rep(TRUE, nrow(data))
  if (!is.null(query$chromosome)) keep <- keep & data$chromosome == query$chromosome
  if (!is.null(query$start)) keep <- keep & data$base_pair_location >= query$start
  if (!is.null(query$end)) keep <- keep & data$base_pair_location <= query$end
  if (!is.null(variant_ids)) {
    values <- unique(as.character(variant_ids))
    rsid_match <- vapply(strsplit(as.character(data$rsid), ";", fixed = TRUE),
                         function(value) any(value %in% values), logical(1L))
    keep <- keep & (data$variant_id %in% values | rsid_match)
  }
  if (!is.null(rsids)) {
    values <- unique(as.character(rsids))
    rsid_match <- vapply(strsplit(as.character(data$rsid), ";", fixed = TRUE),
                         function(value) any(value %in% values), logical(1L))
    keep <- keep & rsid_match
  }
  if (!is.null(variant_indices)) {
    values <- suppressWarnings(as.integer(variant_indices))
    keep <- keep & data$variant_index %in% values
  }
  keep[is.na(keep)] <- FALSE
  keep
}

reference_parquet_query <- function(path, query, variant_ids = NULL, rsids = NULL,
                                    variant_indices = NULL, compact = FALSE,
                                    index_offset = 0L) {
  reader <- arrow::ParquetFileReader$create(path)
  groups <- seq_len(reader$num_row_groups) - 1L
  needs_filter <- !is.null(query$chromosome) || !is.null(query$start) ||
    !is.null(query$end) || !is.null(variant_ids) || !is.null(rsids) ||
    !is.null(variant_indices)
  if (!needs_filter && length(groups)) {
    raw <- as.data.frame(reader$ReadRowGroups(groups))
    out <- if (nrow(raw)) normalise_reference_columns(raw) else empty_reference_table()
    if (!"variant_index" %in% names(out)) {
      out$variant_index <- as.integer(index_offset + seq_len(nrow(out)) - 1L)
    }
    return(out)
  }
  if (!length(groups)) return(empty_reference_table())

  parts <- list()
  row_offset <- as.integer(index_offset)
  for (group in groups) {
    raw <- as.data.frame(reader$ReadRowGroup(group))
    if (!nrow(raw)) next
    compact_block <- if (compact && "position_block" %in% names(raw)) {
      unique(raw$position_block)
    } else NULL
    part <- normalise_reference_columns(raw)
    if (!"variant_index" %in% names(part)) {
      part$variant_index <- as.integer(row_offset + seq_len(nrow(part)) - 1L)
    }
    row_offset <- row_offset + nrow(part)
    # A compact master row group is one delta block. Decode it in full, then
    # retain only the requested rows; no complete chromosome is accumulated.
    if (!is.null(compact_block) && length(compact_block) > 1L) {
      stop("compact reference row group crosses position blocks", call. = FALSE)
    }
    part <- part[reference_query_filter(part, query, variant_ids, rsids,
                                        variant_indices), , drop = FALSE]
    if (nrow(part)) parts[[length(parts) + 1L]] <- part
  }
  if (!length(parts)) empty_reference_table() else do.call(rbind, parts)
}

reference_table <- function(reference, chromosome = NULL, start = NULL, end = NULL,
                            region = NULL, variant_ids = NULL, rsids = NULL,
                            variant_indices = NULL, columns = NULL) {
  if (is.null(reference)) return(NULL)
  query <- reference_query_region(chromosome, start, end, region)
  if (!is.null(variant_ids)) {
    parsed <- parse_canonical_variant_keys(variant_ids)
    parsed_chromosome <- unique(parsed$chromosome[!is.na(parsed$chromosome)])
    if (is.null(query$chromosome) && length(parsed_chromosome) == 1L) {
      query$chromosome <- normalise_chromosome(parsed_chromosome)
    }
    parsed_position <- parsed$base_pair_location[!is.na(parsed$base_pair_location)]
    if (length(parsed_position) && is.null(query$start)) query$start <- min(parsed_position)
    if (length(parsed_position) && is.null(query$end)) query$end <- max(parsed_position)
  }
  resolved <- resolve_reference(reference)
  metadata <- resolved$metadata %||% list(
    id = resolved$id %||% "local", build = resolved$build %||% "GRCh38"
  )
  if (is.data.frame(resolved$variants)) {
    out <- if (!nrow(resolved$variants)) empty_reference_table() else
      normalise_reference_columns(resolved$variants)
    if (!"variant_index" %in% names(out)) {
      out$variant_index <- as.integer(seq_len(nrow(out)) - 1L)
    }
    out <- out[reference_query_filter(out, query, variant_ids, rsids,
                                      variant_indices), , drop = FALSE]
  } else {
    path <- resolved$variants
    if (dir.exists(path)) {
      master <- file.path(path, "variants.parquet")
      if (file.exists(master)) {
        out <- reference_parquet_query(master, query, variant_ids, rsids,
                                       variant_indices, compact = TRUE)
      } else {
        paths <- reference_partition_paths(path, query$chromosome)
        if (!length(paths)) {
          out <- empty_reference_table()
        } else {
          running_offsets <- vapply(paths, function(file) {
            chromosome <- sub("^chr", "", sub("[-]part-[0-9]+$", "",
                                               tools::file_path_sans_ext(basename(file))),
                               ignore.case = TRUE)
            reference_partition_offset(metadata, chromosome)
          }, integer(1L))
          if (length(paths) > 1L) {
            for (i in seq_len(length(paths) - 1L)) {
              rows <- tryCatch(arrow::ParquetFileReader$create(paths[i])$num_rows,
                               error = function(e) 0L)
              running_offsets[i + 1L] <- running_offsets[i] + as.integer(rows)
            }
          }
          parts <- lapply(paths, function(file) {
            offset <- running_offsets[match(file, paths)]
            reference_parquet_query(file, query, variant_ids, rsids,
                                    variant_indices, index_offset = offset)
          })
          parts <- Filter(function(x) nrow(x), parts)
          out <- if (length(parts)) do.call(rbind, parts) else empty_reference_table()
        }
      }
    } else if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
      stop("resolved reference does not contain a readable variants asset", call. = FALSE)
    } else if (grepl("[.]parquet$", path, ignore.case = TRUE)) {
      out <- reference_parquet_query(path, query, variant_ids, rsids,
                                     variant_indices)
    } else if (grepl("[.]bim(?:[.]gz)?$", path, ignore.case = TRUE)) {
      out <- normalise_reference_columns(read_reference_source(path))
      if (!"variant_index" %in% names(out)) {
        out$variant_index <- as.integer(seq_len(nrow(out)) - 1L)
      }
      out <- out[reference_query_filter(out, query, variant_ids, rsids,
                                        variant_indices), , drop = FALSE]
    } else {
      normalized_cache <- reference_normalized_cache_path(resolved)
      if (!is.null(normalized_cache) && file.exists(normalized_cache)) {
        out <- reference_parquet_query(normalized_cache, query, variant_ids, rsids,
                                       variant_indices)
      } else {
        out <- normalise_reference_columns(read_reference_source(path))
        if (!"variant_index" %in% names(out)) {
          out$variant_index <- as.integer(seq_len(nrow(out)) - 1L)
        }
        if (!is.null(normalized_cache)) write_normalized_reference_cache(out, normalized_cache)
        out <- out[reference_query_filter(out, query, variant_ids, rsids,
                                          variant_indices), , drop = FALSE]
      }
    }
  }
  row.names(out) <- NULL
  if (anyNA(out$variant_key) || anyDuplicated(out$variant_key)) {
    stop("canonical reference has missing or duplicated chrom:pos:ref:alt keys", call. = FALSE)
  }
  if (!"variant_index" %in% names(out)) {
    offset <- if (is.null(query$chromosome)) 0L else
      reference_partition_offset(metadata, query$chromosome)
    out$variant_index <- as.integer(offset + seq_len(nrow(out)) - 1L)
  }
  out$variant_index <- as.integer(out$variant_index)
  if (anyDuplicated(out$variant_index) || any(out$variant_index < 0L, na.rm = TRUE)) {
    stop("canonical reference has invalid or duplicated variant_index values", call. = FALSE)
  }
  out$reference_index <- out$variant_index
  if (!is.null(query$chromosome)) metadata$partition_chromosome <- as.character(query$chromosome)
  normalized_cache <- reference_normalized_cache_path(resolved)
  if (!is.null(normalized_cache) && file.exists(normalized_cache)) {
    metadata$normalized_cache_path <- normalizePath(normalized_cache, mustWork = FALSE)
  }
  source_metadata <- attr(out, "reference_source_metadata") %||%
    resolved$metadata$source_metadata %||% NULL
  source_format <- attr(out, "reference_source_format") %||%
    resolved$metadata$source_format %||% NULL
  if (!is.null(source_metadata)) metadata$source_metadata <- source_metadata
  if (!is.null(source_format)) metadata$source_format <- source_format
  attr(out, "reference_metadata") <- metadata
  if (!is.null(metadata$sha256)) attr(out, "reference_sha256") <- metadata$sha256
  if (!is.null(columns)) {
    columns <- unique(as.character(columns))
    missing_columns <- setdiff(columns, names(out))
    if (length(missing_columns)) {
      stop("reference query requested missing columns: ",
           paste(missing_columns, collapse = ", "), call. = FALSE)
    }
    out <- out[columns, drop = FALSE]
    attr(out, "reference_metadata") <- metadata
    if (!is.null(metadata$sha256)) attr(out, "reference_sha256") <- metadata$sha256
  }
  out
}

reference_hash <- function(reference_data) {
  if (is.null(reference_data)) return(NA_character_)
  known_hash <- attr(reference_data, "reference_sha256")
  if (!is.null(known_hash) && length(known_hash) == 1L && nzchar(known_hash)) return(known_hash)
  digest::digest(reference_data, algo = "sha256", serialize = TRUE)
}

.compressor_harmonise_context <- new.env(parent = emptyenv())

normalise_reference_qc <- function(qc = NULL) {
  defaults <- list(
    liftover = TRUE,
    strand = TRUE,
    palindromic = "drop",
    frequency = "none",
    frequency_tolerance = 0.1,
    example_limit = 5L
  )
  if (is.null(qc)) return(defaults)
  if (is.logical(qc) && length(qc) == 1L && !is.na(qc)) {
    qc <- list(strand = qc)
  }
  if (!is.list(qc)) stop("qc must be NULL, logical, or a named list", call. = FALSE)
  aliases <- c(
    frequency_qc = "frequency", eaf = "frequency", eaf_qc = "frequency",
    frequency_check = "frequency", allow_strand = "strand",
    strand_qc = "strand", liftover_qc = "liftover",
    examples = "example_limit", max_examples = "example_limit"
  )
  for (name in names(aliases)) {
    if (!is.null(qc[[name]]) && is.null(qc[[aliases[[name]]]])) {
      qc[[aliases[[name]]]] <- qc[[name]]
    }
  }
  result <- utils::modifyList(defaults, qc[intersect(names(qc), names(defaults))])
  if (length(result$liftover) != 1L || !is.logical(result$liftover) ||
      is.na(result$liftover)) stop("qc$liftover must be TRUE or FALSE", call. = FALSE)
  if (length(result$strand) != 1L || !is.logical(result$strand) || is.na(result$strand)) {
    stop("qc$strand must be TRUE or FALSE", call. = FALSE)
  }
  palindromic <- tolower(as.character(result$palindromic))
  palindromic[palindromic %in% c("keep", "retain", "direct")] <- "allow"
  palindromic[palindromic %in% c("infer", "eaf")] <- "frequency"
  result$palindromic <- match.arg(palindromic, c("drop", "frequency", "allow"))
  frequency <- result$frequency
  if (is.logical(frequency) && length(frequency) == 1L && !is.na(frequency)) {
    frequency <- if (frequency) "report" else "none"
  }
  frequency <- tolower(as.character(frequency))
  frequency[frequency %in% c("check", "warn", "warning")] <- "report"
  frequency[frequency %in% c("fail", "stop")] <- "error"
  result$frequency <- match.arg(frequency, c("none", "report", "drop", "error"))
  if (length(result$frequency_tolerance) != 1L ||
      !is.numeric(result$frequency_tolerance) ||
      is.na(result$frequency_tolerance) || result$frequency_tolerance < 0 ||
      result$frequency_tolerance > 1) {
    stop("qc$frequency_tolerance must be between 0 and 1", call. = FALSE)
  }
  if (length(result$example_limit) != 1L || !is.numeric(result$example_limit) ||
      is.na(result$example_limit) || result$example_limit < 0 ||
      result$example_limit != as.integer(result$example_limit)) {
    stop("qc$example_limit must be a non-negative integer", call. = FALSE)
  }
  result$example_limit <- as.integer(result$example_limit)
  result
}

reference_palindromic_pair <- function(effect, other) {
  is_snp_allele(effect) & is_snp_allele(other) &
    ((effect == "A" & other == "T") | (effect == "T" & other == "A") |
       (effect == "C" & other == "G") | (effect == "G" & other == "C"))
}

alignment_example_rows <- function(data, status, limit = 5L, extra_status = NULL) {
  if (!nrow(data) || !limit) return(list())
  status <- as.character(status)
  if (!is.null(extra_status)) status <- as.character(extra_status)
  wanted <- unique(status[!is.na(status) & status != "aligned"])
  out <- list()
  for (label in wanted) {
    rows <- which(status == label)[seq_len(min(length(which(status == label)), limit))]
    out[[label]] <- lapply(rows, function(i) {
      list(
        row = as.integer(data$.compressor_input_order[i] %||% i),
        variant_id = as.character(data$variant_id[i] %||% NA_character_),
        chromosome = as.character(data$chromosome[i] %||% NA_character_),
        base_pair_location = as.integer(data$base_pair_location[i] %||% NA_integer_),
        effect_allele = as.character(data$effect_allele[i] %||% NA_character_),
        other_allele = as.character(data$other_allele[i] %||% NA_character_)
      )
    })
  }
  out
}

alignment_diagnostics <- function(data, status, counts, limit = 5L,
                                  frequency_status = NULL) {
  examples <- alignment_example_rows(data, status, limit)
  if (!is.null(frequency_status)) {
    frequency_rows <- which(frequency_status == "mismatch")
    if (length(frequency_rows)) {
      examples$frequency_mismatch <- lapply(
        utils::head(frequency_rows, limit), function(i) {
          list(
            row = as.integer(data$.compressor_input_order[i] %||% i),
            variant_id = as.character(data$variant_id[i] %||% NA_character_),
            chromosome = as.character(data$chromosome[i] %||% NA_character_),
            base_pair_location = as.integer(data$base_pair_location[i] %||% NA_integer_),
            effect_allele = as.character(data$effect_allele[i] %||% NA_character_),
            other_allele = as.character(data$other_allele[i] %||% NA_character_)
          )
        }
      )
    }
  }
  list(
    counts = counts,
    examples = examples
  )
}

harmonise_to_reference <- function(data, reference, strict = FALSE, preserve = FALSE,
                                   reference_data = NULL, reference_key = NULL) {
  qc <- .compressor_harmonise_context$qc %||% normalise_reference_qc()
  input_rows <- nrow(data)
  ref <- if (!is.null(reference_data)) reference_data else reference_table(reference)
  if (is.null(ref)) {
    data$harmonisation_status <- rep("unreferenced", input_rows)
    data$harmonisation_flip <- rep(FALSE, input_rows)
    data$input_duplicate <- duplicated(data$variant_id)
    counts <- list(
      input_rows = input_rows, output_rows = nrow(data), aligned_rows = 0L,
      duplicate_input_rows = sum(data$input_duplicate), unmatched_rows = 0L,
      incompatible_rows = 0L, ambiguous_rows = 0L, dropped_duplicates = 0L,
      dropped_unmatched = 0L, dropped_incompatible = 0L, dropped_ambiguous = 0L
    )
    counts$diagnostic_counts <- counts
    counts$diagnostic_examples <- list()
    counts$counts <- counts$diagnostic_counts
    counts$examples <- counts$diagnostic_examples
    counts$diagnostics <- alignment_diagnostics(data, data$harmonisation_status,
                                                counts$diagnostic_counts, qc$example_limit)
    return(list(data = data, reference_hash = NA_character_,
                reference_rows = NA_integer_, reference_metadata = NULL,
                alignment_stats = counts))
  }
  input_key <- paste(data$chromosome, data$base_pair_location,
                     data$other_allele, data$effect_allele, sep = ":")
  input_key[is.na(data$chromosome) | is.na(data$base_pair_location) |
              is.na(data$other_allele) | is.na(data$effect_allele)] <- NA_character_
  input_effect_complement <- complement_allele(data$effect_allele)
  input_other_complement <- complement_allele(data$other_allele)
  complement_key <- paste(data$chromosome, data$base_pair_location,
                          input_other_complement, input_effect_complement, sep = ":")
  reverse_key <- paste(data$chromosome, data$base_pair_location,
                       data$effect_allele, data$other_allele, sep = ":")
  complement_reverse_key <- paste(data$chromosome, data$base_pair_location,
                                  input_effect_complement, input_other_complement, sep = ":")
  ref_key <- reference_key %||% ref$variant_key
  direct_i <- match(input_key, ref_key)
  reverse_i <- match(reverse_key, ref_key)
  complement_i <- if (isTRUE(qc$strand)) match(complement_key, ref_key) else
    rep(NA_integer_, input_rows)
  complement_reverse_i <- if (isTRUE(qc$strand)) match(complement_reverse_key, ref_key) else
    rep(NA_integer_, input_rows)
  candidates <- cbind(!is.na(direct_i), !is.na(reverse_i),
                      !is.na(complement_i), !is.na(complement_reverse_i))
  candidate_count <- rowSums(candidates)
  unmatched <- candidate_count == 0L
  multi <- candidate_count > 1L
  key <- direct_i
  selected <- rep("direct", input_rows)
  replace <- is.na(key) & !is.na(reverse_i)
  key[replace] <- reverse_i[replace]
  selected[replace] <- "reverse"
  replace <- is.na(key) & !is.na(complement_i)
  key[replace] <- complement_i[replace]
  selected[replace] <- "complement"
  replace <- is.na(key) & !is.na(complement_reverse_i)
  key[replace] <- complement_reverse_i[replace]
  selected[replace] <- "complement_reverse"
  palindromic <- reference_palindromic_pair(data$effect_allele, data$other_allele)
  ambiguous <- multi
  if (qc$palindromic != "drop") {
    pal_rows <- palindromic & !unmatched
    # Direct and swapped matches are allele-identifiable; only a strand-only
    # palindromic match needs frequency evidence. Prefer the non-strand path.
    direct_path <- pal_rows & !is.na(direct_i)
    reverse_path <- pal_rows & is.na(direct_i) & !is.na(reverse_i)
    key[direct_path] <- direct_i[direct_path]
    selected[direct_path] <- "direct"
    key[reverse_path] <- reverse_i[reverse_path]
    selected[reverse_path] <- "reverse"
    ambiguous[pal_rows] <- FALSE
    if (qc$palindromic == "frequency") {
      ref_eaf <- rep(NA_real_, input_rows)
      has_key <- !is.na(key)
      ref_eaf[has_key] <- ref$effect_allele_frequency[key[has_key]]
      study_eaf <- data$effect_allele_frequency
      direct_distance <- abs(study_eaf - ref_eaf)
      flipped_distance <- abs((1 - study_eaf) - ref_eaf)
      evidence <- pal_rows & is.finite(study_eaf) & is.finite(ref_eaf)
      choose_flip <- evidence & flipped_distance < direct_distance
      choose_direct <- evidence & direct_distance < flipped_distance
      ambiguous[pal_rows & !evidence] <- TRUE
      ambiguous[pal_rows & evidence & !choose_flip & !choose_direct] <- TRUE
      selected[choose_flip] <- "reverse"
      # The canonical key is the same for a palindromic pair, so retaining the
      # selected path is enough for the later flip calculation.
    }
  }
  incompatible <- !unmatched & !multi & !is.na(key) &
    !(is_reference_allele(data$effect_allele) & is_reference_allele(data$other_allele))
  # Palindromic strand matches cannot be proven from alleles. Population EAF
  # differences are not strand evidence, so ambiguous A/T and C/G rows remain
  # ambiguous and are dropped rather than risking an effect inversion.
  status <- rep("aligned", input_rows)
  status[unmatched] <- "unmatched"
  status[incompatible] <- "incompatible"
  status[ambiguous] <- "ambiguous"
  flip_flag <- selected %in% c("reverse", "complement_reverse") & !is.na(key)
  flip_flag[unmatched | incompatible | ambiguous] <- FALSE
  aligned <- !unmatched & !incompatible & !ambiguous

  # Any second input row targeting the same canonical variant is discarded as
  # ambiguous. Keeping neither row avoids silently choosing a study duplicate.
  duplicate_target <- rep(FALSE, input_rows)
  duplicate_target[aligned] <- duplicated(key[aligned]) | duplicated(key[aligned], fromLast = TRUE)
  status[duplicate_target] <- "duplicate"
  aligned[duplicate_target] <- FALSE
  frequency_status <- rep("not_checked", input_rows)
  if (qc$frequency != "none") {
    ref_eaf <- rep(NA_real_, input_rows)
    has_key <- !is.na(key)
    ref_eaf[has_key] <- ref$effect_allele_frequency[key[has_key]]
    study_eaf <- data$effect_allele_frequency
    comparable <- aligned & is.finite(study_eaf) & is.finite(ref_eaf)
    frequency_status[aligned & !is.finite(study_eaf)] <- "missing_study"
    frequency_status[aligned & is.finite(study_eaf) & !is.finite(ref_eaf)] <- "missing_reference"
    oriented_eaf <- ifelse(flip_flag, 1 - study_eaf, study_eaf)
    frequency_status[comparable] <- ifelse(
      abs(oriented_eaf[comparable] - ref_eaf[comparable]) <= qc$frequency_tolerance,
      "match", "mismatch"
    )
    mismatch <- frequency_status == "mismatch"
    if (qc$frequency %in% c("drop", "error")) {
      status[mismatch] <- "frequency_mismatch"
      aligned[mismatch] <- FALSE
    }
    if (qc$frequency == "error" && any(mismatch)) {
      stop("reference frequency QC failed: mismatched=", sum(mismatch), call. = FALSE)
    }
  }
  if (isTRUE(strict) && any(unmatched | incompatible | ambiguous | duplicate_target |
                             status == "frequency_mismatch")) {
    stop("reference alignment failed: unmatched=", sum(unmatched),
         ", incompatible=", sum(incompatible), ", ambiguous=", sum(ambiguous),
         ", duplicate=", sum(duplicate_target),
         ", frequency_mismatch=", sum(status == "frequency_mismatch"), call. = FALSE)
  }

  if (any(flip_flag)) {
    data$beta[flip_flag] <- -data$beta[flip_flag]
    data$z[flip_flag] <- -data$z[flip_flag]
    data$effect_allele_frequency[flip_flag] <- 1 - data$effect_allele_frequency[flip_flag]
    if ("odds_ratio" %in% names(data)) {
      valid_or <- is.finite(data$odds_ratio[flip_flag]) & data$odds_ratio[flip_flag] > 0
      flipped_or <- data$odds_ratio[flip_flag]
      flipped_or[valid_or] <- 1 / flipped_or[valid_or]
      data$odds_ratio[flip_flag] <- flipped_or
    }
  }
  data$chromosome[aligned] <- ref$chromosome[key[aligned]]
  data$base_pair_location[aligned] <- ref$base_pair_location[key[aligned]]
  data$effect_allele[aligned] <- ref$alternate_allele[key[aligned]]
  data$other_allele[aligned] <- ref$reference_allele[key[aligned]]
  data$variant_id[aligned] <- ref$variant_id[key[aligned]]
  if ("rsid" %in% names(data) && "rsid" %in% names(ref)) {
    data$rsid[aligned] <- ref$rsid[key[aligned]]
  }
  data$variant_key <- ifelse(!is.na(key), ref$variant_key[key], NA_character_)
  data$harmonisation_status <- status
  data$harmonisation_flip <- flip_flag
  data$input_duplicate <- duplicate_target
  if (qc$frequency != "none") data$frequency_qc_status <- frequency_status
  data$.compressor_reference_index <- ifelse(aligned, ref$reference_index[key], NA_integer_)
  bad <- !aligned
  diagnostic_data <- data
  if (!isTRUE(preserve)) data <- data[aligned, , drop = FALSE]
  row.names(data) <- NULL
  counts <- list(
    input_rows = input_rows, output_rows = nrow(data), aligned_rows = sum(aligned),
    duplicate_input_rows = sum(duplicate_target), unmatched_rows = sum(unmatched),
    incompatible_rows = sum(incompatible), ambiguous_rows = sum(ambiguous),
    flipped_rows = sum(flip_flag),
    dropped_duplicates = if (isTRUE(preserve)) 0L else sum(duplicate_target),
    dropped_unmatched = if (isTRUE(preserve)) 0L else sum(unmatched),
    dropped_incompatible = if (isTRUE(preserve)) 0L else sum(incompatible),
    dropped_ambiguous = if (isTRUE(preserve)) 0L else sum(ambiguous),
    retained_unresolved_rows = if (isTRUE(preserve)) sum(bad) else 0L
  )
  if (qc$frequency != "none") {
    counts$frequency_checked_rows <- sum(frequency_status %in% c("match", "mismatch"))
    counts$frequency_mismatch_rows <- sum(frequency_status == "mismatch")
    counts$frequency_missing_rows <- sum(frequency_status %in% c("missing_study", "missing_reference"))
    counts$dropped_frequency_mismatch <- if (isTRUE(preserve) || qc$frequency == "report") 0L else
      sum(frequency_status == "mismatch")
  }
  counts$diagnostic_counts <- counts
  counts$diagnostic_examples <- alignment_example_rows(diagnostic_data, status, qc$example_limit,
                                                       extra_status = NULL)
  counts$counts <- counts$diagnostic_counts
  counts$examples <- counts$diagnostic_examples
  counts$diagnostics <- alignment_diagnostics(diagnostic_data, status,
                                              counts$diagnostic_counts,
                                              qc$example_limit, frequency_status)
  list(data = data, reference_hash = reference_hash(ref), reference_rows = nrow(ref),
       reference_metadata = attr(ref, "reference_metadata"),
       alignment_stats = counts)
}

is_partitioned_reference <- function(reference) {
  resolved <- resolve_reference(reference)
  (is.character(resolved$variants) && length(resolved$variants) == 1L &&
     dir.exists(resolved$variants)) || isTRUE(resolved$metadata$partitioned)
}

combine_alignment_statistics <- function(parts, input_rows, output_rows) {
  first <- parts[[1L]]$alignment_stats
  scalar_names <- names(first)[vapply(first, function(value) {
    is.numeric(value) && length(value) == 1L
  }, logical(1L))]
  stats <- lapply(scalar_names, function(name) {
    sum(vapply(parts, function(part) as.numeric(part$alignment_stats[[name]]), numeric(1)))
  })
  names(stats) <- scalar_names
  stats$input_rows <- as.integer(input_rows)
  stats$output_rows <- as.integer(output_rows)
  examples <- lapply(parts, function(part) part$alignment_stats$diagnostic_examples %||% list())
  example_names <- unique(unlist(lapply(examples, names), use.names = FALSE))
  merged_examples <- list()
  limit <- (.compressor_harmonise_context$qc %||% normalise_reference_qc())$example_limit
  for (name in example_names) {
    values <- unlist(lapply(examples, function(part) part[[name]] %||% list()),
                     recursive = FALSE, use.names = FALSE)
    if (length(values)) merged_examples[[name]] <- utils::head(values, limit)
  }
  stats$diagnostic_counts <- stats
  stats$diagnostic_examples <- merged_examples
  stats$counts <- stats$diagnostic_counts
  stats$examples <- stats$diagnostic_examples
  stats$diagnostics <- list(counts = stats$diagnostic_counts,
                            examples = merged_examples)
  stats
}

harmonise_by_chromosome <- function(data, reference, strict = FALSE, preserve = FALSE,
                                    chrom_threads = 1L) {
  resolved <- resolve_reference(reference)
  partitioned <- (is.character(resolved$variants) &&
                  length(resolved$variants) == 1L &&
                  dir.exists(resolved$variants)) ||
    isTRUE(resolved$metadata$partitioned)
  ref <- if (partitioned) NULL else reference_table(resolved)

  if (!partitioned && (is.null(ref) || nrow(data) < 2L ||
                       as.integer(chrom_threads) <= 1L)) {
    return(harmonise_to_reference(data, resolved, strict = strict, preserve = preserve,
                                  reference_data = ref))
  }
  if (!nrow(data)) {
    return(harmonise_to_reference(data, resolved, strict = strict, preserve = preserve,
                                  reference_data = if (partitioned) empty_reference_table() else ref))
  }

  chrom_threads <- max(1L, as.integer(chrom_threads))
  order_column <- ".compressor_input_order"
  if (order_column %in% names(data)) {
    stop("reserved column name is present: ", order_column, call. = FALSE)
  }
  data[[order_column]] <- seq_len(nrow(data))
  group_label <- as.character(data$chromosome)
  group_label[is.na(group_label) | !nzchar(group_label)] <- "__UNRESOLVED__"
  groups <- split(seq_len(nrow(data)), group_label, drop = TRUE)
  alignment_qc <- .compressor_harmonise_context$qc %||% normalise_reference_qc()
  worker <- function(idx) {
    had_qc <- exists("qc", envir = .compressor_harmonise_context, inherits = FALSE)
    old_qc <- if (had_qc) .compressor_harmonise_context$qc else NULL
    .compressor_harmonise_context$qc <- alignment_qc
    on.exit({
      if (had_qc) .compressor_harmonise_context$qc <- old_qc else
        rm("qc", envir = .compressor_harmonise_context)
    }, add = TRUE)
    ref_part <- if (partitioned) {
      chromosome <- as.character(data$chromosome[idx[1L]])
      if (is.na(chromosome) || !nzchar(chromosome)) empty_reference_table() else
        reference_table(resolved, chromosome = chromosome)
    } else {
      ref
    }
    harmonise_to_reference(
      data[idx, , drop = FALSE], NULL, strict = FALSE, preserve = preserve,
      reference_data = ref_part,
      reference_key = ref_part$variant_key
    )
  }
  if (chrom_threads <= 1L || length(groups) <= 1L) {
    parts <- lapply(groups, worker)
  } else if (.Platform$OS.type == "windows") {
    cluster <- parallel::makeCluster(min(chrom_threads, length(groups)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parts <- parallel::parLapply(cluster, groups, worker)
  } else {
    parts <- parallel::mclapply(groups, worker, mc.cores = min(chrom_threads, length(groups)),
                                mc.preschedule = FALSE)
  }
  failed <- vapply(parts, inherits, logical(1L), what = "try-error")
  if (any(failed)) {
    stop("chromosome-parallel harmonisation worker failed: ",
         as.character(parts[[which(failed)[1L]]]), call. = FALSE)
  }
  combined <- do.call(rbind, lapply(parts, function(x) x$data))
  combined <- combined[order(combined[[order_column]]), , drop = FALSE]
  combined[[order_column]] <- NULL
  row.names(combined) <- NULL
  stats <- combine_alignment_statistics(parts, nrow(data), nrow(combined))
  if (isTRUE(strict) && any(unlist(stats[c(
    "unmatched_rows", "incompatible_rows", "ambiguous_rows",
    "duplicate_input_rows"
  )], use.names = FALSE) > 0)) {
    stop("reference alignment failed: unmatched=", stats$unmatched_rows,
         ", incompatible=", stats$incompatible_rows,
         ", ambiguous=", stats$ambiguous_rows,
         ", duplicate=", stats$duplicate_input_rows, call. = FALSE)
  }
  list(data = combined, reference_hash = parts[[1L]]$reference_hash,
       reference_rows = if (partitioned) resolved$metadata$rows else parts[[1L]]$reference_rows,
       reference_metadata = if (partitioned) resolved$metadata else parts[[1L]]$reference_metadata,
       alignment_stats = stats)
}

# Keep alias lookup on the same bounded query path as coordinate lookup. The
# compact EBI master is row-group framed, so requested rsIDs are decoded one
# group at a time and only matching rows are retained; a whole chromosome is
# no longer materialised as an intermediate alias table.
reference_alias_table <- function(reference, ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids) & ids != "."]
  empty <- data.frame(rsid = character(), variant_id = character(),
                      variant_index = integer(), chromosome = character(),
                      base_pair_location = integer(), stringsAsFactors = FALSE)
  if (!length(ids)) return(empty)
  resolved <- resolve_reference(reference)
  path <- resolved$variants
  if (is.character(path) && length(path) == 1L && dir.exists(path) &&
      file.exists(file.path(path, "variants.parquet"))) {
    ref <- reference_table(resolved, rsids = ids)
    if (!nrow(ref)) return(empty)
    aliases <- reference_alias_rows(ref)
    return(aliases[aliases$rsid %in% ids, , drop = FALSE])
  }
  if (is.character(path) && length(path) == 1L && dir.exists(path)) {
    alias_dir <- resolved$metadata$aliases_dir %||% "aliases"
    if (!grepl("^[/~]", alias_dir)) alias_dir <- file.path(path, alias_dir)
    alias_files <- list.files(alias_dir, pattern = "[.]parquet$",
                              full.names = TRUE, ignore.case = TRUE)
    if (!length(alias_files)) return(empty)
    parts <- lapply(alias_files, function(file) {
      part <- as.data.frame(arrow::read_parquet(file))
      part[part$rsid %in% ids, , drop = FALSE]
    })
    parts <- Filter(function(x) nrow(x), parts)
    if (!length(parts)) return(empty)
    out <- do.call(rbind, parts)
    if (!"chromosome" %in% names(out)) {
      parsed <- parse_canonical_variant_keys(out$variant_id)
      out$chromosome <- parsed$chromosome
      out$base_pair_location <- parsed$base_pair_location
    }
    return(out[!duplicated(paste(out$rsid, out$variant_id, sep = "\r")), , drop = FALSE])
  }
  ref <- reference_table(resolved, rsids = ids)
  if (!nrow(ref)) return(empty)
  aliases <- reference_alias_rows(ref)
  aliases[aliases$rsid %in% ids, , drop = FALSE]
}
