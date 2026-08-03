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

  if (!"variant_id" %in% names(data)) data$variant_id <- NA_character_
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

  if (!"rsid" %in% names(data)) data$rsid <- NA_character_
  data$rsid <- as.character(data$rsid)
  data$rsid[data$rsid %in% c("", ".", "NA")] <- NA_character_
  if (!"effect_allele_frequency" %in% names(data)) data$effect_allele_frequency <- NA_real_
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
  data$variant_id <- data$variant_key
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
    paste("gzip -dc", shQuote(source), "| grep -v '^##'")
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
      args$cmd <- paste("gzip -dc", shQuote(normalizePath(source)))
    } else {
      args$input <- source
    }
    return(do.call(data.table::fread, args))
  }
  if (grepl("[.]vcf(?:[.]gz|[.]bgz)?$", source, ignore.case = TRUE)) {
    return(read_reference_vcf_source(source))
  }
  if (grepl("[.]pvar(?:[.]gz|[.]zst)?$", source, ignore.case = TRUE)) {
    if (grepl("[.]zst$", source, ignore.case = TRUE)) {
      return(data.table::fread(cmd = paste("zstd -dc --quiet", shQuote(normalizePath(source))),
                               data.table = FALSE, showProgress = FALSE, check.names = FALSE))
    }
    if (grepl("[.]gz$", source, ignore.case = TRUE)) {
      return(data.table::fread(cmd = paste("gzip -dc", shQuote(normalizePath(source))),
                               data.table = FALSE, showProgress = FALSE, check.names = FALSE))
    }
    return(data.table::fread(source, data.table = FALSE, showProgress = FALSE, check.names = FALSE))
  }
  if (grepl("[.]gz$", source, ignore.case = TRUE)) {
    return(data.table::fread(cmd = paste("gzip -dc", shQuote(normalizePath(source))),
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
    mapped_rows <- source_rows[one]
    out$chromosome[mapped_rows] <- sub(
      "^chr", "", as.character(GenomicRanges::seqnames(mapped)), ignore.case = TRUE
    )
    out$base_pair_location[mapped_rows] <- GenomicRanges::start(mapped)
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
  data <- canonicalize_reference_data(read_reference_source(source), source_build, chain,
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

reference_table <- function(reference, chromosome = NULL) {
  if (is.null(reference)) return(NULL)
  resolved <- resolve_reference(reference)
  if (is.data.frame(resolved$variants)) {
    out <- if (!nrow(resolved$variants)) empty_reference_table() else
      normalise_reference_columns(resolved$variants)
  } else {
    path <- resolved$variants
    if (dir.exists(path)) {
      master <- file.path(path, "variants.parquet")
      if (file.exists(master)) {
        if (is.null(chromosome)) {
          raw <- arrow::read_parquet(master)
        } else {
          code <- chromosome_sort_key(chromosome)[1L]
          reader <- arrow::ParquetFileReader$create(master)
          groups <- seq_len(reader$num_row_groups) - 1L
          selected <- vapply(groups, function(i) {
            header <- as.data.frame(reader$ReadRowGroup(i, column_indices = 0L))
            isTRUE(nrow(header) > 0L && header$chromosome_code[1L] == code)
          }, logical(1L))
          selected <- groups[selected]
          raw <- if (length(selected)) {
            as.data.frame(reader$ReadRowGroups(selected))
          } else {
            data.frame()
          }
        }
        out <- if (!nrow(raw)) empty_reference_table() else
          normalise_reference_columns(raw)
      } else {
        paths <- reference_partition_paths(path, chromosome)
        out <- if (!length(paths)) {
          empty_reference_table()
        } else {
          normalise_reference_columns(do.call(rbind, lapply(paths, arrow::read_parquet)))
        }
      }
    } else if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
      stop("resolved reference does not contain a readable variants asset", call. = FALSE)
    } else if (grepl("[.]parquet$", path, ignore.case = TRUE)) {
      out <- normalise_reference_columns(arrow::read_parquet(path))
    } else if (grepl("[.]bim(?:[.]gz)?$", path, ignore.case = TRUE)) {
      out <- normalise_reference_columns(read_reference_source(path))
    } else {
      normalized_cache <- reference_normalized_cache_path(resolved)
      if (!is.null(normalized_cache) && file.exists(normalized_cache)) {
        out <- normalise_reference_columns(arrow::read_parquet(normalized_cache))
      } else {
        out <- normalise_reference_columns(read_reference_source(path))
        if (!is.null(normalized_cache)) write_normalized_reference_cache(out, normalized_cache)
      }
    }
  }
  if (!is.null(chromosome)) {
    chromosome <- sub("^chr", "", as.character(chromosome), ignore.case = TRUE)
    out <- out[out$chromosome == chromosome, , drop = FALSE]
    row.names(out) <- NULL
  }
  if (anyNA(out$variant_key) || anyDuplicated(out$variant_key)) {
    stop("canonical reference has missing or duplicated chrom:pos:ref:alt keys", call. = FALSE)
  }
  if (!"variant_index" %in% names(out)) {
    offset <- if (is.null(chromosome)) 0L else
      reference_partition_offset(resolved$metadata %||% list(), chromosome)
    out$variant_index <- as.integer(offset + seq_len(nrow(out)) - 1L)
  }
  out$variant_index <- as.integer(out$variant_index)
  if (anyDuplicated(out$variant_index) || any(out$variant_index < 0L, na.rm = TRUE)) {
    stop("canonical reference has invalid or duplicated variant_index values", call. = FALSE)
  }
  out$reference_index <- out$variant_index
  metadata <- resolved$metadata %||% list(
    id = resolved$id %||% "local",
    build = resolved$build %||% "GRCh38"
  )
  if (!is.null(chromosome)) metadata$partition_chromosome <- as.character(chromosome)
  normalized_cache <- reference_normalized_cache_path(resolved)
  if (!is.null(normalized_cache) && file.exists(normalized_cache)) {
    metadata$normalized_cache_path <- normalizePath(normalized_cache, mustWork = FALSE)
  }
  attr(out, "reference_metadata") <- metadata
  if (!is.null(metadata$sha256)) attr(out, "reference_sha256") <- metadata$sha256
  out
}

reference_hash <- function(reference_data) {
  if (is.null(reference_data)) return(NA_character_)
  known_hash <- attr(reference_data, "reference_sha256")
  if (!is.null(known_hash) && length(known_hash) == 1L && nzchar(known_hash)) return(known_hash)
  digest::digest(reference_data, algo = "sha256", serialize = TRUE)
}

harmonise_to_reference <- function(data, reference, strict = FALSE, preserve = FALSE,
                                   reference_data = NULL, reference_key = NULL) {
  input_rows <- nrow(data)
  ref <- if (!is.null(reference_data)) reference_data else reference_table(reference)
  if (is.null(ref)) {
    data$harmonisation_status <- "unreferenced"
    data$harmonisation_flip <- FALSE
    data$input_duplicate <- duplicated(data$variant_id)
    return(list(data = data, reference_hash = NA_character_,
                reference_rows = NA_integer_, reference_metadata = NULL,
                alignment_stats = list(
                  input_rows = input_rows, output_rows = nrow(data), aligned_rows = 0L,
                  duplicate_input_rows = sum(data$input_duplicate), unmatched_rows = 0L,
                  incompatible_rows = 0L, ambiguous_rows = 0L, dropped_duplicates = 0L,
                  dropped_unmatched = 0L, dropped_incompatible = 0L, dropped_ambiguous = 0L
                )))
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
  complement_i <- match(complement_key, ref_key)
  complement_reverse_i <- match(complement_reverse_key, ref_key)
  candidates <- cbind(!is.na(direct_i), !is.na(reverse_i),
                      !is.na(complement_i), !is.na(complement_reverse_i))
  candidate_count <- rowSums(candidates)
  unmatched <- candidate_count == 0L
  multi <- candidate_count > 1L
  key <- direct_i
  key[is.na(key)] <- reverse_i[is.na(key)]
  key[is.na(key)] <- complement_i[is.na(key)]
  key[is.na(key)] <- complement_reverse_i[is.na(key)]
  ambiguous <- multi
  incompatible <- !unmatched & !multi & !is.na(key) &
    !(is_reference_allele(data$effect_allele) & is_reference_allele(data$other_allele))
  # Palindromic strand matches are not resolvable from alleles alone. If both
  # the study and reference carry EAF, accept only a clearly separated match.
  ref_eaf <- ref$effect_allele_frequency[key]
  same_error <- abs(data$effect_allele_frequency - ref_eaf)
  flip_error <- abs((1 - data$effect_allele_frequency) - ref_eaf)
  informative <- multi & is.finite(same_error) & is.finite(flip_error)
  resolved_multi <- informative & abs(same_error - flip_error) > 0.1
  ambiguous[resolved_multi] <- FALSE
  status <- rep("aligned", input_rows)
  status[unmatched] <- "unmatched"
  status[incompatible] <- "incompatible"
  status[ambiguous] <- "ambiguous"
  flip_flag <- (!is.na(key) & !is.na(reverse_i) & reverse_i == key) |
    (!is.na(key) & !is.na(complement_reverse_i) & complement_reverse_i == key)
  flip_flag[is.na(flip_flag)] <- FALSE
  if (any(resolved_multi)) flip_flag[resolved_multi] <- flip_error[resolved_multi] < same_error[resolved_multi]
  flip_flag[unmatched | incompatible | ambiguous] <- FALSE
  aligned <- !unmatched & !incompatible & !ambiguous

  # Any second input row targeting the same canonical variant is discarded as
  # ambiguous. Keeping neither row avoids silently choosing a study duplicate.
  duplicate_target <- rep(FALSE, input_rows)
  duplicate_target[aligned] <- duplicated(key[aligned]) | duplicated(key[aligned], fromLast = TRUE)
  status[duplicate_target] <- "duplicate"
  aligned[duplicate_target] <- FALSE
  if (isTRUE(strict) && any(unmatched | incompatible | ambiguous | duplicate_target)) {
    stop("reference alignment failed: unmatched=", sum(unmatched),
         ", incompatible=", sum(incompatible), ", ambiguous=", sum(ambiguous),
         ", duplicate=", sum(duplicate_target), call. = FALSE)
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
  data$.compressor_reference_index <- ifelse(aligned, ref$reference_index[key], NA_integer_)
  bad <- !aligned
  if (!isTRUE(preserve)) data <- data[aligned, , drop = FALSE]
  row.names(data) <- NULL
  list(data = data, reference_hash = reference_hash(ref), reference_rows = nrow(ref),
       reference_metadata = attr(ref, "reference_metadata"),
       alignment_stats = list(
         input_rows = input_rows, output_rows = nrow(data), aligned_rows = sum(aligned),
         duplicate_input_rows = sum(duplicate_target), unmatched_rows = sum(unmatched),
         incompatible_rows = sum(incompatible), ambiguous_rows = sum(ambiguous),
         flipped_rows = sum(flip_flag),
         dropped_duplicates = if (isTRUE(preserve)) 0L else sum(duplicate_target),
         dropped_unmatched = if (isTRUE(preserve)) 0L else sum(unmatched),
         dropped_incompatible = if (isTRUE(preserve)) 0L else sum(incompatible),
         dropped_ambiguous = if (isTRUE(preserve)) 0L else sum(ambiguous),
         retained_unresolved_rows = if (isTRUE(preserve)) sum(bad) else 0L
       ))
}

is_partitioned_reference <- function(reference) {
  resolved <- resolve_reference(reference)
  (is.character(resolved$variants) && length(resolved$variants) == 1L &&
     dir.exists(resolved$variants)) || isTRUE(resolved$metadata$partitioned)
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
  worker <- function(idx) {
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
  stats_names <- names(parts[[1L]]$alignment_stats)
  stats <- lapply(stats_names, function(nm) {
    sum(vapply(parts, function(x) as.numeric(x$alignment_stats[[nm]]), numeric(1)))
  })
  names(stats) <- stats_names
  stats$input_rows <- nrow(data)
  stats$output_rows <- nrow(combined)
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
