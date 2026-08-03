#' Write a CompreSSoR store
#'
#' @param input A data.frame or a delimited summary-statistics file.
#' @param output Destination directory for the store.
#' @param reference GRCh38 reference spine, `"GRCh38"` for the immutable
#'   reference configured by `COMPRESSOR_CANONICAL_REFERENCE`, or `NULL` to
#'   skip reference alignment.
#' @param mode One of `"qc"` (default, harmonise and drop unresolved rows),
#'   `"convert"` (minimal conversion without reference QC), `"all"` (alias
#'   for `"qc"`), `"core"`, or `"hm3"`.
#' @param variant_set Panel data.frame or file used by `mode = "core"` or
#'   `mode = "hm3"`. PLINK `.bim`, Parquet and delimited panel files are
#'   supported.
#' @param strict If `TRUE`, fail when variants are absent, incompatible,
#'   ambiguous, duplicated, or unsupported by the selected backend.
#' @param drop_unresolved Whether unmatched, incompatible, ambiguous and
#'   duplicate rows are dropped. The default is `TRUE`.
#' @param input_build Build of the input. Non-GRCh38 input requires `chain`.
#' @param chain Optional GRCh37-to-GRCh38 chain file.
#' @param chrom_threads Number of chromosome workers for harmonisation. Values
#'   above one use chromosome-parallel harmonisation and share the reference.
#' @param profile `"standard"` uses semantic Z9/EAF8/SE6 streams with sparse
#'   float32 exceptions; `"exact"` is available with the Parquet backend. P and
#'   beta are derived rather than stored.
#' @param backend Storage backend. The default, `"pcodec"`, is the compact
#'   block-framed Pcodec format with a self-contained GRCh38 identity key.
#'   `"parquet"` retains the interoperable legacy backend.
#' @param assume_grch38_ref_alt Required for Pcodec `mode = "convert"` unless
#'   the input has explicit REF and ALT columns. Set `TRUE` only when positions
#'   are GRCh38, `other_allele` is REF, `effect_allele` is ALT, and beta/Z/EAF
#'   refer to ALT. Normal QC mode establishes this contract from the reference.
#' @param keep_extras Whether arbitrary non-core input columns should be kept
#'   in an extras sidecar. The default is `FALSE`; harmonisation/QC summaries
#'   remain in the manifest, while row-wise QC labels are retained only when
#'   this is explicitly set to `TRUE`.
#' @param overwrite Whether an existing destination may be replaced.
#' @param cache Whether to build the optional q8 framed cache after writing.
#' @param block_rows Number of rows per Parquet row group; Pcodec uses its
#'   measured fixed key and value frame sizes.
#' @return A `compressor_store` object.
#' @examples
#' \dontrun{
#' example <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
#' store <- compress_sumstats(
#'   example, file.path(tempdir(), "example.cpr"), mode = "convert",
#'   reference = NULL, assume_grch38_ref_alt = TRUE, overwrite = TRUE
#' )
#' read_sumstats(store, columns = c("chromosome", "base_pair_location", "z"))
#' }
#' @export
compress_sumstats <- function(input, output, reference = "GRCh38",
                              profile = c("standard", "exact"),
                              overwrite = FALSE, cache = FALSE, strict = FALSE,
                              keep_extras = FALSE,
                              block_rows = 65536L,
                              mode = c("qc", "convert", "all", "core", "hm3"),
                              variant_set = NULL, chrom_threads = 1L,
                              drop_unresolved = TRUE,
                              input_build = "GRCh38", chain = NULL,
                              backend = c("pcodec", "parquet"),
                              assume_grch38_ref_alt = FALSE) {
  profile <- match.arg(profile)
  backend <- match.arg(backend)
  mode <- match.arg(mode)
  if (identical(backend, "parquet")) {
    require_parquet_backend("backend='parquet'", dplyr = TRUE)
  }
  if (length(assume_grch38_ref_alt) != 1L || !is.logical(assume_grch38_ref_alt) ||
      is.na(assume_grch38_ref_alt)) {
    stop("assume_grch38_ref_alt must be TRUE or FALSE", call. = FALSE)
  }
  if (identical(backend, "pcodec") && !identical(mode, "convert") &&
      !isTRUE(drop_unresolved)) {
    stop("Pcodec canonical stores cannot retain unresolved rows; use the Parquet backend for an audit store",
         call. = FALSE)
  }
  if (length(keep_extras) != 1L || !is.logical(keep_extras) || is.na(keep_extras)) {
    stop("keep_extras must be TRUE or FALSE", call. = FALSE)
  }
  transaction <- stage_store_output(output, overwrite = overwrite)
  completed <- FALSE
  on.exit(if (!completed && dir.exists(transaction$staging)) {
    unlink(transaction$staging, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  output <- transaction$staging

  raw <- import_sumstats(input)
  source_keys <- alias_key(attr(raw, "source_columns") %||% character())
  explicit_ref_alt <- isTRUE(attr(raw, "explicit_ref_alt")) ||
    (any(source_keys %in% c("ref", "referenceallele")) &&
    any(source_keys %in% c("alt", "alternateallele"))
    )
  if (identical(backend, "pcodec") && identical(mode, "convert") &&
      !isTRUE(assume_grch38_ref_alt) && !explicit_ref_alt) {
    stop("Pcodec mode='convert' needs explicit REF/ALT columns or assume_grch38_ref_alt=TRUE; ordinary effect/other alleles do not prove REF/ALT orientation",
         call. = FALSE)
  }
  prepared <- prepare_sumstats_data(raw, reference, mode = mode, variant_set = variant_set,
                                    strict = strict, chrom_threads = chrom_threads,
                                    drop_unresolved = drop_unresolved,
                                    input_build = input_build, chain = chain)
  alignment <- prepared$alignment
  data <- prepared$data
  reference_index <- if (".compressor_reference_index" %in% names(data)) {
    as.integer(data$.compressor_reference_index)
  } else {
    NULL
  }
  data$.compressor_reference_index <- NULL
  validate_sumstats_values(data, require_identity = identical(backend, "pcodec") || isTRUE(drop_unresolved))
  if (identical(backend, "pcodec")) {
    if (!identical(profile, "standard")) {
      stop("backend='pcodec' currently provides the standard semantic profile; "
           , "use backend='parquet' for profile='exact'", call. = FALSE)
    }
    if (isTRUE(cache)) {
      stop("backend='pcodec' is already block-framed for regional access; "
           , "cache=TRUE is not needed", call. = FALSE)
    }
    if (isTRUE(keep_extras)) {
      stop("backend='pcodec' currently stores only the core summary-statistics "
           , "columns; use backend='parquet' to retain extras", call. = FALSE)
    }
    supported_chromosome <- as.character(data$chromosome) %in%
      c(as.character(1:22), "X", "Y")
    bad_allele <- !supported_chromosome |
      is.na(data$effect_allele) | is.na(data$other_allele) |
      is.na(nchar(data$effect_allele)) | is.na(nchar(data$other_allele)) |
      nchar(data$effect_allele) != 1L | nchar(data$other_allele) != 1L |
      !data$effect_allele %in% c("A", "C", "G", "T") |
      !data$other_allele %in% c("A", "C", "G", "T") |
      data$effect_allele == data$other_allele
    unsupported_rows <- sum(bad_allele)
    if (unsupported_rows && isTRUE(strict)) {
      stop("Pcodec identity supports biallelic A/C/G/T SNVs on chromosomes 1-22, X and Y; ",
           unsupported_rows, " unsupported row(s) were found", call. = FALSE)
    }
    if (unsupported_rows) {
      data <- data[!bad_allele, , drop = FALSE]
      stats <- alignment$alignment_stats %||% list()
      stats$unsupported_pcodec_rows <- as.integer(unsupported_rows)
      stats$dropped_unsupported_pcodec <- as.integer(unsupported_rows)
      alignment$alignment_stats <- stats
    }
    if (!nrow(data)) {
      stop("no supported biallelic GRCh38 SNVs remain for the Pcodec store", call. = FALSE)
    }
    input_build_key <- toupper(gsub("[. -]", "", as.character(input_build)))
    input_is_grch38 <- input_build_key %in% c("GRCH38", "HG38", "38")
    if (!input_is_grch38 && is.null(chain)) {
      stop("backend='pcodec' requires a GRCh37-to-GRCh38 chain for non-GRCh38 input", call. = FALSE)
    }
    reference_metadata <- alignment$reference_metadata %||% list(
      id = "none", build = "GRCh38", status = "reference alignment skipped"
    )
    reference_metadata$sha256 <- alignment$reference_hash
    reference_metadata$rows <- alignment$reference_rows
    pcodec_metadata <- list(
      genome_build = "GRCh38",
      reference = reference_metadata,
      harmonisation = list(
        method = if (identical(prepared$effective_mode, "convert") || is.null(reference)) "none" else "reference_spine",
        mode = prepared$requested_mode,
        chrom_threads = as.integer(chrom_threads),
        strict = isTRUE(strict),
        drop_unresolved = isTRUE(drop_unresolved),
        input_build = as.character(input_build),
        liftover_chain = chain_manifest_metadata(chain),
        alignment = alignment$alignment_stats
      ),
      source_columns = attr(raw, "source_columns") %||% names(raw),
      source = attr(raw, "source_provenance") %||% NULL
    )
    pcodec_write_store(data, output, metadata = pcodec_metadata)
    commit_store_output(transaction)
    completed <- TRUE
    return(open_compressor(transaction$target))
  }

  n <- nrow(data)

  variant_names <- intersect(c("row", "chromosome", "base_pair_location", "effect_allele",
                               "other_allele", "variant_id", "rsid"), names(data))
  if (all(c("variant_id", "rsid") %in% variant_names)) {
    same_identity <- (is.na(data$variant_id) & is.na(data$rsid)) |
      (!is.na(data$variant_id) & !is.na(data$rsid) &
         as.character(data$variant_id) == as.character(data$rsid))
    if (all(same_identity)) variant_names <- setdiff(variant_names, "rsid")
  }
  variants <- data[intersect(variant_names, names(data))]
  variants$row <- seq_len(n) - 1L
  variants <- variants[c("row", setdiff(names(variants), "row"))]
  reference_metadata <- alignment$reference_metadata %||% list()
  reusable_reference <- !identical(reference_metadata$id %||% "", "in_memory") &&
    any(vapply(c("local_path", "normalized_cache_path", "source_url"),
               function(field) !is.null(reference_metadata[[field]]) &&
                 nzchar(as.character(reference_metadata[[field]])), logical(1)))
  reference_indexed <- !is.null(reference_index) && reusable_reference
  files <- list(variants = "variants.parquet", values = "values.parquet",
                exceptions = NULL, unmatched = NULL, extras = NULL)
  if (reference_indexed) {
    index_table <- data.frame(
      row = seq_len(n) - 1L,
      reference_index = reference_index,
      stringsAsFactors = FALSE
    )
    arrow::write_parquet(index_table, file.path(output, "variants.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, use_dictionary = TRUE,
                         chunk_size = as.integer(block_rows))
    unmatched <- variants[is.na(reference_index), , drop = FALSE]
    if (nrow(unmatched)) {
      arrow::write_parquet(unmatched, file.path(output, "unmatched.parquet"),
                           compression = "zstd", compression_level = 7,
                           write_statistics = TRUE, use_dictionary = TRUE,
                           chunk_size = as.integer(block_rows))
      files$unmatched <- "unmatched.parquet"
    }
  } else {
    arrow::write_parquet(variants, file.path(output, "variants.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, use_dictionary = TRUE,
                         chunk_size = as.integer(block_rows))
  }

  logical_columns <- c("z", "standard_error", "effect_allele_frequency")
  if (profile == "exact") {
    values <- data[logical_columns]
    values$row <- seq_len(n) - 1L
    values <- values[c("row", logical_columns)]
    arrow::write_parquet(values, file.path(output, "values.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, chunk_size = as.integer(block_rows))
    codec <- list(name = "lossless_semantic_numeric", z_bits = NULL, se_bits = NULL,
                  eaf_bits = NULL, p_storage = "omitted; derived from z")
    files$exceptions <- NULL
  } else {
    encoded <- q_encode(data$beta, data$standard_error, data$effect_allele_frequency,
                        z = data$z, z_bits = 9L, se_bits = 6L, eaf_bits = 8L,
                        block_rows = block_rows)
    arrow::write_parquet(encoded$main, file.path(output, "values.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, chunk_size = as.integer(block_rows))
    arrow::write_parquet(encoded$exceptions, file.path(output, "exceptions.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, chunk_size = as.integer(block_rows))
    codec <- encoded$metadata
    codec$compression <- "zstd"
    files$exceptions <- "exceptions.parquet"
  }

  # Harmonisation/QC is performed before writing and its aggregate counts are
  # retained in the manifest. Row-wise status strings are deliberately not
  # part of the standard durable payload; callers can retain them explicitly
  # with keep_extras = TRUE.
  core_names <- unique(c(variant_names, logical_columns, "beta", "p_value", "odds_ratio", "rsid"))
  extra_names <- if (isTRUE(keep_extras)) setdiff(names(data), core_names) else character()
  extra_names <- extra_names[vapply(data[extra_names], function(x) is.atomic(x) && !is.list(x), logical(1))]
  if (length(extra_names)) {
    dir.create(file.path(output, "extras"), showWarnings = FALSE)
    extras <- data[c(extra_names)]
    extras$row <- seq_len(n) - 1L
    extras <- extras[c("row", extra_names)]
    arrow::write_parquet(extras, file.path(output, "extras", "values.parquet"),
                         compression = "zstd", compression_level = 7,
                         chunk_size = as.integer(block_rows))
    files$extras <- file.path("extras", "values.parquet")
  }

  reference_manifest <- alignment$reference_metadata %||% list(
    id = "none",
    build = "GRCh38",
    status = "reference alignment skipped"
  )
  reference_manifest$sha256 <- alignment$reference_hash
  reference_manifest$rows <- alignment$reference_rows
  manifest <- list(
    format = "CompreSSoR",
    format_version = "0.1.0",
    profile = profile,
    backend = "parquet",
    genome_build = prepared$genome_build,
    reference = reference_manifest,
    harmonisation = list(
      method = if (identical(prepared$effective_mode, "convert") || is.null(reference)) "none" else "reference_spine",
      mode = prepared$requested_mode,
      chrom_threads = as.integer(chrom_threads),
      strict = isTRUE(strict),
      drop_unresolved = isTRUE(drop_unresolved),
      input_build = as.character(input_build),
      liftover_chain = chain_manifest_metadata(chain),
      alignment = alignment$alignment_stats
    ),
    n_rows = n,
    block_rows = as.integer(block_rows),
    logical_columns = logical_columns,
    derived_columns = list(beta = "z * standard_error", p_value = "2 * pnorm(-abs(z))"),
    source_columns = attr(raw, "source_columns") %||% names(raw),
    source = attr(raw, "source_provenance") %||% NULL,
    codec = codec,
    variant_storage = if (reference_indexed) "shared_reference_index" else "inline",
    reference_index_base = if (reference_indexed) 0L else NULL,
    variant_identity = list(
      rsid = if ("rsid" %in% variant_names) "stored" else "derived_from_variant_id"
    ),
    files = files,
    benchmark = benchmark_metadata(),
    benchmark_comparisons = list(
      vcf_tabix = vcf_benchmark_metadata(),
      finngen_end_to_end = finngen_benchmark_metadata(),
      finngen_optimization = finngen_optimization_metadata(),
      modes_edge = mode_benchmark_metadata(),
      release_gate = release_gate_benchmark_metadata(),
      release_gate_followup = release_gate_followup_metadata(),
      storage_size = storage_size_benchmark_metadata(),
      storage_amortization = storage_amortization_benchmark_metadata()
    ),
    created_utc = now_utc(),
    tolerances = if (profile == "standard") list(eaf_abs_max = 0.0002, z_central_range = c(-3.5, 3.5)) else list(exact = TRUE)
  )
  write_manifest(manifest, file.path(output, "manifest.json"))
  store <- open_compressor(output)
  if (isTRUE(cache)) build_cache(store, overwrite = TRUE, block_rows = block_rows)
  commit_store_output(transaction)
  completed <- TRUE
  open_compressor(transaction$target)
}

#' Open a CompreSSoR store
#'
#' @param path Store directory.
#' @return A `compressor_store` object.
#' @export
open_compressor <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) stop("store directory does not exist: ", path, call. = FALSE)
  manifest <- read_manifest(file.path(path, "manifest.json"))
  if (!identical(manifest$format, "CompreSSoR")) stop("not a CompreSSoR store", call. = FALSE)
  if (identical(manifest$backend, "pcodec")) {
    verify_pcodec_manifest(file.path(path, "manifest.json"))
  }
  if (identical(manifest$backend, "pcodec") &&
      !manifest$format_version %in% c("0.2.0-pcodec", "0.3.0-pcodec")) {
    stop("unsupported Pcodec format version: ", manifest$format_version %||% "missing", call. = FALSE)
  }
  structure(list(path = path, manifest = manifest), class = "compressor_store")
}

print.compressor_store <- function(x, ...) {
  m <- x$manifest
  cat("CompreSSoR store\n")
  cat("Path:        ", x$path, "\n", sep = "")
  cat("Rows:        ", m$n_rows, "\n", sep = "")
  cat("Genome:      ", m$genome_build, "\n", sep = "")
  cat("Profile:     ", m$profile, "\n", sep = "")
  cat("Backend:     ", m$backend, "\n", sep = "")
  cat("Codec:       ", m$codec$name, "\n", sep = "")
  cat("Block rows:  ", m$block_rows, "\n", sep = "")
  invisible(x)
}

read_region_bounds <- function(region) {
  if (is.null(region)) return(NULL)
  if (length(region) == 1L && is.character(region)) {
    parts <- regexec("^chr?([^:]+):([0-9]+)-([0-9]+)$", region, ignore.case = TRUE)
    hit <- regmatches(region, parts)[[1L]]
    if (length(hit) != 4L) stop("region must look like chr1:100-200", call. = FALSE)
    return(list(chromosome = sub("^chr", "", hit[2L], ignore.case = TRUE),
                start = as.numeric(hit[3L]), end = as.numeric(hit[4L])))
  }
  if (is.numeric(region) && length(region) == 3L) return(list(chromosome = as.character(region[1L]), start = region[2L], end = region[3L]))
  stop("region must be a string such as chr1:100-200 or c(chr, start, end)", call. = FALSE)
}

read_parquet_region <- function(path, bounds) {
  dataset <- arrow::open_dataset(path, format = "parquet")
  chromosome_value <- bounds$chromosome
  start_value <- bounds$start
  end_value <- bounds$end
  query <- dplyr::filter(
    dataset,
    chromosome == !!chromosome_value,
    base_pair_location >= !!start_value,
    base_pair_location <= !!end_value
  )
  dplyr::collect(query)
}

read_parquet_rows <- function(path, rows) {
  rows <- as.integer(rows)
  if (!length(rows)) return(NULL)
  dataset <- arrow::open_dataset(path, format = "parquet")
  query <- dplyr::filter(dataset, row %in% !!rows)
  dplyr::collect(query)
}

store_reference_descriptor <- function(store) {
  metadata <- store$manifest$reference %||% list()
  if (identical(metadata$id %||% "none", "none")) return(NULL)
  normalized <- metadata$normalized_cache_path %||% NULL
  local <- metadata$local_path %||% NULL
  variants <- if (!is.null(normalized) && file.exists(normalized)) {
    normalized
  } else if (!is.null(local) && file.exists(local)) {
    local
  } else if (!is.null(metadata$source_url) && nzchar(metadata$source_url)) {
    list(url = metadata$source_url, filename = metadata$filename,
         md5 = metadata$md5, sha256 = metadata$sha256)
  } else {
    stop("shared reference is not available; restore the canonical GRCh38 reference or pass its recorded path",
         call. = FALSE)
  }
  list(id = metadata$id %||% "canonical",
       build = metadata$build %||% "GRCh38",
       source = metadata$source %||% NULL,
       source_url = metadata$source_url %||% NULL,
       variants = variants,
       cache_dir = reference_cache_dir())
}

store_reference_table <- function(store) {
  reference <- store_reference_descriptor(store)
  if (is.null(reference)) return(NULL)
  out <- reference_table(reference)
  required <- c("chromosome", "base_pair_location", "reference_allele",
                "alternate_allele", "variant_id")
  missing <- setdiff(required, names(out))
  if (length(missing)) stop("canonical reference is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!"rsid" %in% names(out)) out$rsid <- NA_character_
  out$effect_allele <- out$alternate_allele
  out$other_allele <- out$reference_allele
  out$reference_index <- seq_len(nrow(out)) - 1L
  out
}

read_shared_reference_variants <- function(store, bounds = NULL) {
  index_path <- file.path(store$path, store$manifest$files$variants)
  reference <- store_reference_table(store)
  if (is.null(bounds)) {
    index <- arrow::read_parquet(index_path)
    ref <- reference
  } else {
    keep <- reference$chromosome == bounds$chromosome &
      reference$base_pair_location >= bounds$start &
      reference$base_pair_location <= bounds$end
    ref <- reference[keep, , drop = FALSE]
    wanted <- ref$reference_index
    if (length(wanted)) {
      dataset <- arrow::open_dataset(index_path, format = "parquet")
      index <- dplyr::collect(dplyr::filter(dataset, reference_index %in% !!wanted))
    } else {
      index <- data.frame(row = integer(), reference_index = integer())
    }
  }
  if (nrow(index)) {
    hit <- match(index$reference_index, ref$reference_index)
    matched <- !is.na(hit)
    out <- data.frame(row = index$row[matched],
                      chromosome = ref$chromosome[hit[matched]],
                      base_pair_location = ref$base_pair_location[hit[matched]],
                      effect_allele = ref$alternate_allele[hit[matched]],
                      other_allele = ref$reference_allele[hit[matched]],
                      variant_id = ref$variant_id[hit[matched]],
                      rsid = ref$rsid[hit[matched]],
                      stringsAsFactors = FALSE)
  } else {
    out <- data.frame(row = integer(), chromosome = character(),
                      base_pair_location = integer(), effect_allele = character(),
                      other_allele = character(), variant_id = character(),
                      rsid = character(), stringsAsFactors = FALSE)
  }
  unmatched_file <- store$manifest$files$unmatched %||% NULL
  if (!is.null(unmatched_file) && nzchar(unmatched_file) &&
      file.exists(file.path(store$path, unmatched_file))) {
    unmatched_path <- file.path(store$path, unmatched_file)
    unmatched <- if (is.null(bounds)) arrow::read_parquet(unmatched_path) else {
      read_parquet_region(unmatched_path, bounds)
    }
    if (nrow(unmatched)) {
      keep <- intersect(names(out), names(unmatched))
      out <- rbind(out[keep], unmatched[keep])
    }
  }
  identity <- store$manifest$variant_identity$rsid %||% "stored"
  if (!"rsid" %in% names(out) && identical(identity, "derived_from_variant_id")) {
    out$rsid <- out$variant_id
  }
  out[order(out$row), , drop = FALSE]
}

read_variant_table <- function(store, bounds = NULL) {
  if (identical(store$manifest$variant_storage %||% "inline", "shared_reference_index")) {
    return(read_shared_reference_variants(store, bounds = bounds))
  }
  path <- file.path(store$path, store$manifest$files$variants)
  variants <- if (is.null(bounds)) arrow::read_parquet(path) else read_parquet_region(path, bounds)
  identity <- store$manifest$variant_identity$rsid %||% "stored"
  if (!"rsid" %in% names(variants) && identical(identity, "derived_from_variant_id") &&
      "variant_id" %in% names(variants)) {
    variants$rsid <- variants$variant_id
  }
  variants
}

read_standard_values <- function(store, rows = NULL, include_beta = TRUE, include_p = FALSE) {
  m <- store$manifest
  values_path <- file.path(store$path, m$files$values)
  values <- if (is.null(rows)) arrow::read_parquet(values_path) else read_parquet_rows(values_path, rows)
  if (is.null(values)) {
    if (identical(m$profile, "exact")) {
      return(data.frame(row = integer(), z = numeric(), standard_error = numeric(),
                        effect_allele_frequency = numeric()))
    }
    return(data.frame(row = integer(), z_code = integer(), se_code = integer(),
                      eaf_code = integer()))
  }
  if (identical(m$profile, "exact")) {
    return(values)
  }
  exceptions_path <- file.path(store$path, m$files$exceptions)
  exceptions <- if (is.null(rows)) arrow::read_parquet(exceptions_path) else read_parquet_rows(exceptions_path, rows)
  if (is.null(exceptions)) exceptions <- data.frame(row = integer(), z_value = numeric(), se_value = numeric())
  decoded <- q_decode(values[c("row", "z_code", "se_code", "eaf_code")],
                      exceptions, m$codec, include_beta = include_beta,
                      include_p = include_p)
  data.frame(row = values$row, decoded, check.names = FALSE)
}

#' Read decoded summary statistics
#'
#' @param store A store object or path.
#' @param region Optional genomic region.
#' @param variants Optional variant IDs, zero-based row IDs, or—for Pcodec
#'   stores—canonical `chromosome:position:REF:ALT` keys.
#' @param columns Optional output columns.
#' @param use_cache Use an existing q8 cache for a region when available.
#' @return A data.frame with ordinary summary-statistics columns.
#' @examples
#' \dontrun{
#' store <- open_compressor("gwas.cpr")
#' read_sumstats(store, region = "chr1:1000000-2000000",
#'               columns = c("beta", "standard_error", "p_value"))
#' }
#' @export
read_sumstats <- function(store, region = NULL, variants = NULL, columns = NULL, use_cache = FALSE) {
  store <- if (inherits(store, "compressor_store")) store else pcodec_open_store_cached(store)
  if (identical(store$manifest$backend, "pcodec")) {
    if (isTRUE(use_cache)) warning("use_cache is ignored for the block-framed Pcodec backend", call. = FALSE)
    return(pcodec_read_store(store, region = region, variants = variants, columns = columns))
  }
  require_parquet_backend("reading a Parquet CompreSSoR store", dplyr = TRUE)
  cache_columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                     "variant_id", "rsid", "beta", "standard_error",
                     "z", "effect_allele_frequency", "p_value")
  cache_eligible <- is.null(variants) && (is.null(columns) || all(columns %in% cache_columns))
  if (isTRUE(use_cache) && cache_eligible && !is.null(region) && dir.exists(file.path(store$path, "cache.q8"))) {
    cached <- read_q8_cache(store, region = region)
    if (is.null(columns)) return(cached)
    missing <- setdiff(columns, names(cached))
    if (length(missing)) stop("requested columns are not present: ", paste(missing, collapse = ", "), call. = FALSE)
    return(cached[columns])
  }
  m <- store$manifest
  bounds <- read_region_bounds(region)
  variants_data <- read_variant_table(store, bounds = bounds)
  if (!is.null(variants)) {
    if (is.numeric(variants)) variants_data <- variants_data[variants_data$row %in% variants, , drop = FALSE]
    else variants_data <- variants_data[
      variants_data$variant_id %in% as.character(variants) |
        variants_data$rsid %in% as.character(variants), , drop = FALSE]
  }
  if (!nrow(variants_data)) {
    out <- variants_data
    if (!is.null(columns)) {
      missing <- setdiff(columns, names(out))
      if (length(missing)) stop("requested columns are not present: ", paste(missing, collapse = ", "), call. = FALSE)
      return(out[columns])
    }
    return(out)
  }
  need_beta <- is.null(columns) || "beta" %in% columns
  need_p <- is.null(columns) || "p_value" %in% columns
  values <- read_standard_values(store, rows = variants_data$row,
                                 include_beta = need_beta, include_p = need_p)
  out <- merge(variants_data, values, by = "row", sort = FALSE)
  extras_path <- if (!is.null(m$files$extras)) file.path(store$path, m$files$extras) else NULL
  if (!is.null(extras_path) && file.exists(extras_path)) {
    extras <- read_parquet_rows(extras_path, variants_data$row)
    if (!is.null(extras)) out <- merge(out, extras, by = "row", sort = FALSE)
  }
  out <- out[order(out$row), , drop = FALSE]
  if (need_beta && !"z" %in% names(out) && all(c("beta", "standard_error") %in% names(out))) {
    out$z <- out$beta / out$standard_error
  }
  if (need_beta && !"beta" %in% names(out) && all(c("z", "standard_error") %in% names(out))) {
    out$beta <- out$z * out$standard_error
  }
  if (need_p && !"p_value" %in% names(out)) {
    out$p_value <- 2 * stats::pnorm(-abs(out$z))
  }
  out$row <- NULL
  if (is.null(columns)) return(out)
  missing <- setdiff(columns, names(out))
  if (length(missing)) stop("requested columns are not present: ", paste(missing, collapse = ", "), call. = FALSE)
  out[columns]
}

#' Read canonical variants from several compressed GWAS files
#'
#' Starts the Pcodec runtime once for the complete batch. This is substantially
#' faster than repeated [read_sumstats()] calls when an analysis extracts a
#' small instrument set from several GWAS files.
#'
#' @param stores A non-empty list or character vector of Pcodec stores.
#' @param variants A canonical `chromosome:position:REF:ALT` vector shared by
#'   every store, or one such vector per store in a list.
#' @param columns Output columns requested from every store.
#' @param threads Number of stores to decode concurrently inside the shared
#'   Pcodec process.
#' @return A list of decoded data frames in the same order as `stores`.
#' @examples
#' \dontrun{
#' read_sumstats_batch(
#'   c(exposure = "exposure.cpr", outcome = "outcome.cpr"),
#'   c("1:100000:A:G", "1:200000:C:T"),
#'   columns = c("chromosome", "base_pair_location", "effect_allele",
#'               "other_allele", "beta", "standard_error")
#' )
#' }
#' @export
read_sumstats_batch <- function(
    stores,
    variants,
    columns = c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error"),
    threads = 1L) {
  if (is.character(stores)) stores <- as.list(stores)
  if (!is.list(stores) || !length(stores)) {
    stop("stores must be a non-empty list or character vector", call. = FALSE)
  }
  store_names <- names(stores)
  if (is.character(variants)) {
    variants <- rep(list(variants), length(stores))
  }
  if (!is.list(variants) || length(variants) != length(stores)) {
    stop("variants must be a canonical-key vector or one list element per store",
         call. = FALSE)
  }
  result <- pcodec_read_stores(
    stores, variants, unique(as.character(columns)), threads = threads
  )
  if (!is.null(store_names)) names(result) <- store_names
  result
}

#' Decompress a CompreSSoR store
#'
#' @inheritParams read_sumstats
#' @examples
#' \dontrun{
#' decoded <- decompress_sumstats("gwas.cpr")
#' }
#' @export
decompress_sumstats <- function(store, region = NULL, variants = NULL, columns = NULL, use_cache = FALSE) {
  read_sumstats(store, region = region, variants = variants, columns = columns, use_cache = use_cache)
}

#' Validate a CompreSSoR store
#'
#' @param store A store object or path.
#' @param full For Pcodec stores, decode and semantically check every frame in
#'   addition to verifying every file, frame, index, and manifest checksum.
#' @return A list with `valid`, `errors`, `rows` and `profile`.
#' @examples
#' \dontrun{
#' validate_compressor("gwas.cpr", full = TRUE)
#' }
#' @export
validate_compressor <- function(store, full = FALSE) {
  store_object <- tryCatch(
    if (inherits(store, "compressor_store")) store else open_compressor(store),
    error = function(e) NULL
  )
  if (!is.null(store_object) && identical(store_object$manifest$backend, "pcodec")) {
    return(pcodec_validate_store(store_object, full = full))
  }
  errors <- character()
  result <- tryCatch({
    s <- if (inherits(store, "compressor_store")) store else open_compressor(store)
    m <- s$manifest
    file_names <- unlist(m$files, use.names = FALSE)
    file_names <- file_names[!is.na(file_names) & nzchar(file_names)]
    missing <- file_names[!file.exists(file.path(s$path, file_names))]
    if (length(missing)) errors <<- c(errors, paste("missing", missing))
    variants <- if (!length(errors)) arrow::read_parquet(file.path(s$path, m$files$variants)) else NULL
    values <- if (!length(errors)) read_standard_values(s) else NULL
    if (!is.null(variants) && nrow(variants) != m$n_rows) errors <<- c(errors, "variant row count does not match manifest")
    if (!is.null(values) && nrow(values) != m$n_rows) errors <<- c(errors, "value row count does not match manifest")
    if (!is.null(variants) && anyDuplicated(variants$row)) errors <<- c(errors, "variant row IDs are duplicated")
    if (!is.null(values) && anyDuplicated(values$row)) errors <<- c(errors, "value row IDs are duplicated")
    if (!is.null(variants) && !is.null(values) && !identical(sort(variants$row), sort(values$row))) {
      errors <<- c(errors, "variant and value row IDs do not agree")
    }
    extras_file <- m$files$extras
    if (!is.null(extras_file) && length(extras_file) && !is.na(extras_file) && nzchar(extras_file)) {
      extras <- if (!length(errors)) arrow::read_parquet(file.path(s$path, extras_file)) else NULL
      if (!is.null(extras) && nrow(extras) != m$n_rows) errors <<- c(errors, "extra row count does not match manifest")
      if (!is.null(extras) && anyDuplicated(extras$row)) errors <<- c(errors, "extra row IDs are duplicated")
      if (!is.null(extras) && !is.null(values) && !identical(sort(extras$row), sort(values$row))) {
        errors <<- c(errors, "extra and value row IDs do not agree")
      }
    }
    list(valid = !length(errors), errors = errors, rows = m$n_rows, profile = m$profile)
  }, error = function(e) list(valid = FALSE, errors = conditionMessage(e), rows = NA_integer_, profile = NA_character_))
  result
}
