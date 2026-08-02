#' Write a CompreSSoR store
#'
#' @param input A data.frame or a delimited summary-statistics file.
#' @param output Destination directory for the store.
#' @param reference GRCh38 reference spine, "GRCh38" for the cached/downloaded
#'   default, or NULL to skip reference alignment.
#' @param mode One of `"qc"` (default, harmonise and preserve all rows),
#'   `"convert"` (minimal conversion without reference QC), `"all"` (alias
#'   for `"qc"`), `"core"`, or `"hm3"`.
#' @param variant_set Panel data.frame or file used by `mode = "core"` or
#'   `mode = "hm3"`. PLINK `.bim`, Parquet and delimited panel files are
#'   supported.
#' @param strict If `TRUE`, fail when variants are absent, incompatible, or
#'   ambiguous. The default preserves such rows and records their status in
#'   the extra-column sidecar and manifest.
#' @param chrom_threads Number of chromosome workers for harmonisation. Values
#'   above one use chromosome-parallel harmonisation and share the reference.
#' @param profile `"standard"` uses q9/SE10 compact Parquet; `"exact"`
#'   stores the numeric values without quantisation in Parquet.
#' @param overwrite Whether an existing destination may be replaced.
#' @param cache Whether to build the optional q8 framed cache after writing.
#' @param block_rows Number of rows per Parquet row group.
#' @return A `compressor_store` object.
#' @export
compress_sumstats <- function(input, output, reference = "GRCh38",
                              profile = c("standard", "exact"),
                              overwrite = FALSE, cache = FALSE, strict = FALSE,
                              block_rows = 65536L,
                              mode = c("qc", "convert", "all", "core", "hm3"),
                              variant_set = NULL, chrom_threads = 1L) {
  profile <- match.arg(profile)
  if (length(output) != 1L || !is.character(output)) stop("output must be one directory path", call. = FALSE)
  if (file.exists(output)) {
    if (!overwrite) stop("output already exists; use overwrite = TRUE", call. = FALSE)
    unlink(output, recursive = TRUE, force = TRUE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output, "extras"), showWarnings = FALSE)

  raw <- normalise_sumstats_columns(read_sumstats_input(input))
  prepared <- prepare_sumstats_data(raw, reference, mode = mode, variant_set = variant_set,
                                    strict = strict, chrom_threads = chrom_threads)
  alignment <- prepared$alignment
  data <- prepared$data
  if (!identical(prepared$effective_mode, "convert")) validate_sumstats_values(data)
  n <- nrow(data)

  variant_names <- intersect(c("row", "chromosome", "base_pair_location", "effect_allele",
                               "other_allele", "variant_id", "rsid"), names(data))
  variants <- data[intersect(variant_names, names(data))]
  variants$row <- seq_len(n) - 1L
  variants <- variants[c("row", setdiff(names(variants), "row"))]
  arrow::write_parquet(variants, file.path(output, "variants.parquet"),
                       compression = "zstd", compression_level = 7,
                       write_statistics = TRUE, use_dictionary = TRUE,
                       chunk_size = as.integer(block_rows))

  logical_columns <- c("beta", "standard_error", "effect_allele_frequency", "p_value")
  if (profile == "exact") {
    values <- data[logical_columns]
    values$row <- seq_len(n) - 1L
    values <- values[c("row", logical_columns)]
    arrow::write_parquet(values, file.path(output, "values.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, chunk_size = as.integer(block_rows))
    codec <- list(name = "lossless_numeric", z_bits = NULL, se_bits = NULL, eaf_bits = NULL)
    files <- list(variants = "variants.parquet", values = "values.parquet", exceptions = NULL, extras = NULL)
  } else {
    se_meta <- q_profile_metadata(data$standard_error)
    encoded <- q_encode(data$beta, data$standard_error, data$effect_allele_frequency,
                        z_bits = 9L, se_bits = 10L,
                        se_lo = se_meta$se_log_min, se_hi = se_meta$se_log_max)
    encoded$main$p_value <- data$p_value
    arrow::write_parquet(encoded$main, file.path(output, "values.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, chunk_size = as.integer(block_rows))
    arrow::write_parquet(encoded$exceptions, file.path(output, "exceptions.parquet"),
                         compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, chunk_size = as.integer(block_rows))
    codec <- c(list(name = "q9_z_se10_eaf12", compression = "zstd"), encoded$metadata)
    files <- list(variants = "variants.parquet", values = "values.parquet", exceptions = "exceptions.parquet", extras = NULL)
  }

  core_names <- unique(c(variant_names, logical_columns))
  extra_names <- setdiff(names(data), core_names)
  extra_names <- extra_names[vapply(data[extra_names], function(x) is.atomic(x) && !is.list(x), logical(1))]
  if (length(extra_names)) {
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
      alignment = alignment$alignment_stats
    ),
    n_rows = n,
    block_rows = as.integer(block_rows),
    logical_columns = logical_columns,
    source_columns = names(raw),
    codec = codec,
    files = files,
    benchmark = benchmark_metadata(),
    benchmark_comparisons = list(
      vcf_tabix = vcf_benchmark_metadata(),
      finngen_end_to_end = finngen_benchmark_metadata(),
      finngen_optimization = finngen_optimization_metadata(),
      modes_edge = mode_benchmark_metadata(),
      release_gate = release_gate_benchmark_metadata(),
      release_gate_followup = release_gate_followup_metadata()
    ),
    created_utc = now_utc(),
    tolerances = if (profile == "standard") list(eaf_abs_max = 0.0002, z_central_range = c(-3.5, 3.5)) else list(exact = TRUE)
  )
  write_manifest(manifest, file.path(output, "manifest.json"))
  store <- open_compressor(output)
  if (isTRUE(cache)) build_cache(store, overwrite = TRUE, block_rows = block_rows)
  store
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

read_variant_table <- function(store, bounds = NULL) {
  path <- file.path(store$path, store$manifest$files$variants)
  if (is.null(bounds)) arrow::read_parquet(path) else read_parquet_region(path, bounds)
}

read_standard_values <- function(store, rows = NULL) {
  m <- store$manifest
  values_path <- file.path(store$path, m$files$values)
  values <- if (is.null(rows)) arrow::read_parquet(values_path) else read_parquet_rows(values_path, rows)
  if (is.null(values)) {
    if (identical(m$profile, "exact")) {
      return(data.frame(row = integer(), beta = numeric(), standard_error = numeric(),
                        effect_allele_frequency = numeric(), p_value = numeric()))
    }
    return(data.frame(row = integer(), z_code = integer(), se_code = integer(),
                      eaf_code = integer(), p_value = numeric()))
  }
  if (identical(m$profile, "exact")) {
    names(values)[names(values) == "beta"] <- "beta"
    return(values)
  }
  exceptions_path <- file.path(store$path, m$files$exceptions)
  exceptions <- if (is.null(rows)) arrow::read_parquet(exceptions_path) else read_parquet_rows(exceptions_path, rows)
  if (is.null(exceptions)) exceptions <- data.frame(row = integer(), z_value = numeric(), se_value = numeric())
  decoded <- q_decode(values[c("row", "z_code", "se_code", "eaf_code")], exceptions, m$codec)
  data.frame(row = values$row, decoded, p_value = values$p_value, check.names = FALSE)
}

#' Read decoded summary statistics
#'
#' @param store A store object or path.
#' @param region Optional genomic region.
#' @param variants Optional variant IDs or zero-based row IDs.
#' @param columns Optional output columns.
#' @param use_cache Use an existing q8 cache for a region when available.
#' @return A data.frame with ordinary summary-statistics columns.
#' @export
read_sumstats <- function(store, region = NULL, variants = NULL, columns = NULL, use_cache = FALSE) {
  store <- if (inherits(store, "compressor_store")) store else open_compressor(store)
  cache_columns <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                     "variant_id", "rsid", "beta", "standard_error",
                     "effect_allele_frequency", "p_value")
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
  values <- read_standard_values(store, rows = variants_data$row)
  out <- merge(variants_data, values, by = "row", sort = FALSE)
  extras_path <- if (!is.null(m$files$extras)) file.path(store$path, m$files$extras) else NULL
  if (!is.null(extras_path) && file.exists(extras_path)) {
    extras <- read_parquet_rows(extras_path, variants_data$row)
    if (!is.null(extras)) out <- merge(out, extras, by = "row", sort = FALSE)
  }
  out <- out[order(out$row), , drop = FALSE]
  out$row <- NULL
  if (is.null(columns)) return(out)
  missing <- setdiff(columns, names(out))
  if (length(missing)) stop("requested columns are not present: ", paste(missing, collapse = ", "), call. = FALSE)
  out[columns]
}

#' Decompress a CompreSSoR store
#'
#' @inheritParams read_sumstats
#' @export
decompress_sumstats <- function(store, region = NULL, variants = NULL, columns = NULL, use_cache = FALSE) {
  read_sumstats(store, region = region, variants = variants, columns = columns, use_cache = use_cache)
}

#' Validate a CompreSSoR store
#'
#' @param store A store object or path.
#' @return A list with `valid`, `errors`, `rows` and `profile`.
#' @export
validate_compressor <- function(store) {
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
