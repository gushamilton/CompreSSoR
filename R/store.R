compressor_manifest_contract <- function(path, build, selection, profile,
                                         backend, threads, row_policy, source,
                                         preparation = NULL, panel = NULL) {
  manifest_path <- file.path(path, "manifest.json")
  manifest <- read_manifest(manifest_path)
  manifest$genome_build <- build
  manifest$input_build <- build
  manifest$stored_build <- build
  manifest$assembly <- build
  manifest$selection_scope <- selection
  manifest$profile <- manifest$profile %||% profile
  manifest$backend <- manifest$backend %||% backend
  manifest$codec <- manifest$codec %||% list(name = profile)
  manifest$tolerances <- manifest$tolerances %||% list()
  writer <- manifest$writer %||% list()
  effective_writer_workers <- as.integer(writer$effective_workers %||% 1L)
  manifest$threads <- list(
    requested = as.integer(threads),
    effective = effective_writer_workers,
    compression = as.integer(threads),
    compression_writer_parallel = effective_writer_workers > 1L
  )
  manifest$row_policy <- row_policy
  manifest$provenance <- list(
    input = source %||% NULL,
    input_build = build,
    stored_build = build,
    reference = list(status = "not_used"),
    chain = list(status = "not_used"),
    panel = panel %||% NULL,
    preparation = preparation %||% list(method = "strict_prepared_input")
  )
  if (is.null(manifest$selection)) {
    manifest$selection <- list(name = selection, method = "full_store")
  }
  manifest$metadata <- c(
    manifest$metadata %||% list(),
    list(api = list(selection = selection, store_build = build,
                    input_build = build, threads = as.integer(threads),
                    row_policy = row_policy))
  )
  write_manifest(manifest, manifest_path)
  if (identical(backend, "pcodec")) seal_pcodec_manifest(manifest_path)
  invisible(manifest)
}

record_commit_timing <- function(path, seconds) {
  manifest_path <- file.path(path, "manifest.json")
  manifest <- read_manifest(manifest_path)
  timings <- manifest$timings %||% list(unit = "seconds", phases = list())
  timings$unit <- timings$unit %||% "seconds"
  timings$phases <- timings$phases %||% list()
  timings$phases$commit <- as.numeric(seconds)
  manifest$timings <- timings
  write_manifest(manifest, manifest_path)
  if (identical(manifest$backend, "pcodec")) seal_pcodec_manifest(manifest_path)
  invisible(manifest)
}

write_selection_regions <- function(output, selection) {
  if (is.null(selection)) return(NULL)
  regions <- selection$regions_table
  if (!is.data.frame(regions)) {
    stop("selection region table is missing or malformed", call. = FALSE)
  }
  file_name <- if (identical(selection$tag %||% "core", "core_plus")) {
    "core_plus_regions.json"
  } else {
    "core_regions.json"
  }
  payload <- list(
    format = "CompreSSoR-core-regions",
    version = 1L,
    tag = selection$tag %||% "core",
    method = selection$method %||% "pvalue_regions",
    pvalue_threshold = selection$pvalue_threshold,
    padding_bp = selection$padding_bp,
    regions = lapply(seq_len(nrow(regions)), function(index) {
      list(chromosome = as.character(regions$chromosome[[index]]),
           start = as.integer(regions$start[[index]]),
           end = as.integer(regions$end[[index]]),
           seed_snps = as.integer(regions$seed_snps[[index]]))
    })
  )
  jsonlite::write_json(payload, file.path(output, file_name),
                       auto_unbox = TRUE, pretty = FALSE, digits = 17)
  selection$regions_table <- NULL
  selection$file <- file_name
  selection
}

#' Write a CompreSSoR store
#'
#' @param input A data.frame or a delimited summary-statistics file.
#' @param output Destination directory for the store.
#' @param selection Storage scope: `"full"`, `"core"`, `"hm3"`, or
#'   `"core_plus"`. Selection is applied after the strict input contract and
#'   structural QC; it never changes identity or allele orientation.
#' @param store_build Build encoded by the store. It must normalize to the
#'   same build as `input_build`; cross-build conversion is rejected.
#' @param threads Positive worker count requested for compression.
#' @param row_policy `"report"` drops structurally invalid rows and records
#'   aggregate counts; `"error"` rejects them. `strict = TRUE` is an alias.
#' @param variant_set Optional prevalidated panel used only for deterministic
#'   post-import selection; it is never used as a reference.
#' @param strict If `TRUE`, fail when variants are absent, incompatible,
#'   ambiguous, duplicated, or unsupported by the selected backend.
#' @param input_build Input build, explicitly GRCh37/hg19 or GRCh38/hg38.
#' @param pvalue_threshold Strict p-value threshold for `selection = "pvalue_regions"`
#'   or `selection = "core_plus"`; p-values are derived from canonical Z.
#' @param region_padding Number of base pairs added on each side of significant
#'   SNPs for `selection = "pvalue_regions"` or `selection = "core_plus"`.
#' @param profile `"standard"` uses semantic Z9/EAF8/SE6 streams with sparse
#'   float32 exceptions; `"exact"` is available with the Parquet backend. P and
#'   beta are derived rather than stored.
#' @param backend Storage backend. The default, `"pcodec"`, is the compact
#'   block-framed Pcodec format with a self-contained build-aware identity key.
#'   `"parquet"` retains the interoperable legacy backend.
#' @param keep_extras Whether arbitrary non-core input columns should be kept
#'   in an extras sidecar. The default is `FALSE`; aggregate QC counts remain
#'   in the manifest.
#' @param overwrite Whether an existing destination may be replaced.
#' @param cache Whether to build the optional q8 framed cache after writing.
#' @param block_rows Number of rows per Parquet row group; Pcodec uses its
#'   measured fixed key and value frame sizes.
#' @return A `compressor_store` object.
#' @examples
#' \dontrun{
#' example <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
#' store <- compress_sumstats(
#'   example, file.path(tempdir(), "example.cpr"),
#'   input_build = "GRCh38", store_build = "GRCh38", overwrite = TRUE
#' )
#' read_sumstats(store, columns = c("chromosome", "base_pair_location", "z"))
#' }
#' @export
compress_sumstats <- function(input, output,
                              profile = c("standard", "exact"),
                              overwrite = FALSE, cache = FALSE, strict = NULL,
                              keep_extras = FALSE,
                              block_rows = 65536L,
                              variant_set = NULL, input_build = "GRCh38",
                              backend = c("pcodec", "parquet"),
                              pvalue_threshold = 1e-5,
                              region_padding = 50000L,
                              store_build = "GRCh38",
                              selection = c("full", "core", "hm3", "core_plus"),
                              threads = NULL,
                              row_policy = c("report", "error")) {
  builds <- core_builds(input_build, store_build)
  input_build <- builds$input_build
  store_build <- builds$store_build
  if (missing(selection)) selection <- NULL
  if (is.null(threads)) {
    threads <- 1L
  } else {
    threads <- pcodec_validate_threads(threads)
  }
  if (missing(row_policy) || is.null(row_policy)) {
    row_policy <- if (isTRUE(strict)) "error" else "report"
  } else {
    row_policy <- match.arg(row_policy, c("report", "error"))
    if (!is.null(strict) &&
        (!is.logical(strict) || length(strict) != 1L || is.na(strict) ||
         isTRUE(strict) != identical(row_policy, "error"))) {
      stop("row_policy and strict disagree; use row_policy alone", call. = FALSE)
    }
  }
  if (is.null(strict)) strict <- identical(row_policy, "error")
  if (!is.null(selection)) {
    if (length(selection) != 1L || is.na(selection) ||
        !selection %in% c("full", "core", "hm3", "core_plus", "pvalue_regions")) {
      stop("selection must be one of full, core, hm3, core_plus, or pvalue_regions", call. = FALSE)
    }
  } else {
    selection <- "full"
  }
  profile <- match.arg(profile)
  backend <- match.arg(backend)
  if (identical(backend, "pcodec") && isTRUE(keep_extras)) {
    stop("backend='pcodec' currently stores only the core summary-statistics "
         , "columns; use backend='parquet' to retain extras", call. = FALSE)
  }
  if (selection %in% c("core", "hm3", "core_plus") && is.null(variant_set)) {
    variant_set <- if (selection == "hm3") "hm3" else "core"
  }
  selection_scope <- selection
  if (identical(backend, "parquet")) {
    require_parquet_backend("backend='parquet'", dplyr = TRUE)
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

  raw <- import_sumstats_impl(
    input, input_build = input_build,
    project_columns = identical(backend, "pcodec") || !isTRUE(keep_extras),
    core_only = identical(backend, "pcodec") || !isTRUE(keep_extras),
    allow_p_to_se = FALSE,
    run_qc = FALSE
  )
  source_columns <- attr(raw, "source_columns")
  source_columns_read <- attr(raw, "source_columns_read")
  source_provenance <- attr(raw, "source_provenance")
  phase_timings <- attr(raw, "phase_timings") %||%
    list(unit = "seconds", phases = list())
  phase_timings$unit <- phase_timings$unit %||% "seconds"
  phase_timings$phases <- phase_timings$phases %||% list()
  qc_started <- phase_clock()
  validate_core_schema(raw, source_columns = source_columns,
                       input_build = input_build, store_build = store_build)
  raw <- validate_core_orientation(raw, row_policy = row_policy)
  structural <- apply_structural_qc(
    raw, input_build = input_build, strict = strict,
    row_policy = row_policy, require_statistics = TRUE,
    drop_duplicates = TRUE, check_duplicates = TRUE, detail = "compact"
  )
  phase_timings$phases$qc <- phase_seconds(qc_started)
  raw <- structural$data
  if (structural$report$input_rows > 0L && !nrow(raw)) {
    failure <- format_structural_qc_failure(structural$report)
    orientation_rejections <- structural$report$rejection_counts[["orientation_mismatch"]] %||% 0L
    if (orientation_rejections > 0L) {
      failure <- paste(
        "effect_allele and other_allele are inconsistent with explicit REF/ALT;",
        failure
      )
    }
    stop(failure, call. = FALSE)
  }
  identity_started <- phase_clock()
  raw <- canonicalize_core_identity(raw, build = store_build)
  phase_timings$phases$identity_sort <- phase_seconds(identity_started)
  attr(raw, "source_columns") <- source_columns
  attr(raw, "source_columns_read") <- source_columns_read
  attr(raw, "source_provenance") <- source_provenance
  prepared <- prepare_core_sumstats_data(
    raw, selection = selection, variant_set = variant_set,
    pvalue_threshold = pvalue_threshold, region_padding = region_padding,
    build = store_build
  )
  preparation <- prepared$preparation
  preparation$structural_qc <-
    compact_structural_qc_report(structural$report)
  data <- prepared$data
  selection <- prepared$selection
  validate_sumstats_values(data, require_identity = TRUE)
  if (identical(backend, "pcodec")) {
    if (!identical(profile, "standard")) {
      stop("backend='pcodec' currently provides the standard semantic profile; "
           , "use backend='parquet' for profile='exact'", call. = FALSE)
    }
    if (isTRUE(cache)) {
      stop("backend='pcodec' is already block-framed for regional access; "
           , "cache=TRUE is not needed", call. = FALSE)
    }
    supported_chromosome <- as.character(data$chromosome) %in%
      names(compressor_chromosome_lengths(store_build))
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
      stats <- preparation %||% list()
      stats$unsupported_pcodec_rows <- as.integer(unsupported_rows)
      stats$dropped_unsupported_pcodec <- as.integer(unsupported_rows)
      preparation <- stats
    }
    if (!nrow(data)) {
      stop("no supported biallelic primary-chromosome SNVs remain for the Pcodec store", call. = FALSE)
    }
    reference_metadata <- list(id = "none", build = store_build, status = "not_used",
                               external_reference_required = FALSE)
    pcodec_metadata <- list(
      block_rows = as.integer(block_rows),
      genome_build = store_build,
      identity = compressor_identity_manifest(
        input_build = input_build, stored_build = store_build
      ),
      reference = reference_metadata,
      preparation = list(
        method = "strict_prepared_input",
        threads_requested = as.integer(threads),
        row_policy = row_policy,
        strict = isTRUE(strict),
        input_build = as.character(input_build),
        preparation = preparation
      ),
      source_columns = attr(raw, "source_columns") %||% names(raw),
      source_columns_read = attr(raw, "source_columns_read") %||% names(raw),
      source = attr(raw, "source_provenance") %||% NULL,
      selection = selection,
      timings = phase_timings
    )
    pcodec_write_store(data, output, metadata = pcodec_metadata)
    compressor_manifest_contract(
      output, build = store_build, selection = selection_scope,
      profile = profile, backend = backend,
      threads = threads, row_policy = row_policy,
      source = attr(raw, "source_provenance") %||% NULL,
      preparation = preparation,
      panel = preparation$variant_set %||% NULL
    )
    commit_started <- phase_clock()
    commit_store_output(transaction)
    commit_seconds <- phase_seconds(commit_started)
    completed <- TRUE
    record_commit_timing(transaction$target, commit_seconds)
    return(open_compressor(transaction$target))
  }

  encode_started <- phase_clock()
  n <- nrow(data)
  if (!is.null(selection)) {
    selection$kept_rows <- as.integer(n)
    selection$dropped_rows <- as.integer(selection$input_rows - n)
    selection <- write_selection_regions(output, selection)
  }

  variant_names <- intersect(c("row", "chromosome", "base_pair_location",
                               "reference_allele", "alternate_allele",
                               "effect_allele", "other_allele", "variant_id", "rsid"), names(data))
  if (all(c("variant_id", "rsid") %in% variant_names)) {
    same_identity <- (is.na(data$variant_id) & is.na(data$rsid)) |
      (!is.na(data$variant_id) & !is.na(data$rsid) &
         as.character(data$variant_id) == as.character(data$rsid))
    if (all(same_identity)) variant_names <- setdiff(variant_names, "rsid")
  }
  variants <- data[intersect(variant_names, names(data))]
  variants$row <- seq_len(n) - 1L
  variants <- variants[c("row", setdiff(names(variants), "row"))]
  files <- list(variants = "variants.parquet", values = "values.parquet",
                exceptions = NULL, unmatched = NULL, extras = NULL,
                selection_regions = if (is.null(selection)) NULL else selection$file)
  arrow::write_parquet(variants, file.path(output, "variants.parquet"),
                       compression = "zstd", compression_level = 7,
                       write_statistics = TRUE, use_dictionary = TRUE,
                       chunk_size = as.integer(block_rows))

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

  # Structural QC is performed before writing and its aggregate counts are
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
  phase_timings$phases$encode <- phase_seconds(encode_started)

  reference_manifest <- list(id = "none", build = store_build, status = "not_used",
                             external_reference_required = FALSE)
  manifest <- list(
    format = "CompreSSoR",
    format_version = "0.1.0",
    profile = profile,
    backend = "parquet",
    genome_build = store_build,
    identity = compressor_identity_manifest(
      input_build = input_build, stored_build = store_build
    ),
    reference = reference_manifest,
    preparation = list(
      method = "strict_prepared_input",
      threads_requested = as.integer(threads),
      row_policy = row_policy,
      strict = isTRUE(strict),
      input_build = as.character(input_build),
      preparation = preparation
    ),
    n_rows = n,
    block_rows = as.integer(block_rows),
    logical_columns = logical_columns,
    derived_columns = list(beta = "z * standard_error", p_value = "2 * pnorm(-abs(z))"),
    source_columns = attr(raw, "source_columns") %||% names(raw),
    source_columns_read = attr(raw, "source_columns_read") %||% names(raw),
    source = attr(raw, "source_provenance") %||% NULL,
    codec = codec,
    variant_storage = "inline",
    variant_identity = list(
      rsid = if ("rsid" %in% variant_names) "stored" else "derived_from_variant_id"
    ),
    selection = selection,
    files = files,
    timings = phase_timings,
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
  compressor_manifest_contract(
    output, build = store_build, selection = selection_scope,
    profile = profile, backend = backend,
    threads = threads, row_policy = row_policy,
    source = attr(raw, "source_provenance") %||% NULL,
    preparation = preparation,
    panel = preparation$variant_set %||% NULL
  )
  store <- open_compressor(output)
  if (isTRUE(cache)) build_cache(store, overwrite = TRUE, block_rows = block_rows)
  commit_started <- phase_clock()
  commit_store_output(transaction)
  commit_seconds <- phase_seconds(commit_started)
  completed <- TRUE
  record_commit_timing(transaction$target, commit_seconds)
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
      !isTRUE(manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS)) {
    stop("unsupported or archived Pcodec format version: ",
         manifest$format_version %||% "missing",
         "; this build reads native 0.4 stores only", call. = FALSE)
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

read_variant_table <- function(store, bounds = NULL) {
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
#' @param columns Optional output columns. Pcodec stores additionally expose
#'   `global_position` and `substitution`, the compact identity fields used by
#'   native analytical readers; requesting these avoids chromosome/REF/ALT
#'   reconstruction.
#' @param threads Number of workers for Pcodec reads. The native full-store
#'   reader performs indexed stream scans and block decompression in C++; the
#'   regional, canonical-key and fallback paths use Unix worker processes.
#'   Set `threads` explicitly, or use `options(CompreSSoR.pcodec.threads = n)`
#'   for the non-native paths. On Windows, reads remain serial.
#' @param use_cache Use an existing q8 cache for a region when available.
#' @return A data.frame with ordinary summary-statistics columns.
#' @examples
#' \dontrun{
#' store <- open_compressor("gwas.cpr")
#' read_sumstats(store, region = "chr1:1000000-2000000",
#'               columns = c("beta", "standard_error", "p_value"))
#' }
#' @export
read_sumstats <- function(store, region = NULL, variants = NULL, columns = NULL,
                          use_cache = FALSE, threads = NULL) {
  store <- if (inherits(store, "compressor_store")) store else pcodec_open_store_cached(store)
  if (identical(store$manifest$backend, "pcodec")) {
    if (isTRUE(use_cache)) warning("use_cache is ignored for the block-framed Pcodec backend", call. = FALSE)
    return(pcodec_read_store(store, region = region, variants = variants,
                             columns = columns, threads = threads))
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
#' Reads the requested canonical keys from several GWAS files. With
#' `threads > 1` on Unix-like systems, independent stores are decoded in
#' parallel; this is substantially faster than repeated [read_sumstats()]
#' calls when an analysis extracts a small instrument set from several files.
#'
#' @param stores A non-empty list or character vector of Pcodec stores.
#' @param variants A canonical `chromosome:position:REF:ALT` vector shared by
#'   every store, or one such vector per store in a list.
#' @param columns Output columns requested from every store.
#' @param threads Number of independent Pcodec store readers to run in
#'   parallel on Unix-like systems. The default is one. Windows uses serial
#'   reads for portability.
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
decompress_sumstats <- function(store, region = NULL, variants = NULL, columns = NULL,
                                use_cache = FALSE, threads = NULL) {
  read_sumstats(store, region = region, variants = variants, columns = columns,
                threads = threads, use_cache = use_cache)
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
