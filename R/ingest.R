#' Import GWAS summary statistics
#'
#' Reads a data frame, delimited text file, gzip-compressed text file, or
#' single-ALT VCF and maps common GWAS column names onto CompreSSoR's canonical
#' in-memory schema. Unrecognised columns are retained. It does not change
#' genome build or allele orientation; the compression core requires callers
#' to provide both explicitly.
#'
#' @param input A data.frame or an existing summary-statistics file.
#' @param strict Whether structural QC failures should stop import. The
#'   default is FALSE, which keeps canonical rows and attaches a bounded
#'   structural_qc_report attribute for the caller to apply.
#' @param row_policy Either report (attach a report and continue) or error
#'   (stop on any structural rejection). strict = TRUE is an alias for error.
#' @param input_build Build used for coordinate-bound preflight, either
#'   GRCh37/hg19 or GRCh38/hg38.
#' @param allow_p_to_se Explicitly opt in to conflict-checked p-value-to-SE
#'   conversion at the boundary. The strict compression core always leaves
#'   this disabled.
#' @return A data.frame containing canonical columns such as `chromosome`,
#'   `base_pair_location`, `reference_allele`, `alternate_allele`,
#'   `effect_allele`, `other_allele`, `beta`, `standard_error`, `z`,
#'   `effect_allele_frequency`, and `p_value`.
#' @export
import_sumstats <- function(input, strict = FALSE,
                            row_policy = c("report", "error"),
                            input_build = "GRCh38", allow_p_to_se = FALSE) {
  import_sumstats_impl(
    input, strict = strict, row_policy = row_policy,
    input_build = input_build, project_columns = FALSE,
    allow_p_to_se = allow_p_to_se
  )
}

import_sumstats_impl <- function(input, strict = FALSE,
                                 row_policy = c("report", "error"),
                                 input_build = "GRCh38",
                                 project_columns = FALSE,
                                 core_only = FALSE,
                                 allow_p_to_se = FALSE,
                                 run_qc = TRUE,
                                 prepared_core = FALSE,
                                 construct_variant_id = TRUE,
                                 qc_detail = c("full", "compact"),
                                 include_p_value = FALSE,
                                 hash_threads = 1L,
                                 allele_columns = NULL) {
  if (length(strict) != 1L || !is.logical(strict) || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  row_policy <- match.arg(row_policy)
  qc_detail <- match.arg(qc_detail)
  if (isTRUE(strict)) row_policy <- "error"
  read_started <- phase_clock()
  allele_columns <- validate_allele_columns(allele_columns)
  raw <- read_sumstats_input(
    input,
    parse_policy = if (identical(row_policy, "error")) "error" else "report",
    project_columns = project_columns,
    core_only = core_only,
    allow_p_to_se = allow_p_to_se,
    include_p_value = include_p_value,
    extra_columns = unname(allele_columns)
  )
  read_seconds <- phase_seconds(read_started)
  if (!nrow(raw)) {
    stop("input contains zero rows; CompreSSoR requires at least one prepared summary-statistics row",
         call. = FALSE)
  }
  source_columns <- attr(raw, "source_columns") %||% names(raw)
  source_columns_read <- attr(raw, "source_columns_read") %||% names(raw)
  input_read_metadata <- attr(raw, "input_read_metadata") %||% list()
  provenance <- if (is.data.frame(input)) {
    list(
      kind = "data.frame",
      rows = nrow(input),
      # Content hash of the columns read (source_columns_read), not of the
      # whole serialized object; see columns_content_hash().
      hash_algorithm = COMPRESSOR_COLUMNS_HASH_ALGORITHM,
      hash = columns_content_hash(raw, threads = hash_threads),
      columns_before = as.integer(length(source_columns)),
      columns_read = as.integer(length(source_columns_read)),
      projected = isTRUE(input_read_metadata$projected),
      projection_elapsed_seconds = as.numeric(
        input_read_metadata$projection_elapsed_seconds %||% 0
      ),
      read_elapsed_seconds = as.numeric(input_read_metadata$elapsed_seconds %||% NA_real_)
    )
  } else {
    path <- normalizePath(input, mustWork = TRUE)
    list(
      kind = "file",
      file = basename(path),
      bytes = unname(file.info(path)$size),
      sha256 = digest::digest(path, algo = "sha256", file = TRUE),
      columns_before = as.integer(length(source_columns)),
      columns_read = as.integer(length(source_columns_read)),
      projected = isTRUE(input_read_metadata$projected),
      projection_elapsed_seconds = as.numeric(
        input_read_metadata$projection_elapsed_seconds %||% 0
      ),
      read_elapsed_seconds = as.numeric(input_read_metadata$elapsed_seconds %||% NA_real_)
    )
  }
  if (!is.null(allele_columns)) {
    raw <- apply_allele_columns(raw, allele_columns)
    provenance$allele_columns <- list(
      reference_allele = unname(allele_columns[["ref"]]),
      alternate_allele = unname(allele_columns[["alt"]])
    )
  }
  normalise_started <- phase_clock()
  out <- if (isTRUE(prepared_core)) {
    normalise_prepared_core_columns(raw, input_build = input_build)
  } else {
    normalise_sumstats_columns(
      raw, parse_policy = if (identical(row_policy, "error")) "error" else "report",
      allow_p_to_se = allow_p_to_se,
      construct_variant_id = construct_variant_id
    )
  }
  normalise_seconds <- phase_seconds(normalise_started)
  provenance$resolution <- attr(out, "resolution_provenance") %||%
    sumstats_resolution_contract(allow_p_to_se = allow_p_to_se)
  provenance$p_value <- list(
    column_present = isTRUE(attr(out, "p_value_source_present")),
    source_alias = attr(out, "p_value_source_alias") %||% "absent",
    selection_policy = "finite_supplied_authoritative_else_pre_encoding_z_fallback"
  )
  provenance$eaf <- eaf_coverage_metadata(out$effect_allele_frequency)
  qc_report <- NULL
  qc_seconds <- NULL
  if (isTRUE(run_qc)) {
    qc_started <- phase_clock()
    qc_report <- structural_qc_report(out, input_build = input_build,
                                      require_statistics = TRUE,
                                      detail = qc_detail)
    qc_seconds <- phase_seconds(qc_started)
    if (identical(row_policy, "error") && qc_report$rejected_rows > 0L) {
      stop(format_structural_qc_failure(qc_report), call. = FALSE)
    }
  }
  attr(out, "source_columns") <- source_columns
  attr(out, "source_columns_read") <- source_columns_read
  if (!is.null(allele_columns)) {
    attr(out, "allele_columns") <- c(reference_allele = unname(allele_columns[["ref"]]),
                                     alternate_allele = unname(allele_columns[["alt"]]))
  }
  attr(out, "source_provenance") <- provenance
  attr(out, "eaf_provenance") <- provenance$eaf
  attr(out, "input_build") <- normalise_build_name(input_build)
  attr(out, "structural_qc_report") <- qc_report
  attr(out, "phase_timings") <- list(
    unit = "seconds",
    phases = c(
      list(
        read = read_seconds,
        projection = as.numeric(input_read_metadata$projection_elapsed_seconds %||% 0),
        statistic_resolution = normalise_seconds,
        normalise = normalise_seconds
      ),
      if (is.null(qc_seconds)) list() else list(qc = qc_seconds)
    )
  )
  out
}

#' Preflight arbitrary summary statistics
#'
#' Imports and canonicalises a GWAS, then applies the mandatory structural QC
#' layer. In report mode the
#' returned list contains rows that can proceed and a structured report;
#' duplicate groups retain their first canonical copy and later copies are
#' rejected. In error mode no rows are returned when any row is rejected.
#'
#' @param input A data.frame or a summary-statistics file accepted by
#'   import_sumstats().
#' @param input_build Build used for coordinate-bound validation.
#' @param strict Whether any rejected row should stop the operation.
#' @param row_policy Either report or error.
#' @param require_statistics Whether finite Z and positive standard error are
#'   required after canonical derivation.
#' @param max_examples Maximum row numbers retained per rejection reason.
#' @return A list with data and report components.
#' @noRd
preflight_sumstats <- function(input, input_build = "GRCh38", strict = FALSE,
                               row_policy = c("report", "error"),
                               require_statistics = TRUE, max_examples = 5L) {
  if (length(strict) != 1L || !is.logical(strict) || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  row_policy <- match.arg(row_policy)
  if (isTRUE(strict)) row_policy <- "error"
  imported <- import_sumstats_impl(
    input, strict = FALSE, row_policy = "report",
    input_build = input_build, run_qc = FALSE
  )
  result <- apply_structural_qc(imported, input_build = input_build,
                                row_policy = row_policy,
                                require_statistics = require_statistics,
                                max_examples = max_examples, detail = "full")
  attr(result$data, "source_columns") <- attr(imported, "source_columns")
  attr(result$data, "source_columns_read") <- attr(imported, "source_columns_read")
  attr(result$data, "source_provenance") <- attr(imported, "source_provenance")
  attr(result$data, "resolution_provenance") <- attr(imported, "resolution_provenance")
  attr(result$data, "input_build") <- attr(imported, "input_build")
  attr(result$data, "phase_timings") <- attr(imported, "phase_timings")
  result
}

# allele_columns = c(ref = "<column>", alt = "<column>") declares which input
# columns hold REF and ALT when the input has no explicit REF/ALT columns (for
# example a GWAS-SSF table, whose effect_allele is ALT and other_allele REF).
validate_allele_columns <- function(allele_columns) {
  if (is.null(allele_columns)) return(NULL)
  if (!is.character(allele_columns) || length(allele_columns) != 2L ||
      is.null(names(allele_columns)) ||
      !setequal(names(allele_columns), c("ref", "alt")) ||
      anyNA(allele_columns) || any(!nzchar(trimws(allele_columns))) ||
      identical(trimws(allele_columns[["ref"]]), trimws(allele_columns[["alt"]]))) {
    stop("allele_columns must be NULL or two distinct column names named ref and alt, ",
         "e.g. c(ref = \"other_allele\", alt = \"effect_allele\")", call. = FALSE)
  }
  out <- trimws(allele_columns[c("ref", "alt")])
  names(out) <- c("ref", "alt")
  out
}

# Add reference_allele/alternate_allele columns that refer to the declared
# input columns (no copy). An input that already has an explicit REF or ALT
# column cannot also declare one.
apply_allele_columns <- function(raw, allele_columns) {
  missing <- setdiff(unname(allele_columns), names(raw))
  if (length(missing)) {
    stop("allele_columns refers to column(s) not in the input: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  alias_map <- sumstats_column_alias_map()
  explicit <- names(raw)[alias_key(names(raw)) %in% alias_key(c(
    alias_map$reference_allele, alias_map$alternate_allele
  ))]
  if (length(explicit)) {
    stop("allele_columns cannot be combined with explicit REF/ALT column(s): ",
         paste(explicit, collapse = ", "), call. = FALSE)
  }
  ref <- raw[[allele_columns[["ref"]]]]
  alt <- raw[[allele_columns[["alt"]]]]
  raw$reference_allele <- ref
  raw$alternate_allele <- alt
  raw
}

# Provenance hash of an in-memory input table (manifest source$hash, recorded
# with source$hash_algorithm). Each column that was read is hashed by its
# content with XXH64 over a platform-independent byte stream (see
# compressor_column_xxh64() in src/write_native.cpp), the columns in parallel;
# the per-column digests are then combined with SHA-256 over a text record of
# the row count and each column's name, type and digest. The value does not
# depend on the thread count or platform. Stores written before this
# definition carry source$sha256 = object_sha256(input) instead; neither value
# is needed to open or validate a store.
COMPRESSOR_COLUMNS_HASH_ALGORITHM <- "xxh64_per_column_sha256_combined_v1"

columns_content_hash <- function(data, threads = 1L) {
  columns <- unclass(data)
  attributes(columns) <- NULL
  column_names <- names(data) %||% rep("", length(columns))
  types <- character(length(columns))
  digests <- character(length(columns))
  native <- logical(length(columns))
  for (j in seq_along(columns)) {
    x <- columns[[j]]
    class_tag <- oldClass(x)
    if (is.factor(x)) {
      # A factor is hashed by its labels.
      x <- as.character(x)
      columns[j] <- list(x)
    }
    types[[j]] <- paste0(typeof(x), if (length(class_tag)) {
      paste0("<", paste(class_tag, collapse = ","), ">")
    } else {
      ""
    })
    native[[j]] <- typeof(x) %in% c("logical", "integer", "double", "character")
    if (!native[[j]]) {
      types[[j]] <- paste0("serialized:", types[[j]])
      digests[[j]] <- object_sha256(x)
    }
  }
  if (any(native)) {
    digests[native] <- .Call("compressor_column_xxh64", columns[native],
                             as.integer(max(1L, threads)), PACKAGE = "CompreSSoR")
  }
  record <- paste0(
    "CompreSSoR column content hash v1\nrows\t", nrow(data), "\n",
    paste0(enc2utf8(column_names), "\t", types, "\t", digests, "\n", collapse = "")
  )
  .Call("compressor_sha256_raw", charToRaw(enc2utf8(record)), PACKAGE = "CompreSSoR")
}

# digest::digest(object, algo = "sha256", serialize = TRUE), computed by
# streaming the serialization through SHA-256 natively rather than first
# materialising the whole serialized object in memory. digest serializes with
# its configured version in binary XDR form and skips the 14-byte header; any
# other configuration falls back to digest itself. With threads > 1 a second
# thread hashes while the main thread serializes.
object_sha256 <- function(object, threads = 1L) {
  version <- tryCatch(
    utils::getFromNamespace(".getSerializeVersion", "digest")(),
    error = function(e) NULL
  )
  no_sharing <- tryCatch(
    isTRUE(utils::getFromNamespace(".hasNoSharing", "digest")()),
    error = function(e) TRUE
  )
  if (!no_sharing && length(version) == 1L && !is.na(version) &&
      as.integer(version) %in% c(2L, 3L) &&
      is.loaded("compressor_serialize_sha256", PACKAGE = "CompreSSoR")) {
    return(.Call("compressor_serialize_sha256", object, as.integer(version), 14,
                 as.integer(threads), PACKAGE = "CompreSSoR"))
  }
  digest::digest(object, algo = "sha256", serialize = TRUE)
}
