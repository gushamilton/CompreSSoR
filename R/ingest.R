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
                                 qc_detail = c("full", "compact")) {
  if (length(strict) != 1L || !is.logical(strict) || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  row_policy <- match.arg(row_policy)
  qc_detail <- match.arg(qc_detail)
  if (isTRUE(strict)) row_policy <- "error"
  read_started <- phase_clock()
  raw <- read_sumstats_input(
    input,
    parse_policy = if (identical(row_policy, "error")) "error" else "report",
    project_columns = project_columns,
    core_only = core_only,
    allow_p_to_se = allow_p_to_se
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
      sha256 = digest::digest(input, algo = "sha256", serialize = TRUE),
      columns_before = as.integer(length(source_columns)),
      columns_read = as.integer(length(source_columns_read)),
      projected = isTRUE(input_read_metadata$projected),
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
      read_elapsed_seconds = as.numeric(input_read_metadata$elapsed_seconds %||% NA_real_)
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
  attr(out, "source_provenance") <- provenance
  attr(out, "input_build") <- normalise_build_name(input_build)
  attr(out, "structural_qc_report") <- qc_report
  attr(out, "phase_timings") <- list(
    unit = "seconds",
    phases = c(
      list(read = read_seconds, normalise = normalise_seconds),
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
