#' Import GWAS summary statistics
#'
#' Reads a data frame, delimited text file, gzip-compressed text file, or
#' single-ALT VCF and maps common GWAS column names onto CompreSSoR's canonical
#' in-memory schema. Unrecognised columns are retained. This function does not
#' change genome build or allele orientation; use [liftover_sumstats()] and
#' [harmonise_sumstats()] for those operations.
#'
#' @param input A data.frame or an existing summary-statistics file.
#' @param strict Whether structural QC failures should stop import. The
#'   default is FALSE, which keeps canonical rows and attaches a bounded
#'   structural_qc_report attribute for the caller to apply.
#' @param row_policy Either report (attach a report and continue) or error
#'   (stop on any structural rejection). strict = TRUE is an alias for error.
#' @param input_build Build used for coordinate-bound preflight, either
#'   GRCh37/hg19 or GRCh38/hg38.
#' @return A data.frame containing canonical columns such as `chromosome`,
#'   `base_pair_location`, `effect_allele`, `other_allele`, `beta`,
#'   `standard_error`, `z`, `effect_allele_frequency`, and `p_value`.
#' @export
import_sumstats <- function(input, strict = FALSE,
                            row_policy = c("report", "error"),
                            input_build = "GRCh38") {
  if (length(strict) != 1L || !is.logical(strict) || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  row_policy <- match.arg(row_policy)
  if (isTRUE(strict)) row_policy <- "error"
  raw <- read_sumstats_input(
    input, parse_policy = if (identical(row_policy, "error")) "error" else "report"
  )
  source_columns <- attr(raw, "source_columns") %||% names(raw)
  explicit_ref_alt <- isTRUE(attr(raw, "explicit_ref_alt"))
  provenance <- if (is.data.frame(input)) {
    list(
      kind = "data.frame",
      rows = nrow(input),
      sha256 = digest::digest(input, algo = "sha256", serialize = TRUE)
    )
  } else {
    path <- normalizePath(input, mustWork = TRUE)
    list(
      kind = "file",
      file = basename(path),
      bytes = unname(file.info(path)$size),
      sha256 = digest::digest(path, algo = "sha256", file = TRUE)
    )
  }
  out <- normalise_sumstats_columns(
    raw, parse_policy = if (identical(row_policy, "error")) "error" else "report"
  )
  qc_report <- structural_qc_report(out, input_build = input_build,
                                     require_statistics = TRUE)
  if (identical(row_policy, "error") && length(qc_report$invalid_rows)) {
    stop(format_structural_qc_failure(qc_report), call. = FALSE)
  }
  attr(out, "source_columns") <- source_columns
  attr(out, "source_provenance") <- provenance
  attr(out, "explicit_ref_alt") <- explicit_ref_alt
  attr(out, "input_build") <- normalise_build_name(input_build)
  attr(out, "structural_qc_report") <- qc_report
  out
}

#' Preflight arbitrary summary statistics
#'
#' Imports and canonicalises a GWAS, then applies the mandatory structural QC
#' layer before optional reference-backed harmonisation. In report mode the
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
  imported <- import_sumstats(input, strict = FALSE, row_policy = "report",
                              input_build = input_build)
  result <- apply_structural_qc(imported, input_build = input_build,
                                row_policy = row_policy,
                                require_statistics = require_statistics,
                                max_examples = max_examples)
  attr(result$data, "source_columns") <- attr(imported, "source_columns")
  attr(result$data, "source_provenance") <- attr(imported, "source_provenance")
  attr(result$data, "explicit_ref_alt") <- attr(imported, "explicit_ref_alt")
  attr(result$data, "input_build") <- attr(imported, "input_build")
  result
}

#' Lift GWAS summary statistics to GRCh38
#'
#' Uses a caller-supplied UCSC chain through Bioconductor's `rtracklayer` and
#' preserves an explicit row-level mapping status. Reverse-strand mappings
#' reverse-complement both alleles. Zero- and multi-mapped rows are dropped by
#' default and always counted in the `liftover_stats` attribute.
#'
#' @param input A data.frame or summary-statistics file accepted by
#'   [import_sumstats()].
#' @param input_build Source genome build, normally `"GRCh37"`. GRCh38 inputs
#'   are returned with status `"not_needed"`.
#' @param chain A local GRCh37-to-GRCh38 chain file. Required unless
#'   `input_build` already denotes GRCh38.
#' @param drop_unmapped Drop rows that do not map exactly once.
#' @return A canonical data.frame with a `liftover_status` column and
#'   `liftover_stats` and `genome_build` attributes.
#' @export
liftover_sumstats <- function(input, input_build = "GRCh37", chain = NULL,
                              drop_unmapped = TRUE) {
  if (length(drop_unmapped) != 1L || !is.logical(drop_unmapped) || is.na(drop_unmapped)) {
    stop("drop_unmapped must be TRUE or FALSE", call. = FALSE)
  }
  data <- import_sumstats(input)
  lifted <- lift_table_to_grch38(data, input_build = input_build, chain = chain)
  status <- lifted$.compressor_liftover_status %||% rep("not_needed", nrow(lifted))
  stats <- as.list(table(status, useNA = "ifany"))
  if (isTRUE(drop_unmapped)) {
    keep <- status %in% c("mapped", "not_needed")
    lifted <- lifted[keep, , drop = FALSE]
    status <- status[keep]
  }
  lifted$liftover_status <- status
  if (".compressor_source_variant_id" %in% names(lifted)) {
    lifted$source_variant_id <- lifted$.compressor_source_variant_id
    lifted$.compressor_source_variant_id <- NULL
  }
  lifted$.compressor_liftover_status <- NULL
  lifted$.compressor_liftover_row <- NULL
  row.names(lifted) <- NULL
  attr(lifted, "liftover_stats") <- stats
  attr(lifted, "genome_build") <- "GRCh38"
  lifted
}
