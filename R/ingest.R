#' Import GWAS summary statistics
#'
#' Reads a data frame, delimited text file, gzip-compressed text file, or
#' single-ALT VCF and maps common GWAS column names onto CompreSSoR's canonical
#' in-memory schema. Unrecognised columns are retained. This function does not
#' change genome build or allele orientation; use [liftover_sumstats()] and
#' [harmonise_sumstats()] for those operations.
#'
#' @param input A data.frame or an existing summary-statistics file.
#' @return A data.frame containing canonical columns such as `chromosome`,
#'   `base_pair_location`, `effect_allele`, `other_allele`, `beta`,
#'   `standard_error`, `z`, `effect_allele_frequency`, and `p_value`.
#' @export
import_sumstats <- function(input) {
  raw <- read_sumstats_input(input)
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
  out <- normalise_sumstats_columns(raw)
  attr(out, "source_columns") <- source_columns
  attr(out, "source_provenance") <- provenance
  attr(out, "explicit_ref_alt") <- explicit_ref_alt
  out
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
