#' Harmonise summary statistics to a GRCh38 reference
#'
#' Normalises common summary-statistics column aliases, matches variants to a
#' reference spine by variant ID or chromosome and position, aligns alleles,
#' and flips beta and effect-allele frequency when required. Strand-complement
#' alignment is supported for unambiguous variants; palindromic ambiguity fails
#' closed unless informative effect-allele frequencies resolve it.
#'
#' @param input A data.frame or a delimited summary-statistics file.
#' @param reference A GRCh38 reference data.frame, list(variants = ...), a
#'   readable Parquet/TSV file, or "GRCh38" to use the cached/downloaded
#'   default reference.
#' @param mode One of `"qc"` (default, harmonise and preserve all rows),
#'   `"convert"` (minimal conversion without reference QC), `"all"` (alias
#'   for `"qc"`), `"core"`, or `"hm3"`.
#' @param variant_set Panel data.frame or file used by `mode = "core"` or
#'   `mode = "hm3"`. PLINK `.bim`, Parquet and delimited panel files are
#'   supported.
#' @param strict If `TRUE`, fail when variants are absent, incompatible,
#'   ambiguous, or duplicated. The default preserves such rows and records
#'   their status instead of dropping them.
#' @param chrom_threads Number of chromosome workers. Values above one use
#'   chromosome-parallel harmonisation and share the loaded reference.
#' @return A normalized data.frame with attributes `reference_hash`,
#'   `reference_rows`, and `genome_build`.
#' @export
harmonise_sumstats <- function(input, reference = "GRCh38",
                               strict = FALSE,
                               mode = c("qc", "convert", "all", "core", "hm3"),
                               variant_set = NULL,
                               chrom_threads = 1L) {
  raw <- normalise_sumstats_columns(read_sumstats_input(input))
  prepared <- prepare_sumstats_data(raw, reference, mode = mode, variant_set = variant_set,
                                    strict = strict, chrom_threads = chrom_threads)
  if (!identical(prepared$effective_mode, "convert")) validate_sumstats_values(prepared$data)
  out <- prepared$data
  attr(out, "reference_hash") <- prepared$alignment$reference_hash
  attr(out, "reference_rows") <- prepared$alignment$reference_rows
  attr(out, "alignment_stats") <- prepared$alignment$alignment_stats
  attr(out, "mode") <- prepared$requested_mode
  attr(out, "genome_build") <- prepared$genome_build
  out
}
