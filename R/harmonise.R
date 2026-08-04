#' Harmonise summary statistics to a GRCh38 reference
#'
#' Normalises common summary-statistics column aliases, matches variants to the
#' canonical `chrom:pos:ref:alt` key, aligns alleles, and flips beta/Z and
#' effect-allele frequency when required. Unresolved, duplicate, incompatible
#' and unmatched rows are dropped by default.
#'
#' @param input A data.frame or a delimited summary-statistics file.
#' @param reference A GRCh38 reference data.frame, list(variants = ...), a
#'   readable Parquet/TSV file, or `"GRCh38"` to use the immutable reference
#'   configured by `COMPRESSOR_CANONICAL_REFERENCE`.
#' @param mode One of `"qc"` (default, harmonise and drop unresolved rows),
#'   `"convert"` (minimal conversion without reference QC), `"all"` (alias
#'   for `"qc"`), `"core"`, or `"hm3"`.
#' @param variant_set Panel data.frame or file used by `mode = "core"` or
#'   `mode = "hm3"`. PLINK `.bim`, Parquet and delimited panel files are
#'   supported.
#' @param strict If `TRUE`, fail when variants are absent, incompatible,
#'   ambiguous, or duplicated.
#' @param drop_unresolved Whether unresolved rows are dropped; defaults to `TRUE`.
#' @param input_build Build of the incoming summary statistics.
#' @param chain Optional GRCh37-to-GRCh38 chain file for non-GRCh38 input.
#' @param chrom_threads Number of chromosome workers. Values above one use
#'   chromosome-parallel harmonisation and share the loaded reference.
#' @return A normalized data.frame with attributes `reference_hash`,
#'   `reference_rows`, and `genome_build`.
#' @examples
#' study <- data.frame(
#'   chromosome = "1", position = 100L, effect_allele = "A",
#'   other_allele = "C", beta = 0.2, standard_error = 0.1, eaf = 0.3
#' )
#' reference <- data.frame(
#'   chromosome = "1", base_pair_location = 100L,
#'   reference_allele = "A", alternate_allele = "C"
#' )
#' harmonise_sumstats(study, reference)
#' @export
harmonise_sumstats <- function(input, reference = "GRCh38",
                               strict = FALSE,
                               mode = c("qc", "convert", "all", "core", "hm3"),
                               variant_set = NULL,
                               chrom_threads = 1L,
                               drop_unresolved = TRUE,
                               input_build = "GRCh38", chain = NULL) {
  raw <- import_sumstats(input)
  prepared <- prepare_sumstats_data(raw, reference, mode = mode, variant_set = variant_set,
                                    strict = strict, chrom_threads = chrom_threads,
                                    drop_unresolved = drop_unresolved,
                                    input_build = input_build, chain = chain)
  validate_sumstats_values(prepared$data, require_identity = isTRUE(drop_unresolved))
  out <- prepared$data
  out$.compressor_reference_index <- NULL
  attr(out, "reference_hash") <- prepared$alignment$reference_hash
  attr(out, "reference_rows") <- prepared$alignment$reference_rows
  attr(out, "alignment_stats") <- prepared$alignment$alignment_stats
  attr(out, "mode") <- prepared$requested_mode
  attr(out, "genome_build") <- prepared$genome_build
  out
}
