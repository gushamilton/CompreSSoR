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
#'   for `"qc"`), `"core"`, `"hm3"`, `"pvalue_regions"`, or `"core_plus"`.
#' @param variant_set Optional panel data.frame, file, or chromosome-shard
#'   directory used by `mode = "core"`, `mode = "hm3"`, or `mode = "core_plus"`.
#'   Named panels resolve through the configured environment variables. Staged
#'   `*_by_chrom` directories are filtered to chromosomes present after
#'   harmonisation.
#' @param strict If `TRUE`, fail when variants are absent, incompatible,
#'   ambiguous, or duplicated.
#' @param drop_unresolved Whether unresolved rows are dropped; defaults to `TRUE`.
#' @param input_build Build of the incoming summary statistics.
#' @param chain Optional GRCh37-to-GRCh38 chain file for non-GRCh38 input.
#' @param pvalue_threshold Strict p-value threshold for `mode = "pvalue_regions"`
#'   or `mode = "core_plus"`. P-values are derived from the canonical Z score.
#' @param region_padding Number of base pairs added on each side of each
#'   significant SNP for `mode = "pvalue_regions"` or `mode = "core_plus"`.
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
                               mode = c("qc", "convert", "all", "core", "hm3", "pvalue_regions", "core_plus"),
                               variant_set = NULL,
                               chrom_threads = 1L,
                               drop_unresolved = TRUE,
                               input_build = "GRCh38", chain = NULL,
                               pvalue_threshold = 1e-5,
                               region_padding = 50000L) {
  raw <- import_sumstats(input)
  source_keys <- alias_key(attr(raw, "source_columns") %||% character())
  explicit_ref_alt <- isTRUE(attr(raw, "explicit_ref_alt")) ||
    (any(source_keys %in% c("ref", "referenceallele")) &&
     any(source_keys %in% c("alt", "alternateallele")))
  prepared <- prepare_sumstats_data(raw, reference, mode = mode, variant_set = variant_set,
                                    strict = strict, chrom_threads = chrom_threads,
                                    drop_unresolved = drop_unresolved,
                                    input_build = input_build, chain = chain,
                                    pvalue_threshold = pvalue_threshold,
                                    region_padding = region_padding)
  validate_sumstats_values(prepared$data, require_identity = isTRUE(drop_unresolved))
  out <- prepared$data
  out$.compressor_reference_index <- NULL
  attr(out, "reference_hash") <- prepared$alignment$reference_hash
  attr(out, "reference_rows") <- prepared$alignment$reference_rows
  attr(out, "alignment_stats") <- prepared$alignment$alignment_stats
  attr(out, "mode") <- prepared$requested_mode
  attr(out, "selection") <- prepared$selection
  attr(out, "genome_build") <- prepared$genome_build
  attr(out, "reference_metadata") <- prepared$alignment$reference_metadata
  attr(out, "explicit_ref_alt") <- explicit_ref_alt
  attr(out, "compressor_identity_verified") <- explicit_ref_alt ||
    (!is.null(reference) && !identical(prepared$effective_mode, "convert"))
  out
}
