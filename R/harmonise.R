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
#' @param input_build Build of the input summary statistics.
#' @param target_build Build produced by harmonisation. Reference-backed
#'   harmonisation currently targets GRCh38; already aligned GRCh37 input may
#'   be passed through with `reference = NULL`.
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
#' @param qc Optional list controlling reference QC. Supported fields are
#'   `liftover`, `strand`, `palindromic` (`"drop"`, `"frequency"`, or
#'   `"allow"`), `frequency` (`"none"`, `"report"`, `"drop"`, or
#'   `"error"`), `frequency_tolerance`, and `example_limit`.
#' @param liftover Optional scalar override for `qc$liftover`.
#' @param strand Optional scalar override for `qc$strand`.
#' @param palindromic Optional override for `qc$palindromic`.
#' @param frequency_qc Optional override for `qc$frequency`.
#' @param frequency_tolerance Optional override for `qc$frequency_tolerance`.
#' @param max_examples Optional override for `qc$example_limit`.
#' @return A normalized data.frame with attributes `reference_hash`,
#'   `reference_rows`, `genome_build`, `source_provenance`, and bounded audit
#'   diagnostics.
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
                               region_padding = 50000L,
                               target_build = "GRCh38", qc = NULL,
                               liftover = NULL, strand = NULL, palindromic = NULL,
                               frequency_qc = NULL, frequency_tolerance = NULL,
                               max_examples = NULL) {
  input_build <- compressor_normalize_build(input_build)
  target_build <- compressor_normalize_build(target_build)
  if (!is.null(reference) && !identical(target_build, "GRCh38")) {
    stop("reference-backed harmonisation currently targets GRCh38", call. = FALSE)
  }
  raw <- import_sumstats(input)
  source_provenance <- attr(raw, "source_provenance")
  source_columns <- attr(raw, "source_columns")
  imported_explicit_ref_alt <- isTRUE(attr(raw, "explicit_ref_alt"))
  structural <- apply_structural_qc(
    raw, input_build = input_build, strict = strict,
    # Harmonisation has historically failed closed for malformed statistics or
    # alleles; only duplicate handling is deferred to reference alignment so
    # the reference layer can report its canonical-target diagnostics.
    row_policy = "error",
    require_statistics = TRUE, drop_duplicates = FALSE,
    check_duplicates = is.null(reference) && isTRUE(strict)
  )
  raw <- structural$data
  attr(raw, "source_provenance") <- source_provenance
  attr(raw, "source_columns") <- source_columns
  attr(raw, "explicit_ref_alt") <- imported_explicit_ref_alt
  source_keys <- alias_key(attr(raw, "source_columns") %||% character())
  explicit_ref_alt <- imported_explicit_ref_alt ||
    (any(source_keys %in% c("ref", "referenceallele")) &&
     any(source_keys %in% c("alt", "alternateallele")))
  qc_input <- if (is.null(qc)) list() else qc
  if (is.logical(qc_input) && length(qc_input) == 1L && !is.na(qc_input)) {
    qc_input <- list(strand = qc_input)
  }
  if (!is.list(qc_input)) {
    stop("qc must be NULL, logical, or a named list", call. = FALSE)
  }
  if (!is.null(liftover)) qc_input$liftover <- liftover
  if (!is.null(strand)) qc_input$strand <- strand
  if (!is.null(palindromic)) qc_input$palindromic <- palindromic
  if (!is.null(frequency_qc)) qc_input$frequency <- frequency_qc
  if (!is.null(frequency_tolerance)) qc_input$frequency_tolerance <- frequency_tolerance
  if (!is.null(max_examples)) qc_input$example_limit <- max_examples
  qc_config <- normalise_reference_qc(qc_input)
  had_qc <- exists("qc", envir = .compressor_harmonise_context, inherits = FALSE)
  old_qc <- if (had_qc) .compressor_harmonise_context$qc else NULL
  .compressor_harmonise_context$qc <- qc_config
  on.exit({
    if (had_qc) .compressor_harmonise_context$qc <- old_qc else
      rm("qc", envir = .compressor_harmonise_context)
  }, add = TRUE)
  if (!isTRUE(qc_config$liftover) && !identical(input_build, target_build)) {
    stop("qc$liftover = FALSE requires input_build and target_build to match",
         call. = FALSE)
  }
  effective_input_build <- if (isTRUE(qc_config$liftover)) input_build else target_build
  prepared <- prepare_sumstats_data(raw, reference, mode = mode, variant_set = variant_set,
                                    strict = strict, chrom_threads = chrom_threads,
                                    drop_unresolved = drop_unresolved,
                                    input_build = effective_input_build, chain = chain,
                                    pvalue_threshold = pvalue_threshold,
                                    region_padding = region_padding,
                                    target_build = target_build)
  validate_sumstats_values(prepared$data, require_identity = isTRUE(drop_unresolved))
  out <- prepared$data
  out$.compressor_reference_index <- NULL
  alignment_stats <- prepared$alignment$alignment_stats
  alignment_stats$structural_qc <- compact_structural_qc_report(structural$report)
  alignment_stats$source_provenance <- source_provenance
  alignment_stats$input_build <- input_build
  alignment_stats$target_build <- prepared$genome_build
  attr(out, "reference_hash") <- prepared$alignment$reference_hash
  attr(out, "reference_rows") <- prepared$alignment$reference_rows
  attr(out, "alignment_stats") <- alignment_stats
  attr(out, "mode") <- prepared$requested_mode
  attr(out, "selection") <- prepared$selection
  attr(out, "genome_build") <- prepared$genome_build
  attr(out, "reference_metadata") <- prepared$alignment$reference_metadata
  attr(out, "explicit_ref_alt") <- explicit_ref_alt
  attr(out, "source_columns") <- source_columns
  attr(out, "source_provenance") <- source_provenance
  attr(out, "harmonisation_qc") <- qc_config
  attr(out, "diagnostics") <- alignment_stats$diagnostics %||% list(
    counts = alignment_stats$counts %||% list(),
    examples = alignment_stats$examples %||% list()
  )
  attr(out, "audit") <- list(
    source = source_provenance,
    reference = prepared$alignment$reference_metadata,
    counts = alignment_stats$diagnostic_counts %||% alignment_stats,
    examples = alignment_stats$diagnostic_examples %||% list(),
    qc = qc_config
  )
  if (!is.null(alignment_stats$liftover)) attr(out, "liftover_stats") <- alignment_stats$liftover
  attr(out, "compressor_identity_verified") <- explicit_ref_alt ||
    (!is.null(reference) && !identical(prepared$effective_mode, "convert"))
  out
}
