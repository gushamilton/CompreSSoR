## Strict prepared-input contract for the installed package.

core_builds <- function(input_build, store_build) {
  if (is.null(input_build) || is.null(store_build)) {
    stop("input_build and store_build are required; choose the same GRCh37/hg19 " ,
         "or GRCh38/hg38 build explicitly", call. = FALSE)
  }
  input_build <- compressor_normalize_build(input_build)
  store_build <- compressor_normalize_build(store_build)
  if (!identical(input_build, store_build)) {
    stop("input_build and store_build must be the same build for the compression " ,
         "core; build conversion is outside CompreSSoR (received ", input_build, " -> ",
         store_build, ")", call. = FALSE)
  }
  list(input_build = input_build, store_build = store_build)
}

core_source_has_alias <- function(source_columns, aliases) {
  any(alias_key(source_columns %||% character()) %in% alias_key(aliases))
}

validate_core_schema <- function(data, source_columns = names(data),
                                 input_build, store_build) {
  core_builds(input_build, store_build)
  source_columns <- source_columns %||% names(data)
  if (!core_source_has_alias(source_columns,
                             c("reference_allele", "reference", "ref", "REF")) ||
      !core_source_has_alias(source_columns,
                             c("alternate_allele", "alternate", "alt", "ALT"))) {
    stop("strict compression core requires explicit REF and ALT columns " ,
         "(accepted aliases include reference_allele/alternate_allele or ref/alt); " ,
         "effect_allele and other_allele do not establish REF/ALT orientation", call. = FALSE)
  }
  if (!core_source_has_alias(source_columns,
                             c("effect_allele", "ea", "EA", "a1", "A1",
                               "ALLELE1", "alleleB", "ALLELEB")) ||
      !core_source_has_alias(source_columns,
                             c("other_allele", "oa", "NEA", "nea", "a2", "A2",
                               "ALLELE0", "alleleA", "ALLELEA"))) {
    stop("strict compression core requires explicit effect_allele and other_allele columns",
         call. = FALSE)
  }
  missing <- setdiff(required_sumstats_columns(), names(data))
  if (length(missing)) {
    stop("strict compression core is missing canonical resolved fields: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

validate_core_orientation <- function(data, row_policy = c("error", "report")) {
  row_policy <- match.arg(row_policy)
  ref <- toupper(trimws(as.character(data$reference_allele)))
  alt <- toupper(trimws(as.character(data$alternate_allele)))
  effect <- toupper(trimws(as.character(data$effect_allele)))
  other <- toupper(trimws(as.character(data$other_allele)))
  bases <- c("A", "C", "G", "T")
  invalid_identity <- is.na(ref) | is.na(alt) | !ref %in% bases | !alt %in% bases | ref == alt
  if (any(invalid_identity) && identical(row_policy, "error")) {
    stop("explicit REF and ALT must be distinct single A/C/G/T alleles; " ,
         sum(invalid_identity), " invalid row(s) found", call. = FALSE)
  }
  inconsistent <- is.na(effect) | is.na(other) | effect != alt | other != ref
  if (any(inconsistent) && identical(row_policy, "error")) {
    stop("effect_allele and other_allele are inconsistent with explicit REF/ALT; " ,
         "the core requires effect_allele = ALT and other_allele = REF, with no silent flipping (" ,
         sum(inconsistent), " row(s))", call. = FALSE)
  }
  data$reference_allele <- ref
  data$alternate_allele <- alt
  data$effect_allele <- effect
  data$other_allele <- other
  data
}

canonicalize_core_identity <- function(data, build) {
  data$variant_id <- compressor_variant_key(
    data$chromosome, data$base_pair_location,
    data$reference_allele, data$alternate_allele, build = build
  )
  attr(data, "genome_build") <- build
  attr(data, "compressor_identity_verified") <- TRUE
  data
}

prepare_core_sumstats_data <- function(raw, selection = "full", variant_set = NULL,
                                       pvalue_threshold = 1e-5,
                                       region_padding = 50000L, build = "GRCh38") {
  build <- compressor_normalize_build(build)
  selected <- select_variant_rows(
    raw, selection = selection, variant_set = variant_set,
    pvalue_threshold = pvalue_threshold, region_padding = region_padding,
    build = build
  )
  data <- selected$data
  selection_metadata <- selected$metadata %||% list()
  selection_result <- NULL
  if (selection %in% c("pvalue_regions", "core_plus")) {
    selection_result <- selection_metadata
    selection_result$regions_table <- selected$regions
  }
  preparation_stats <- list(
    method = "strict_prepared_input",
    input_rows = as.integer(nrow(raw)),
    output_rows = as.integer(nrow(data)),
    selected_rows = as.integer(nrow(data)),
    reference_lookup = "not_used",
    build_conversion = "not_used",
    rsid_resolution = "not_used",
    variant_set = selection_metadata
  )
  list(
    data = data,
    preparation = preparation_stats,
    selection = selection_result, selection_result = selected,
    genome_build = build
  )
}
