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

# The qc='none' path is intentionally not a structural-QC path.  It still
# needs a cheap guard before canonicalize_core_identity(), because the native
# identity encoder cannot represent rows outside its primary-chromosome,
# whole-coordinate, biallelic A/C/G/T domain.
filter_pcodec_identity_safety <- function(data, build = "GRCh38") {
  build <- compressor_normalize_build(build)
  required <- required_sumstats_columns()
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("qc='none' identity safety needs canonical columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  n <- nrow(data)
  chromosome <- as.character(data$chromosome)
  position <- suppressWarnings(as.numeric(as.character(data$base_pair_location)))
  reference <- toupper(trimws(as.character(data$reference_allele)))
  alternate <- toupper(trimws(as.character(data$alternate_allele)))
  lengths <- compressor_chromosome_lengths(build)
  known_chromosome <- !is.na(chromosome) & chromosome %in% names(lengths)
  finite_position <- is.finite(position)
  integer_position <- finite_position & position == floor(position)
  chromosome_limit <- rep(NA_real_, n)
  chromosome_limit[known_chromosome] <- unname(lengths[chromosome[known_chromosome]])

  invalid_primary_chromosome <- !known_chromosome
  nonfinite_coordinate <- known_chromosome & !finite_position
  noninteger_coordinate <- known_chromosome & finite_position & !integer_position
  coordinate_out_of_range <- known_chromosome & finite_position & integer_position &
    (position < 1 | position > chromosome_limit)
  valid_alleles <- !is.na(reference) & !is.na(alternate) &
    reference %in% c("A", "C", "G", "T") &
    alternate %in% c("A", "C", "G", "T")
  invalid_allele <- !valid_alleles
  same_alleles <- valid_alleles & reference == alternate
  unsupported <- invalid_primary_chromosome | nonfinite_coordinate |
    noninteger_coordinate | coordinate_out_of_range | invalid_allele | same_alleles
  unsupported[is.na(unsupported)] <- TRUE

  counts <- c(
    invalid_primary_chromosome = sum(invalid_primary_chromosome),
    nonfinite_coordinate = sum(nonfinite_coordinate),
    noninteger_coordinate = sum(noninteger_coordinate),
    coordinate_out_of_range = sum(coordinate_out_of_range),
    invalid_allele = sum(invalid_allele),
    same_alleles = sum(same_alleles)
  )
  report <- list(
    mode = "pre_canonicalization",
    input_rows = as.integer(n),
    kept_rows = as.integer(sum(!unsupported)),
    dropped_rows = as.integer(sum(unsupported)),
    dropped_unsupported_identity = as.integer(sum(unsupported)),
    counts = stats::setNames(as.list(as.integer(counts)), names(counts)),
    rejection_counts = stats::setNames(as.list(as.integer(counts)), names(counts))
  )
  list(data = data[!unsupported, , drop = FALSE],
       keep = !unsupported, report = report)
}

# This is the statistic-safety boundary for qc='none'.  It validates only the
# canonical numerical core and does not run structural QC, row reporting, or
# duplicate detection.  Missing EAF/P are representable; required beta/SE and
# supplied non-missing numerical values remain fail-closed.
validate_qc_none_statistics <- function(data) {
  numeric_value <- function(name) {
    suppressWarnings(as.numeric(as.character(data[[name]])))
  }
  beta <- numeric_value("beta")
  if (any(!is.finite(beta))) {
    stop("qc='none' requires finite beta values", call. = FALSE)
  }
  standard_error <- numeric_value("standard_error")
  if (any(!is.finite(standard_error) | standard_error <= 0)) {
    stop("qc='none' requires positive finite standard_error values", call. = FALSE)
  }
  for (field in c("z", "effect_allele_frequency", "p_value")) {
    if (!field %in% names(data)) next
    value <- numeric_value(field)
    if (field == "z" && any(!is.na(value) & !is.finite(value))) {
      stop("qc='none' requires finite z values when supplied", call. = FALSE)
    }
    if (field == "z") {
      supplied <- attr(data, "z_supplied")
      if (is.null(supplied) || length(supplied) != length(value)) {
        supplied <- !is.na(value)
      }
      supplied <- !is.na(supplied) & as.logical(supplied)
      expected <- beta / standard_error
      comparable <- supplied & is.finite(value) & is.finite(expected)
      tolerance <- 1e-6 + 0.01 * pmax(abs(value), abs(expected))
      conflict <- comparable & abs(value - expected) > tolerance
      if (any(conflict)) {
        stop(
          "qc='none' supplied z is inconsistent with beta / standard_error at row(s): ",
          paste(utils::head(which(conflict), 5L), collapse = ", "), call. = FALSE
        )
      }
    }
    if (field == "effect_allele_frequency" && any(
      !is.na(value) & (!is.finite(value) | value < 0 | value > 1)
    )) {
      stop("qc='none' requires effect_allele_frequency values in [0, 1]", call. = FALSE)
    }
    if (field == "p_value" && any(
      !is.na(value) & (!is.finite(value) | value < 0 | value > 1)
    )) {
      stop("qc='none' requires p_value values in [0, 1]", call. = FALSE)
    }
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

canonicalize_core_identity <- function(data, build, include_variant_id = TRUE) {
  build <- compressor_normalize_build(build)
  identity <- compressor_encode_variant_identity(
    data$chromosome, data$base_pair_location,
    data$reference_allele, data$alternate_allele, build = build
  )
  data$chromosome <- identity$chromosome
  data$base_pair_location <- identity$position
  data$reference_allele <- identity$reference_allele
  data$alternate_allele <- identity$alternate_allele
  data$effect_allele <- identity$alternate_allele
  data$other_allele <- identity$reference_allele
  if (isTRUE(include_variant_id)) {
    data$variant_id <- compressor_variant_key(
      identity$chromosome, identity$position,
      identity$reference_allele, identity$alternate_allele, build = build
    )
  }
  attr(data, "compressor_identity") <- list(
    global_position = identity$global_position,
    substitution = identity$substitution,
    code = compressor_identity_code(identity$global_position, identity$substitution)
  )
  attr(data, "genome_build") <- build
  attr(data, "compressor_identity_verified") <- TRUE
  data
}

prepare_core_sumstats_data <- function(raw, selection = "full", variant_set = NULL,
                                       pvalue_threshold = 1e-5,
                                       region_padding = 50000L, build = "GRCh38",
                                       input_rows = NULL, identity_safety = NULL) {
  build <- compressor_normalize_build(build)
  input_rows <- input_rows %||% nrow(raw)
  if (length(input_rows) != 1L || is.na(input_rows) ||
      input_rows < nrow(raw) || input_rows != floor(input_rows)) {
    stop("prepared input row count is inconsistent with identity filtering", call. = FALSE)
  }
  selected <- select_variant_rows(
    raw, selection = selection, variant_set = variant_set,
    pvalue_threshold = pvalue_threshold, region_padding = region_padding,
    build = build
  )
  data <- selected$data
  identity <- attr(raw, "compressor_identity")
  if (is.list(identity) && length(identity$code) == nrow(raw) &&
      length(selected$keep) == nrow(raw)) {
    identity <- lapply(identity, function(values) values[selected$keep])
    attr(data, "compressor_identity") <- identity
  }
  selection_metadata <- selected$metadata %||% list()
  selection_result <- NULL
  if (selection %in% c("pvalue_regions", "core_plus")) {
    selection_result <- selection_metadata
    selection_result$regions_table <- selected$regions
  }
  if (!is.null(identity_safety)) {
    selection_metadata$source_input_rows <- as.integer(input_rows)
    selection_metadata$identity_safety_dropped_rows <-
      as.integer(identity_safety$dropped_rows %||% 0L)
    if (!is.null(selection_result)) {
      selection_result$source_input_rows <- as.integer(input_rows)
      selection_result$identity_safety_dropped_rows <-
        as.integer(identity_safety$dropped_rows %||% 0L)
    }
  }
  preparation_stats <- list(
    method = "strict_prepared_input",
    input_rows = as.integer(input_rows),
    rows_after_identity_safety = as.integer(nrow(raw)),
    output_rows = as.integer(nrow(data)),
    selected_rows = as.integer(nrow(data)),
    reference_lookup = "not_used",
    build_conversion = "not_used",
    rsid_resolution = "not_used",
    variant_set = selection_metadata,
    identity_safety = identity_safety
  )
  list(
    data = data,
    preparation = preparation_stats,
    selection = selection_result, selection_metadata = selection_metadata,
    selection_result = selected,
    genome_build = build
  )
}
