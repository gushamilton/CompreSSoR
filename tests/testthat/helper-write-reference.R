# Reference implementations from CompreSSoR a27b32d (before the write-path
# speedups), used to check that the faster versions give identical results.

reference_structural_qc_report <- function(data, input_build = "GRCh38",
                                 require_statistics = TRUE, max_examples = 5L,
                                 detail = c("full", "compact")) {
  detail <- match.arg(detail)
  data <- clean_input_names(data)
  n <- nrow(data)
  report <- new_structural_qc_report(
    n, input_build, max_examples = max_examples, detail = detail
  )
  report$missing_columns <- attr(data, "missing_columns") %||% character()
  if (!n) {
    report$accepted_rows <- 0L
    report$rejected_rows <- 0L
    report$dropped_rows <- 0L
    report$valid <- TRUE
    if (identical(detail, "full")) {
      report$row_status$structurally_valid <- logical()
      report$row_status$accepted <- logical()
    }
    return(report)
  }

  reasons <- if (identical(detail, "full")) list() else NULL
  rejection_counts <- integer()
  examples <- list()
  invalid <- rep(FALSE, n)
  reason_count <- integer(n)
  add_reason <- function(name, mask) {
    mask <- as.logical(mask)
    if (anyNA(mask)) mask[is.na(mask)] <- FALSE
    if (length(mask) != n) stop("structural QC mask has the wrong number of rows",
                                call. = FALSE)
    if (identical(detail, "full")) reasons[[name]] <<- mask
    count <- sum(mask)
    rejection_counts[name] <<- as.integer(count)
    examples[[name]] <<- utils::head(which(mask), max(0L, as.integer(max_examples)))
    # An all-FALSE mask leaves the running invalid/reason_count state
    # unchanged, so skip the two full-length passes for reasons that no row
    # triggers (the common case).
    if (count) {
      invalid <<- invalid | mask
      reason_count <<- reason_count + as.integer(mask)
    }
  }
  values <- function(name, default = NA_real_) {
    if (name %in% names(data)) return(data[[name]])
    rep(default, n)
  }
  # Elementwise text cleaning, evaluated once per distinct value.
  text_values <- function(name, upper = FALSE) {
    map_unique_values(values(name, NA_character_), function(value) {
      out <- trimws(as.character(value))
      out[is.na(value) | out %in% c("", ".", "NA", "N/A")] <- NA_character_
      if (upper) toupper(out) else out
    })
  }

  chromosome <- normalise_chromosome(values("chromosome", NA_character_))
  position <- suppressWarnings(as.numeric(values("base_pair_location", NA_real_)))
  reference <- text_values("reference_allele", upper = TRUE)
  alternate <- text_values("alternate_allele", upper = TRUE)
  effect <- text_values("effect_allele", upper = TRUE)
  other <- text_values("other_allele", upper = TRUE)
  rsid <- text_values("rsid")
  variant_id <- text_values("variant_id")
  alias_available <- (!is.na(rsid) & nzchar(rsid)) |
    map_unique_values(variant_id, function(value) {
      grepl("^rs", value, ignore.case = TRUE) |
        grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$",
              value, ignore.case = TRUE)
    })
  lengths <- sumstats_chromosome_lengths(input_build)
  primary <- sumstats_primary_chromosomes()
  # lengths[chromosome] without materialising a names attribute per row.
  chromosome_length <- unname(lengths)[match(chromosome, names(lengths))]

  add_reason("missing_chromosome",
             (is.na(chromosome) | !nzchar(chromosome)) & !alias_available)
  add_reason("unsupported_contig", !is.na(chromosome) & nzchar(chromosome) &
               !(chromosome %in% primary))
  add_reason("missing_coordinate", is.na(position) & !alias_available)
  parse_failures <- attr(data, "parse_failures") %||% list()
  add_reason("malformed_coordinate",
             is.na(position) &
               row_index_mask(n, parse_failures$base_pair_location %||% integer()) &
               !alias_available)
  add_reason("nonpositive_coordinate", !is.na(position) & is.finite(position) & position < 1)
  add_reason("noninteger_coordinate", !is.na(position) & is.finite(position) &
               position != trunc(position))
  known_chromosome <- !is.na(chromosome) & chromosome %in% names(lengths)
  add_reason("coordinate_out_of_range", known_chromosome & is.finite(position) &
               position >= 1 & position > chromosome_length)

  add_reason("missing_reference_allele", is.na(reference) | !nzchar(reference))
  add_reason("missing_alternate_allele", is.na(alternate) | !nzchar(alternate))
  add_reason("missing_effect_allele", is.na(effect) | !nzchar(effect))
  add_reason("missing_other_allele", is.na(other) | !nzchar(other))
  # Each allele predicate is elementwise, so evaluate it once per distinct
  # allele string and index the result back to the rows.
  allele_values <- lapply(list(reference, alternate, effect, other), function(value) {
    if (n >= 64L && is.character(value) && is.null(attributes(value)) &&
        is.loaded("compressor_unique_strings", PACKAGE = "CompreSSoR")) {
      distinct <- .Call("compressor_unique_strings", value, PACKAGE = "CompreSSoR")
      list(value = distinct[[1L]], index = distinct[[2L]])
    } else {
      list(value = value, index = NULL)
    }
  })
  any_allele <- function(predicate) {
    Reduce(`|`, Map(function(allele) {
      present <- !is.na(allele$value)
      hit <- present & predicate(allele$value)
      if (is.null(allele$index)) hit else hit[allele$index]
    }, allele_values), init = rep(FALSE, n))
  }
  add_reason("multiallelic", any_allele(function(value) grepl(",", value, fixed = TRUE)))
  add_reason("symbolic_allele", any_allele(function(value) grepl("^<|^\\*", value)))
  add_reason("indel", any_allele(function(value) nchar(value) != 1L))
  add_reason("invalid_allele", any_allele(function(value) !grepl("^[ACGT]$", value)))
  add_reason("same_alleles", (
    !is.na(reference) & !is.na(alternate) & reference == alternate
  ) | (
    !is.na(effect) & !is.na(other) & effect == other
  ))
  orientation_defined <- !is.na(reference) & !is.na(alternate) &
    !is.na(effect) & !is.na(other)
  add_reason("orientation_mismatch", orientation_defined &
               (effect != alternate | other != reference))

  numeric_fields <- c("beta", "standard_error", "z", "effect_allele_frequency",
                      "p_value", "odds_ratio", "sample_size", "info")
  for (field in numeric_fields) {
    value <- suppressWarnings(as.numeric(values(field, NA_real_)))
    malformed <- row_index_mask(n, parse_failures[[field]] %||% integer())
    add_reason(paste0("malformed_", field), malformed)
    add_reason(paste0("nonfinite_", field), !is.na(value) & !is.finite(value))
  }
  beta <- suppressWarnings(as.numeric(values("beta", NA_real_)))
  se <- suppressWarnings(as.numeric(values("standard_error", NA_real_)))
  z <- suppressWarnings(as.numeric(values("z", NA_real_)))
  eaf <- suppressWarnings(as.numeric(values("effect_allele_frequency", NA_real_)))
  p_value <- suppressWarnings(as.numeric(values("p_value", NA_real_)))
  odds_ratio <- suppressWarnings(as.numeric(values("odds_ratio", NA_real_)))
  sample_size <- suppressWarnings(as.numeric(values("sample_size", NA_real_)))
  info <- suppressWarnings(as.numeric(values("info", NA_real_)))
  add_reason("invalid_standard_error", !is.na(se) & is.finite(se) & se <= 0)
  add_reason("invalid_effect_allele_frequency", !is.na(eaf) & is.finite(eaf) &
               (eaf < 0 | eaf > 1))
  add_reason("invalid_p_value", !is.na(p_value) & is.finite(p_value) &
               (p_value < 0 | p_value > 1))
  add_reason("invalid_odds_ratio", !is.na(odds_ratio) & is.finite(odds_ratio) & odds_ratio <= 0)
  add_reason("invalid_sample_size", !is.na(sample_size) & is.finite(sample_size) & sample_size <= 0)
  add_reason("invalid_info", !is.na(info) & is.finite(info) & (info < 0 | info > 1))
  if ("minus_log10_p" %in% names(data)) {
    lp <- suppressWarnings(as.numeric(data$minus_log10_p))
    add_reason("invalid_minus_log10_p", !is.na(lp) & is.finite(lp) & lp < 0)
    add_reason("nonfinite_minus_log10_p", !is.na(lp) & !is.finite(lp))
  }
  if (isTRUE(require_statistics)) {
    add_reason("missing_statistics", !is.finite(beta) | !is.finite(z) |
                 !is.finite(se) | se <= 0)
  }

  valid_key <- known_chromosome & is.finite(position) & position >= 1 &
    position == floor(position) & position <= chromosome_length &
    !is.na(other) & !is.na(effect) & other %in% c("A", "C", "G", "T") &
    effect %in% c("A", "C", "G", "T") & other != effect
  key <- rep(NA_real_, n)
  if (any(valid_key)) {
    identity <- compressor_encode_variant_identity(
      chromosome[valid_key], position[valid_key], other[valid_key], effect[valid_key],
      build = input_build
    )
    key[valid_key] <- compressor_identity_code(identity$global_position,
                                               identity$substitution)
  }
  # Rows whose key occurs more than once: equal to
  # duplicated(key) | duplicated(key, fromLast = TRUE) with one hashing pass.
  repeated <- duplicated(key)
  duplicate <- !is.na(key) & (repeated | key %in% key[repeated])
  add_reason("duplicate_variant", duplicate)
  report$rejection_counts <- sort(rejection_counts, decreasing = TRUE)
  report$counts <- report$rejection_counts
  report$rejections <- report$rejection_counts[report$rejection_counts > 0L]
  report$examples <- examples
  report$rejected_rows <- as.integer(sum(invalid))
  report$valid <- !any(invalid)
  if (identical(detail, "full")) {
    reason_names <- names(reasons)
    reason_text <- rep("", n)
    for (name in reason_names) {
      rows <- reasons[[name]]
      reason_text[rows] <- ifelse(
        nzchar(reason_text[rows]), paste0(reason_text[rows], ",", name), name
      )
    }
    report$row_status$structurally_valid <- !invalid
    report$row_status$reasons <- reason_text
    report$canonical_key <- key
    report$invalid_rows <- which(invalid)
    report$duplicate_rows <- which(duplicate)
    report$structurally_valid_rows <- which(!invalid)
  } else {
    # Private, short-lived state consumed by apply_structural_qc(). It is
    # removed before the compact report is attached to a store or returned.
    report$internal <- list(
      invalid = invalid,
      duplicate = duplicate,
      canonical_key = key,
      reason_count = reason_count
    )
  }
  report
}


reference_pvalue_resolve <- function(data) {
  n <- nrow(data)
  z <- if ("z" %in% names(data)) suppressWarnings(as.numeric(data$z)) else {
    rep(NA_real_, n)
  }
  derived <- 2 * stats::pnorm(-abs(z))
  source_present <- attr(data, "p_value_source_present", exact = TRUE)
  source_present <- if (is.null(source_present)) "p_value" %in% names(data) else
    isTRUE(source_present)
  supplied <- if (source_present && "p_value" %in% names(data)) {
    suppressWarnings(as.numeric(data$p_value))
  } else {
    rep(NA_real_, n)
  }
  supplied_valid <- source_present & is.finite(supplied) & supplied >= 0 & supplied <= 1
  supplied_missing <- source_present & is.na(supplied)
  supplied_invalid <- source_present & !supplied_missing & !supplied_valid
  derived_valid <- !supplied_valid & is.finite(derived)
  p_value <- derived
  p_value[supplied_valid] <- supplied[supplied_valid]
  fallback <- source_present & !supplied_valid
  unresolved <- !is.finite(p_value)
  source <- if (!source_present || !any(supplied_valid)) {
    "derived_from_z"
  } else if (any(fallback)) {
    "supplied_with_z_fallback"
  } else {
    "supplied"
  }
  list(
    p_value = p_value,
    source = source,
    supplied_rows = as.integer(sum(supplied_valid)),
    derived_rows = as.integer(sum(derived_valid)),
    fallback_rows = as.integer(sum(fallback)),
    missing_rows = as.integer(sum(supplied_missing)),
    invalid_rows = as.integer(sum(supplied_invalid)),
    unresolved_rows = as.integer(sum(unresolved))
  )
}


reference_position_gaps <- function(position, block_rows = PCODEC_NATIVE_BLOCK_ROWS) {
  block_rows <- as.integer(block_rows)
  if (length(block_rows) != 1L || is.na(block_rows) || block_rows < 1L) {
    stop("native Pcodec position block_rows must be positive", call. = FALSE)
  }
  position <- as.numeric(position)
  n <- length(position)
  if (n && any(!is.finite(position) | position < 0 | position > 4294967295 |
               position != floor(position))) {
    stop("native Pcodec positions must be uint32 values", call. = FALSE)
  }
  if (n > 1L && any(diff(position) < 0)) {
    stop("native Pcodec positions must be sorted", call. = FALSE)
  }
  gaps <- numeric(n)
  if (n) {
    blocks <- ceiling(n / block_rows)
    for (block in seq_len(blocks)) {
      start <- (block - 1L) * block_rows + 1L
      stop <- min(n, block * block_rows)
      gaps[start] <- 0
      if (stop > start) gaps[(start + 1L):stop] <- diff(position[start:stop])
    }
  }
  gaps
}

