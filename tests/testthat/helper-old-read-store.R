# Reference copy of pcodec_native_read_store before the vectorised key matching (perf T1).
pcodec_native_read_store_reference <- function(store, region = NULL, variants = NULL,
                                      columns = NULL, threads = NULL) {
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  threads <- pcodec_native_default_threads(region = region, variants = variants,
                                            threads = threads)
  manifest <- store$manifest
  build <- compressor_normalize_build(manifest$genome_build %||% "GRCh38")
  index <- pcodec_native_read_index(store)
  n <- as.integer(manifest$n_rows %||% manifest$rows)
  if (!is.null(columns) && !length(columns)) stop("columns must contain at least one column name", call. = FALSE)
  requested <- if (is.null(columns)) {
    c("chromosome", "base_pair_location", "effect_allele", "other_allele",
      "z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  } else unique(as.character(columns))
  allowed <- c("global_position", "substitution", "chromosome", "base_pair_location",
               "reference_allele", "alternate_allele", "effect_allele", "other_allele",
               "z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  unknown <- setdiff(requested, allowed)
  if (length(unknown)) stop("unknown output column(s): ", paste(unknown, collapse = ", "), call. = FALSE)

  row_targets <- NULL
  key_targets <- NULL
  if (!is.null(variants)) {
    if (is.character(variants)) {
      key_targets <- pcodec_native_target_keys(unique(trimws(variants)), build = build)
    } else {
      row_targets <- unique(as.integer(variants))
      if (anyNA(row_targets) || any(row_targets < 0L | row_targets >= n)) {
        stop("variants must be valid zero-based row IDs", call. = FALSE)
      }
    }
  }
  identity_needed <- is.null(columns) || any(c("global_position", "substitution",
                                                "chromosome", "base_pair_location",
                                                "reference_allele", "alternate_allele",
                                                "effect_allele", "other_allele") %in% requested) ||
    !is.null(region) || !is.null(key_targets)
  need_z <- is.null(columns) || any(c("z", "beta", "p_value") %in% requested)
  need_se <- is.null(columns) || any(c("standard_error", "beta") %in% requested)
  need_eaf <- is.null(columns) || "effect_allele_frequency" %in% requested || need_se
  needed <- c(if (need_z) "z", if (need_eaf) "eaf", if (need_se) "se")

  lower <- upper <- NULL
  if (!is.null(region)) {
    bounds <- read_region_bounds(region)
    chromosome <- toupper(sub("^CHR", "", as.character(bounds$chromosome), ignore.case = TRUE))
    lengths <- compressor_chromosome_lengths(build)
    if (!chromosome %in% names(lengths)) stop("unsupported region chromosome", call. = FALSE)
    offsets <- pcodec_native_offsets(build)
    lower <- offsets[match(chromosome, names(lengths))] + bounds$start - 1
    upper <- offsets[match(chromosome, names(lengths))] + bounds$end - 1
  }
  if (is.null(region) && is.null(variants)) {
    output <- pcodec_native_full_read(store, index, requested, identity_needed,
                                      need_z, need_se, need_eaf, threads = threads,
                                      build = build)
    attr(output, "source_bytes_read") <- NA_real_
    return(output)
  }

  source_bytes <- 0
  value_blocks <- pcodec_native_index_blocks(index, "value")
  selected_rows <- selected_position <- selected_substitution <- NULL
  if (identity_needed) {
    key_blocks <- pcodec_native_index_blocks(index, "key")
    key_candidates <- seq_along(key_blocks)
    if (!is.null(lower)) {
      key_candidates <- key_candidates[vapply(key_blocks, function(block) {
        as.numeric(block$last_position) >= lower &&
          as.numeric(block$first_position) <= upper
      }, logical(1))]
    }
    if (!is.null(key_targets)) {
      target_positions <- as.numeric(key_targets$position)
      key_candidates <- key_candidates[vapply(key_blocks, function(block) {
        any(target_positions >= as.numeric(block$first_position) &
              target_positions <= as.numeric(block$last_position))
      }, logical(1))]
    }
    if (!is.null(row_targets)) {
      key_candidates <- intersect(key_candidates,
                                  pcodec_native_block_ids_for_rows(index, row_targets, "key"))
    }
    read_key_candidate <- function(key_block) {
      meta <- key_blocks[[key_block]]
      rows <- as.integer(meta$row_start):(as.integer(meta$row_stop) - 1L)
      position <- pcodec_native_read_stream_block(store, index, "position", key_block)
      substitution <- pcodec_native_read_stream_block(store, index, "substitution", key_block)
      bytes <-
        as.numeric(index$streams$position$blocks[[key_block]]$length) +
        as.numeric(index$streams$substitution$blocks[[key_block]]$length)
      keep <- rep(TRUE, length(rows))
      if (!is.null(lower)) keep <- keep & position >= lower & position <= upper
      if (!is.null(row_targets)) keep <- keep & rows %in% row_targets
      if (!is.null(key_targets)) {
        keep <- keep & paste(position, substitution, sep = ":") %in%
          paste(key_targets$position, key_targets$substitution, sep = ":")
      }
      part <- if (any(keep)) {
        data.frame(
          row = rows[keep], position = position[keep],
          substitution = substitution[keep], stringsAsFactors = FALSE)
      } else NULL
      list(part = part, source_bytes = bytes)
    }
    key_results <- pcodec_parallel_lapply(key_candidates, read_key_candidate,
                                          threads = threads)
    source_bytes <- source_bytes + sum(vapply(key_results,
                                              function(result) result$source_bytes,
                                              numeric(1)))
    key_parts <- lapply(key_results, `[[`, "part")
    key_parts <- key_parts[!vapply(key_parts, is.null, logical(1))]
    if (!length(key_parts)) return(pcodec_native_empty_result(columns))
    selected <- do.call(rbind, key_parts)
    selected <- selected[order(selected$row), , drop = FALSE]
    selected_rows <- as.integer(selected$row)
    selected_position <- as.numeric(selected$position)
    selected_substitution <- as.integer(selected$substitution)
    candidate_blocks <- pcodec_native_block_ids_for_rows(index, selected_rows, "value")
  } else {
    candidate_blocks <- seq_along(value_blocks)
    if (!is.null(row_targets)) {
      candidate_blocks <- pcodec_native_block_ids_for_rows(index, row_targets, "value")
    }
  }
  if (!length(candidate_blocks)) return(pcodec_native_empty_result(columns))

  decode_value_block <- function(block) {
    meta <- value_blocks[[block]]
    row_start <- as.integer(meta$row_start)
    row_stop <- as.integer(meta$row_stop)
    rows <- row_start:(row_stop - 1L)
    block_source_bytes <- 0
    if (identity_needed) {
      keep <- rows %in% selected_rows
    } else {
      keep <- if (is.null(row_targets)) rep(TRUE, length(rows)) else rows %in% row_targets
    }
    if (!any(keep)) return(list(part = NULL, source_bytes = block_source_bytes))
    value_codes <- list()
    for (stream in intersect(c("z", "eaf", "se"), needed)) {
      value_codes[[stream]] <- pcodec_native_read_stream_block(store, index, stream, block)
      block_source_bytes <- block_source_bytes +
        as.numeric(index$streams[[stream]]$blocks[[block]]$length)
    }
    exceptions <- pcodec_native_read_exception_block(store, index, block)
    if (nrow(exceptions)) {
      block_source_bytes <- block_source_bytes +
        as.numeric(index$exceptions$blocks[[block]]$length)
    }
    centre_id <- floor(row_start / as.integer(
      manifest$semantic_codec$se_center_block_rows %||% PCODEC_NATIVE_SE_CENTER_ROWS
    )) + 1L
    decoded <- pcodec_native_decode_values(
      value_codes, exceptions, centre_id,
      as.numeric(unlist(manifest$semantic_codec$block_centers_log2_residual)),
      row_start, length(rows), needed, manifest$semantic_codec
    )
    decoded <- lapply(decoded, function(value) value[keep])
    part <- data.frame(row = rows[keep], stringsAsFactors = FALSE)
    if (identity_needed) {
      selected_index <- match(rows[keep], selected_rows)
      identity_part <- pcodec_native_key_columns(
        selected_position[selected_index], selected_substitution[selected_index],
        build = build)
      part <- cbind(part, as.data.frame(identity_part, stringsAsFactors = FALSE))
    }
    if ("z" %in% names(decoded)) part$z <- decoded$z
    if ("se" %in% names(decoded)) part$standard_error <- decoded$se
    if ("eaf" %in% names(decoded)) part$effect_allele_frequency <- decoded$eaf
    if ("beta" %in% requested) part$beta <- part$z * part$standard_error
    if ("p_value" %in% requested) part$p_value <- 2 * stats::pnorm(-abs(part$z))
    list(part = part, source_bytes = block_source_bytes)
  }
  block_results <- pcodec_parallel_lapply(candidate_blocks, decode_value_block,
                                          threads = threads)
  source_bytes <- source_bytes + sum(vapply(block_results,
                                            function(result) result$source_bytes,
                                            numeric(1)))
  results <- lapply(block_results, `[[`, "part")
  results <- results[!vapply(results, is.null, logical(1))]
  if (!length(results)) return(pcodec_native_empty_result(columns))
  output <- do.call(rbind, results)
  output <- output[order(output$row), , drop = FALSE]
  row.names(output) <- NULL
  if (is.null(columns)) {
    output <- output[c("row", setdiff(c("chromosome", "base_pair_location",
      "reference_allele", "alternate_allele",
      "effect_allele", "other_allele", "z", "beta", "standard_error",
      "effect_allele_frequency", "p_value"), ""))]
  } else {
    missing <- setdiff(requested, names(output))
    if (length(missing)) stop("requested columns are not present: ", paste(missing, collapse = ", "), call. = FALSE)
    output <- output[c("row", requested)]
  }
  attr(output, "source_bytes_read") <- source_bytes
  output
}
environment(pcodec_native_read_store_reference) <- asNamespace("CompreSSoR")
