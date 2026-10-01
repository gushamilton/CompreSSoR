# Threshold/region candidate extraction for native Pcodec stores.
#
# Candidate membership is defined as the full-read convention
# `p_value <= threshold`, where p is the p-value reconstructed from the lossy
# stored Z (the native decoder's erfc table for central bins and erfc of the
# exact float32 Z for exceptions).  Three strategies reach that set at
# different cost; see read_candidates().

candidates_semantic <- function(store) {
  semantic <- store$manifest$semantic_codec %||% list()
  z_range <- as.numeric(unlist(semantic$z_range %||% c(-3.5, 3.5)))
  se_range <- as.numeric(unlist(semantic$se_residual_range %||% c(-1, 1)))
  list(
    z_range = z_range,
    z_count = as.integer(semantic$z_count %||% 510L),
    se_count = as.integer(semantic$se_count %||% PCODEC_NATIVE_SE_COUNT),
    eaf_count = as.integer(semantic$eaf_count %||% 255L),
    z_bits = as.integer(semantic$z_bits %||% 9L),
    se_bits = as.integer(semantic$se_bits %||% PCODEC_NATIVE_SE_BITS),
    eaf_bits = as.integer(semantic$eaf_bits %||% 8L),
    se_range = se_range
  )
}

# p-values via the native decoder's own erfc code path, so that membership
# decisions are bit-identical to a full read.  `codes` are Z codes (central
# bins use the decoder's p table); `exception_z` are exact Z values that go
# through the exception branch.  Returns list(central, exception).
candidates_native_p <- function(semantic, codes = integer(), exception_z = numeric()) {
  n <- length(codes) + length(exception_z)
  z_codes <- c(as.integer(codes), rep.int(semantic$z_count + 1L, length(exception_z)))
  decoded <- .Call(
    "compressor_decode_native", z_codes, integer(n), integer(n),
    semantic$z_range[1], semantic$z_range[2], semantic$z_count,
    semantic$se_count, semantic$eaf_count, semantic$z_bits, semantic$se_bits,
    semantic$eaf_bits, PCODEC_NATIVE_SE_CENTER_ROWS, numeric(),
    seq_along(exception_z) + length(codes) - 1L, as.numeric(exception_z),
    numeric(length(exception_z)), numeric(length(exception_z)),
    rep.int(1L, length(exception_z)), FALSE, TRUE,
    semantic$se_range[1], semantic$se_range[2],
    PACKAGE = "CompreSSoR"
  )
  decoded$p_value
}

# Smallest reconstructed p of any non-exception row.  Strictly greater than
# every Z-exception row's p, because exceptions have |z| outside the central
# range and the outermost central bin centre lies half a step inside it.
candidates_central_edge_p <- function(semantic) {
  step <- (semantic$z_range[2] - semantic$z_range[1]) / semantic$z_count
  z_edge <- min(abs(semantic$z_range)) - step / 2
  list(z_edge = z_edge,
       p = candidates_native_p(semantic, codes = integer(),
                               exception_z = z_edge))
}

candidates_validate_threshold <- function(x) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x < 0 || x > 1) {
    stop("pvalue_threshold must be one finite number between 0 and 1",
         call. = FALSE)
  }
  as.numeric(x)
}

candidates_open <- function(store) {
  store <- if (inherits(store, "compressor_store")) store else
    pcodec_open_store_cached(store)
  if (!identical(store$manifest$backend, "pcodec")) {
    stop("read_candidates requires a native Pcodec store", call. = FALSE)
  }
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  if (!store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS) {
    stop("this CompreSSoR build reads native 0.4 stores only", call. = FALSE)
  }
  store
}

# Returns list(rows (sorted zero-based), strategy).
candidates_select_rows <- function(store, index, threshold, threads) {
  flag <- (store$manifest$domains %||% list())$pvalue_flag
  if (is.list(flag) && isTRUE(as.numeric(flag$threshold) == threshold)) {
    rows <- candidates_read_flag_rows(store, threads)
    # Reconstructed p for flagged rows, taken from the exception sidecar when
    # the row is a Z exception (the usual case); NA otherwise.
    e <- pcodec_native_read_all_exceptions(store, index)
    e <- e[bitwAnd(as.integer(e$flags), 1L) != 0L, , drop = FALSE]
    p <- rep.int(NA_real_, length(rows))
    at <- match(rows, e$row)
    if (any(!is.na(at))) {
      p[!is.na(at)] <- candidates_native_p(candidates_semantic(store),
                                           exception_z = e$z[at[!is.na(at)]])
    }
    return(list(rows = rows, p = p, strategy = "pvalue_flag"))
  }
  semantic <- candidates_semantic(store)
  edge <- candidates_central_edge_p(semantic)
  # Margin guards against ulp-level differences between the analytic bin edge
  # and the decoder table; falling back is always safe.
  exceptions_only <- is.finite(edge$p) && threshold < edge$p * (1 - 1e-9) &&
    identical(semantic$z_range[1], -semantic$z_range[2])
  if (exceptions_only) {
    e <- pcodec_native_read_all_exceptions(store, index)
    zexc <- bitwAnd(as.integer(e$flags), 1L) != 0L
    e <- e[zexc, , drop = FALSE]
    p <- candidates_native_p(semantic, exception_z = e$z)
    hit <- !is.na(p) & p <= threshold
    o <- base::order(e$row[hit])
    return(list(rows = as.integer(e$row[hit][o]), p = p[hit][o],
                strategy = "z_exceptions"))
  }
  codes <- pcodec_native_read_native_codes(store, index, "z", threads = threads)
  p_table <- candidates_native_p(semantic, codes = seq_len(semantic$z_count) - 1L)
  hit_code <- c(!is.na(p_table) & p_table <= threshold, FALSE, FALSE)
  central <- which(hit_code[codes$z + 1L])
  rows <- central - 1L
  p <- p_table[codes$z[central] + 1L]
  ex <- codes$exceptions
  zexc <- bitwAnd(as.integer(ex$flags), 1L) != 0L
  if (any(zexc)) {
    pe <- candidates_native_p(semantic, exception_z = ex$z[zexc])
    hit <- !is.na(pe) & pe <= threshold
    rows <- c(rows, ex$row[zexc][hit])
    p <- c(p, pe[hit])
  }
  o <- base::order(rows)
  list(rows = as.integer(rows[o]), p = p[o], strategy = "z_stream")
}

# Flag rows with one open connection per worker (same validation as
# read_pvalue_flag(): binary, row-aligned).
candidates_read_flag_rows <- function(store, threads) {
  domain <- pcodec_native_pvalue_flag_domain(store)
  blocks <- domain$blocks
  path <- file.path(store$path, domain$file)
  starts <- vapply(blocks, function(b) as.numeric(b$row_start), numeric(1))
  parts <- pcodec_parallel_lapply(
    candidates_split(seq_along(blocks), threads), function(chunk) {
      con <- file(path, open = "rb")
      on.exit(close(con), add = TRUE)
      lapply(chunk, function(b) {
        loc <- blocks[[b]]
        seek(con, where = as.numeric(loc$offset), origin = "start")
        blob <- readBin(con, raw(), n = as.integer(loc$length), endian = "little")
        if (length(blob) != as.integer(loc$length)) {
          stop("native Pcodec p-value flag payload is truncated", call. = FALSE)
        }
        v <- pcodec_native_decompress(blob, as.integer(loc$values), "u8")
        if (any(v > 1L)) {
          stop("native Pcodec p-value flag payload is not a binary row-aligned stream",
               call. = FALSE)
        }
        as.integer(starts[b] + which(v != 0L) - 1)
      })
    }, threads = threads)
  as.integer(unlist(parts, use.names = FALSE))
}

candidates_exact_ranks <- function(store, rows) {
  domain <- pcodec_native_pvalue_order_domain(store)
  blocks <- domain$blocks
  stops <- vapply(blocks, function(b) as.numeric(b$row_stop), numeric(1))
  starts <- vapply(blocks, function(b) as.numeric(b$row_start), numeric(1))
  ids <- findInterval(rows, stops) + 1L
  ranks <- integer(length(rows))
  for (b in unique(ids)) {
    loc <- blocks[[b]]
    blob <- pcodec_native_read_blob(file.path(store$path, domain$file),
                                    loc$offset, loc$length)
    values <- pcodec_native_decompress(blob, as.integer(loc$values), "u32")
    sel <- ids == b
    ranks[sel] <- as.integer(values[rows[sel] - starts[b] + 1L])
  }
  ranks
}

#' Read threshold candidates from a native Pcodec store
#'
#' Returns the rows whose reconstructed p-value (`2 * pnorm(-abs(Z))` from
#' the stored lossy Z, evaluated exactly as `read_sumstats()` does for a full
#' read) satisfies `p <= pvalue_threshold`, optionally restricted to a region,
#' without decoding the whole store when that is avoidable.
#'
#' The strategy is chosen automatically and recorded in
#' `attr(x, "candidate_strategy")`:
#'
#' 1. `"pvalue_flag"` when the store has a flag domain whose threshold equals
#'    `pvalue_threshold` exactly. Membership is then the writer-time flag
#'    (which follows a supplied p-value when one was available).
#' 2. `"z_exceptions"` when `pvalue_threshold` is below the p-value of the
#'    outermost central Z bin centre (about 4.8e-4 for the standard profile).
#'    Every row beyond the central Z range is stored as an exact float32 Z
#'    exception, so candidates are precisely the Z-exception rows meeting the
#'    threshold and only the small exception sidecar is decoded.
#' 3. `"z_stream"` otherwise: only the Z stream is decoded and masked through
#'    a per-code p lookup.
#'
#' @param store A native Pcodec store object or path.
#' @param pvalue_threshold Inclusive threshold in `[0, 1]`.
#' @param region Optional `"chr:start-end"` string or `c(chr, start, end)`
#'   (inclusive base-pair range).
#' @param columns Output columns (as for [read_sumstats()]). The default is
#'   chromosome, position, alleles and `p_value`. The zero-based native `row`
#'   id is always returned first, and `exact_rank` is appended when
#'   `order = "exact"`.
#' @param order `"none"` (native row order), `"reconstructed"` (ascending
#'   reconstructed p, ties by row; approximate, see
#'   `attr(x, "candidate_order")`) or `"exact"` (ascending exact rank from the
#'   `pvalue_order` domain; an error if the domain is absent or covers a
#'   smaller threshold). Candidates without an exact rank (possible only
#'   because reconstructed and exact p differ) get `NA` and sort last.
#' @param threads Decoder threads.
#' @return A data frame with attributes `candidate_strategy`,
#'   `candidate_threshold` and `candidate_order`.
#' @export
read_candidates <- function(store, pvalue_threshold, region = NULL,
                            columns = NULL,
                            order = c("none", "reconstructed", "exact"),
                            threads = 1L) {
  order <- match.arg(order)
  threshold <- candidates_validate_threshold(pvalue_threshold)
  threads <- pcodec_validate_threads(threads)
  store <- candidates_open(store)
  if (!is.null(columns) && (!length(columns) || anyNA(columns))) {
    stop("columns must contain at least one column name", call. = FALSE)
  }
  out_columns <- if (is.null(columns)) {
    c("chromosome", "base_pair_location", "effect_allele", "other_allele",
      "p_value")
  } else unique(as.character(columns))
  if (identical(order, "exact")) {
    domain <- pcodec_native_pvalue_order_domain(store)
    if (threshold > as.numeric(domain$threshold)) {
      stop("the exact p-value ordering domain covers only p <= ",
           format(as.numeric(domain$threshold), scientific = TRUE),
           "; requested threshold ", format(threshold, scientific = TRUE),
           call. = FALSE)
    }
  }
  index <- pcodec_native_read_index(store)
  picked <- candidates_select_rows(store, index, threshold, threads)
  rows <- picked$rows
  allowed <- c("global_position", "substitution", "chromosome", "base_pair_location",
               "reference_allele", "alternate_allele", "effect_allele", "other_allele",
               "z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  unknown <- setdiff(out_columns, allowed)
  if (length(unknown)) {
    stop("unknown output column(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  fetch <- unique(c(out_columns, if (identical(order, "reconstructed")) "p_value"))
  build <- compressor_normalize_build(store$manifest$genome_build %||% "GRCh38")
  range <- candidates_region_range(region, build)
  # p_value is attached from the selection step; Z is needed only when the
  # selection (the flag) did not already provide the reconstructed p.
  wanted <- setdiff(fetch, "p_value")
  if ("p_value" %in% fetch && anyNA(picked$p)) wanted <- unique(c(wanted, "z"))
  data <- if (length(rows)) {
    candidates_fetch(store, index, rows, range, build, wanted, threads)
  } else NULL
  if (is.null(data)) {
    data <- pcodec_native_empty_result(unique(c(wanted, "z")))
    data$p_value <- numeric()
  }
  if (nrow(data) && "p_value" %in% fetch) {
    # Use the full-read decoder's erfc p (the sparse reader's 2 * pnorm can
    # differ in the last ulp). Flag-selected central-bin rows (rare: only for
    # flag thresholds above ~4.8e-4) fall back to erfc of the decoded Z,
    # which can differ from the decoder table by an ulp.
    p <- picked$p[match(data$row, rows)]
    if (anyNA(p)) {
      miss <- is.na(p)
      p[miss] <- candidates_native_p(candidates_semantic(store),
                                     exception_z = data$z[miss])
    }
    data$p_value <- p
  }
  if (identical(order, "exact")) {
    data$exact_rank <- if (nrow(data)) {
      r <- candidates_exact_ranks(store, data$row)
      ifelse(r > 0L, r, NA_integer_)
    } else integer()
    o <- base::order(is.na(data$exact_rank), data$exact_rank, data$row,
                     method = "radix")
    data <- data[o, , drop = FALSE]
  } else if (identical(order, "reconstructed")) {
    data <- data[base::order(data$p_value, data$row, method = "radix"), ,
                 drop = FALSE]
  }
  keep <- c("row", out_columns, if (identical(order, "exact")) "exact_rank")
  data <- data[keep]
  row.names(data) <- NULL
  attr(data, "candidate_strategy") <- picked$strategy
  attr(data, "candidate_threshold") <- threshold
  attr(data, "candidate_order") <- switch(order,
    none = "native_row", reconstructed = "reconstructed_p_approximate",
    exact = "exact_rank")
  data
}

#' Read threshold candidates from several native Pcodec stores
#'
#' Applies [read_candidates()] to each store.
#'
#' @inheritParams read_candidates
#' @param stores A non-empty list or character vector of native Pcodec stores.
#' @param pvalue_threshold One threshold, or one per store.
#' @param bind If `FALSE` (default, as for [read_sumstats_batch()]) return a
#'   list of data frames in the order of `stores` (names preserved). If `TRUE`
#'   return one data frame with a leading `store` column holding the store
#'   names (or positions when unnamed).
#' @param threads Number of stores decoded in parallel on Unix-like systems;
#'   each store is read with one decoder thread when several stores are given.
#' @return A list or data frame; `attr(x, "candidate_strategy")` is the
#'   per-store strategy vector.
#' @export
read_candidates_batch <- function(stores, pvalue_threshold, region = NULL,
                                  columns = NULL,
                                  order = c("none", "reconstructed", "exact"),
                                  threads = 1L, bind = FALSE) {
  order <- match.arg(order)
  if (is.character(stores)) stores <- as.list(stores)
  if (!is.list(stores) || !length(stores)) {
    stop("stores must be a non-empty list or character vector", call. = FALSE)
  }
  if (!length(pvalue_threshold) %in% c(1L, length(stores))) {
    stop("pvalue_threshold must have length 1 or one value per store",
         call. = FALSE)
  }
  thresholds <- rep_len(pvalue_threshold, length(stores))
  threads <- pcodec_validate_threads(threads)
  inner <- if (length(stores) > 1L) 1L else threads
  pieces <- pcodec_parallel_lapply(seq_along(stores), function(i) {
    tryCatch(
      read_candidates(stores[[i]], thresholds[[i]], region = region,
                      columns = columns, order = order, threads = inner),
      error = function(e) e
    )
  }, threads = threads)
  failed <- vapply(pieces, function(x) !is.data.frame(x), logical(1))
  if (any(failed)) {
    first <- which(failed)[1L]
    msg <- if (inherits(pieces[[first]], "error")) conditionMessage(pieces[[first]])
           else "worker failed"
    stop("failed to read candidates from store ", first, ": ", msg, call. = FALSE)
  }
  labels <- names(stores)
  if (is.null(labels)) labels <- as.character(seq_along(stores))
  strategy <- stats::setNames(
    vapply(pieces, attr, character(1), "candidate_strategy"), labels)
  if (!isTRUE(bind)) {
    names(pieces) <- names(stores)
    attr(pieces, "candidate_strategy") <- strategy
    return(pieces)
  }
  combined <- do.call(rbind, lapply(seq_along(pieces), function(i) {
    x <- pieces[[i]]
    cbind(data.frame(store = rep.int(labels[i], nrow(x)),
                     stringsAsFactors = FALSE), x)
  }))
  row.names(combined) <- NULL
  attr(combined, "candidate_strategy") <- strategy
  attr(combined, "candidate_threshold") <- thresholds
  combined
}

# ---------------------------------------------------------------------------
# Sparse row fetch.  Equivalent to pcodec_native_read_store(variants = rows)
# but each worker opens every payload file once and only touches the key and
# value blocks that contain candidate rows (and overlap the region).

candidates_split <- function(x, parts) {
  if (!length(x)) return(list())
  parts <- max(1L, min(as.integer(parts), length(x)))
  unname(split(x, ceiling(seq_along(x) * parts / length(x))))
}

candidates_region_range <- function(region, build) {
  if (is.null(region)) return(NULL)
  bounds <- read_region_bounds(region)
  chromosome <- toupper(sub("^CHR", "", as.character(bounds$chromosome),
                            ignore.case = TRUE))
  lengths <- compressor_chromosome_lengths(build)
  if (!chromosome %in% names(lengths)) {
    stop("unsupported region chromosome", call. = FALSE)
  }
  offset <- pcodec_native_offsets(build)[match(chromosome, names(lengths))]
  c(offset + bounds$start - 1, offset + bounds$end - 1)
}

candidates_fetch <- function(store, index, rows, range, build, wanted,
                             threads) {
  manifest <- store$manifest
  semantic <- manifest$semantic_codec %||% list()
  key_blocks <- pcodec_native_index_blocks(index, "key")
  value_blocks <- pcodec_native_index_blocks(index, "value")
  num <- function(blocks, field) vapply(blocks, function(b) as.numeric(b[[field]]),
                                        numeric(1))
  key_stops <- num(key_blocks, "row_stop")
  key_starts <- num(key_blocks, "row_start")
  key_first <- num(key_blocks, "first_position")
  key_last <- num(key_blocks, "last_position")
  delta <- identical(index$position_encoding, "delta_u32_within_block")
  identity <- c("global_position", "substitution", "chromosome",
                "base_pair_location", "reference_allele", "alternate_allele",
                "effect_allele", "other_allele")
  need_identity <- !is.null(range) || any(identity %in% wanted)
  key_id <- findInterval(rows, key_stops) + 1L
  todo <- if (need_identity) unique(key_id) else integer()
  if (!is.null(range)) todo <- todo[key_last[todo] >= range[1] & key_first[todo] <= range[2]]
  open_stream <- function(spec) file(file.path(store$path, spec$file), open = "rb")
  read_at <- function(con, location) {
    seek(con, where = as.numeric(location$offset), origin = "start")
    blob <- readBin(con, raw(), n = as.integer(location$length), endian = "little")
    if (length(blob) != as.integer(location$length)) {
      stop("native Pcodec stream is truncated", call. = FALSE)
    }
    blob
  }
  keys <- pcodec_parallel_lapply(candidates_split(todo, threads), function(chunk) {
    pos_con <- open_stream(index$streams$position)
    sub_con <- open_stream(index$streams$substitution)
    on.exit({ close(pos_con); close(sub_con) }, add = TRUE)
    lapply(chunk, function(b) {
      mine <- rows[key_id == b]
      local <- mine - key_starts[b] + 1
      pos_loc <- index$streams$position$blocks[[b]]
      position <- pcodec_native_decompress(read_at(pos_con, pos_loc),
                                           as.integer(pos_loc$values), "u32")
      position <- if (delta) cumsum(position) + key_first[b] else as.numeric(position)
      position <- position[local]
      keep <- if (is.null(range)) rep.int(TRUE, length(mine)) else
        position >= range[1] & position <= range[2]
      if (!any(keep)) return(NULL)
      sub_loc <- index$streams$substitution$blocks[[b]]
      substitution <- pcodec_native_decompress(read_at(sub_con, sub_loc),
                                               as.integer(sub_loc$values), "u8")
      list(row = mine[keep], position = position[keep],
           substitution = as.integer(substitution[local[keep]]))
    })
  }, threads = threads)
  if (need_identity) {
    keys <- unlist(keys, recursive = FALSE, use.names = FALSE)
    keys <- keys[!vapply(keys, is.null, logical(1))]
    if (!length(keys)) return(NULL)
    sel_row <- unlist(lapply(keys, `[[`, "row"), use.names = FALSE)
    sel_position <- unlist(lapply(keys, `[[`, "position"), use.names = FALSE)
    sel_substitution <- unlist(lapply(keys, `[[`, "substitution"), use.names = FALSE)
    o <- order(sel_row)
    sel_row <- as.integer(sel_row[o])
    sel_position <- sel_position[o]
    sel_substitution <- sel_substitution[o]
  } else {
    sel_row <- as.integer(rows)
  }

  out <- list(row = sel_row)
  if (any(identity %in% wanted)) {
    cols <- pcodec_native_key_columns(sel_position, sel_substitution, build = build)
    cols$global_position <- sel_position
    cols$substitution <- sel_substitution
    for (nm in intersect(identity, wanted)) out[[nm]] <- cols[[nm]]
  }
  need_z <- any(c("z", "beta") %in% wanted)
  need_se <- any(c("standard_error", "beta") %in% wanted)
  need_eaf <- "effect_allele_frequency" %in% wanted || need_se
  needed <- c(if (need_z) "z", if (need_eaf) "eaf", if (need_se) "se")
  if (length(needed)) {
    value_stops <- num(value_blocks, "row_stop")
    value_starts <- num(value_blocks, "row_start")
    value_id <- findInterval(sel_row, value_stops) + 1L
    centre_rows <- as.integer(semantic$se_center_block_rows %||%
                                PCODEC_NATIVE_SE_CENTER_ROWS)
    centres <- as.numeric(unlist(semantic$block_centers_log2_residual))
    codec <- index$exceptions$codec %||% "raw"
    parts <- pcodec_parallel_lapply(
      candidates_split(unique(value_id), threads), function(chunk) {
        cons <- lapply(stats::setNames(needed, needed), function(s) open_stream(index$streams[[s]]))
        exc_con <- open_stream(index$exceptions)
        on.exit({ lapply(cons, close); close(exc_con) }, add = TRUE)
        lapply(chunk, function(b) {
          pick <- value_id == b
          local <- sel_row[pick] - value_starts[b] + 1
          n_block <- value_stops[b] - value_starts[b]
          codes <- lapply(stats::setNames(needed, needed), function(s) {
            loc <- index$streams[[s]]$blocks[[b]]
            pcodec_native_decompress(read_at(cons[[s]], loc), as.integer(loc$values),
                                     if (identical(s, "z")) "u16" else "u8")
          })
          exc <- index$exceptions$blocks[[b]]
          exceptions <- if (as.integer(exc$count)) {
            blob <- read_at(exc_con, exc)
            if (identical(codec, "zstd")) {
              blob <- pcodec_native_zstd_decompress(
                blob, as.integer(exc$raw_length %||% (as.integer(exc$count) * 17L)))
            } else if (!identical(codec, "raw")) {
              stop("unsupported native exception codec: ", codec, call. = FALSE)
            }
            pcodec_native_read_exception_bytes(blob, as.integer(exc$count))
          } else {
            data.frame(row = integer(), z = numeric(), log2se = numeric(),
                       eaf = numeric(), flags = integer())
          }
          decoded <- pcodec_native_decode_values(
            codes, exceptions, floor(value_starts[b] / centre_rows) + 1L,
            centres, value_starts[b], n_block, needed, semantic)
          lapply(decoded, function(v) v[local])
        })
      }, threads = threads)
    parts <- unlist(parts, recursive = FALSE, use.names = FALSE)
    pull <- function(field) unlist(lapply(parts, `[[`, field), use.names = FALSE)
    if (need_z) out$z <- pull("z")
    if (need_se) out$standard_error <- pull("se")
    if (need_eaf) out$effect_allele_frequency <- pull("eaf")
    if ("beta" %in% wanted) out$beta <- out$z * out$standard_error
  }
  out <- as.data.frame(out[c("row", intersect(c(identity, "z", "beta",
    "standard_error", "effect_allele_frequency"), names(out)))],
    stringsAsFactors = FALSE)
  out
}
