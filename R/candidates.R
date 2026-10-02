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
  pcodec_native_p_from_codes(semantic, z_codes,
                             seq_along(exception_z) + length(codes) - 1L,
                             as.numeric(exception_z))
}

# THE p-value reconstruction used by every reader (full read, selective read
# and read_candidates): the native decoder's erfc table for central Z codes and
# erfc of the exact float32 Z for Z-exception rows.  `exc_row0` are zero-based
# positions into `z_codes` whose Z is overridden by `exc_z`.
pcodec_native_p_from_codes <- function(semantic, z_codes, exc_row0, exc_z) {
  n <- length(z_codes)
  k <- length(exc_z)
  decoded <- .Call(
    "compressor_decode_native", as.integer(z_codes), integer(n), integer(n),
    semantic$z_range[1], semantic$z_range[2], semantic$z_count,
    semantic$se_count, semantic$eaf_count, semantic$z_bits, semantic$se_bits,
    semantic$eaf_bits, PCODEC_NATIVE_SE_CENTER_ROWS, numeric(),
    as.integer(exc_row0), as.numeric(exc_z),
    numeric(k), numeric(k), rep.int(1L, k), FALSE, TRUE,
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

# Reconstructed p for the given rows, taken from the exception sidecar when
# the row is a Z exception (the usual case); NA otherwise.
candidates_exception_p <- function(store, index, rows) {
  e <- pcodec_native_read_all_exceptions(store, index)
  e <- e[bitwAnd(as.integer(e$flags), 1L) != 0L, , drop = FALSE]
  p <- rep.int(NA_real_, length(rows))
  at <- match(rows, e$row)
  if (any(!is.na(at))) {
    p[!is.na(at)] <- candidates_native_p(candidates_semantic(store),
                                         exception_z = e$z[at[!is.na(at)]])
  }
  p
}

# Returns list(rows (sorted zero-based), strategy).
candidates_select_rows <- function(store, index, threshold, threads,
                                   strategy = "auto") {
  flag <- (store$manifest$domains %||% list())$pvalue_flag
  if (identical(strategy, "exact_order")) {
    domain <- pcodec_native_pvalue_order_domain(store)
    if (!isTRUE(as.numeric(domain$threshold) == threshold)) {
      stop("strategy = \"exact_order\" requires pvalue_threshold equal to the ",
           "pvalue_order domain threshold (",
           format(as.numeric(domain$threshold), scientific = TRUE),
           "); ranks encode order and membership at that threshold only",
           call. = FALSE)
    }
    ranks <- pcodec_native_read_pvalue_order(store, threads = threads)
    rows <- as.integer(which(ranks > 0L) - 1L)
    return(list(rows = rows, p = candidates_exception_p(store, index, rows),
                strategy = "exact_order"))
  }
  if (identical(strategy, "pvalue_flag")) {
    if (!is.list(flag) || !isTRUE(as.numeric(flag$threshold) == threshold)) {
      stop("strategy = \"pvalue_flag\" requires a pvalue_flag domain whose ",
           "threshold equals pvalue_threshold", call. = FALSE)
    }
    rows <- candidates_read_flag_rows(store, threads)
    return(list(rows = rows, p = candidates_exception_p(store, index, rows),
                strategy = "pvalue_flag"))
  }
  semantic <- candidates_semantic(store)
  edge <- candidates_central_edge_p(semantic)
  # Margin guards against ulp-level differences between the analytic bin edge
  # and the decoder table; falling back is always safe.
  exceptions_only <- is.finite(edge$p) && threshold < edge$p * (1 - 1e-9) &&
    identical(semantic$z_range[1], -semantic$z_range[2])
  if (identical(strategy, "exceptions") && !exceptions_only) {
    stop("strategy = \"exceptions\" requires pvalue_threshold below the ",
         "outermost central Z bin p-value (", format(edge$p, digits = 3), ")",
         call. = FALSE)
  }
  if (identical(strategy, "exceptions") ||
      (identical(strategy, "auto") && exceptions_only)) {
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
#' 1. `"z_exceptions"` when `pvalue_threshold` is below the p-value of the
#'    outermost central Z bin centre (about 4.8e-4 for the standard profile).
#'    Every row beyond the central Z range is stored as an exact float32 Z
#'    exception, so candidates are precisely the Z-exception rows meeting the
#'    threshold and only the small exception sidecar is decoded.
#' 2. `"z_stream"` otherwise: only the Z stream is decoded and masked through
#'    a per-code p lookup.
#'
#' Membership is always the reconstructed-p definition unless
#' `strategy = "pvalue_flag"` or `strategy = "exact_order"` is requested
#' explicitly. `"exact_order"` defines membership by the full-precision
#' source (supplied) p-value recorded at write time in the `pvalue_order`
#' domain (every row with source p <= the domain threshold, which
#' `pvalue_threshold` must equal), and so can differ from the reconstructed-p
#' set by a few rows near the threshold (about 0.5% of rows at p <= 0.01 in
#' benchmarks). The returned `p_value` column is still the reconstructed p;
#' `candidate_strategy` is `"exact_order"` and `exact_rank` gives the source
#' order.
#'
#' @param store A native Pcodec store object or path.
#' @param pvalue_threshold Inclusive threshold in `[0, 1]`.
#' @param region Optional `"chr:start-end"` string or `c(chr, start, end)`
#'   (inclusive base-pair range).
#' @param columns Output columns: any [read_sumstats()] column (identity,
#'   `z`, `beta`, `standard_error`, `effect_allele_frequency`, `p_value`) plus
#'   `key`, the canonical variant key
#'   (`compressor_variant_key(chromosome, position, other_allele,
#'   effect_allele)`). The default is chromosome, position, alleles and
#'   `p_value`. The zero-based native `row` id is always returned first, and
#'   `exact_rank` is appended when `order = "exact"`. Only the candidate rows'
#'   key and value blocks are decoded, in one pass; values and `p_value` are
#'   bit-identical to `read_sumstats(store, variants = rows)` and to the
#'   full read.
#' @param order `"none"` (native row order), `"reconstructed"` (ascending
#'   reconstructed p, ties by row; approximate, see
#'   `attr(x, "candidate_order")`) or `"exact"` (ascending exact rank from the
#'   `pvalue_order` domain; an error if the domain is absent or covers a
#'   smaller threshold). The `pvalue_order` blocks are row-aligned pcodec
#'   frames, so ranks are read by decoding only the blocks that contain
#'   candidate rows (not the whole domain). Candidates without an exact rank (possible only
#'   because reconstructed and exact p differ) get `NA` and sort last.
#' @param threads Decoder threads.
#' @param strategy `"auto"` (default; exact, never uses the flag),
#'   `"exceptions"` (exceptions-only; error unless the threshold is below the
#'   exceptions edge), `"z_stream"`, or `"pvalue_flag"` (explicit opt-in: the
#'   writer-time flag membership, which follows a supplied p-value when one was
#'   available and so can differ from the reconstructed-p definition; requires
#'   a flag domain whose threshold equals `pvalue_threshold`), or
#'   `"exact_order"` (explicit opt-in: source-p membership from the
#'   `pvalue_order` domain, whose threshold must equal `pvalue_threshold`).
#' @return A data frame with attributes `candidate_strategy`,
#'   `candidate_threshold` and `candidate_order`.
#' @export
read_candidates <- function(store, pvalue_threshold, region = NULL,
                            columns = NULL,
                            order = c("none", "reconstructed", "exact"),
                            threads = 1L,
                            strategy = c("auto", "exceptions", "z_stream",
                                         "pvalue_flag", "exact_order")) {
  order <- match.arg(order)
  strategy <- match.arg(strategy)
  threads <- pcodec_validate_threads(threads)
  ctx <- candidates_prepare(store, pvalue_threshold, region, columns, order,
                            threads, strategy)
  candidates_finish(ctx, threads)
}

# Stage 1: validate, open and select candidate rows (no value/key decode).
candidates_prepare <- function(store, pvalue_threshold, region, columns, order,
                               threads, strategy) {
  threshold <- candidates_validate_threshold(pvalue_threshold)
  store <- candidates_open(store)
  if (!is.null(columns) && (!length(columns) || anyNA(columns))) {
    stop("columns must contain at least one column name", call. = FALSE)
  }
  out_columns <- if (is.null(columns)) {
    c("chromosome", "base_pair_location", "effect_allele", "other_allele",
      "p_value")
  } else unique(as.character(columns))
  allowed <- c("global_position", "substitution", "chromosome", "base_pair_location",
               "reference_allele", "alternate_allele", "effect_allele", "other_allele",
               "z", "beta", "standard_error", "effect_allele_frequency", "p_value",
               "key")
  unknown <- setdiff(out_columns, allowed)
  if (length(unknown)) {
    stop("unknown output column(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
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
  picked <- candidates_select_rows(store, index, threshold, threads, strategy)
  fetch <- unique(c(out_columns, if (identical(order, "reconstructed")) "p_value"))
  build <- compressor_normalize_build(store$manifest$genome_build %||% "GRCh38")
  # p_value is attached from the selection step; Z is needed only when the
  # selection (the flag) did not already provide the reconstructed p.
  wanted <- setdiff(fetch, c("p_value", "key"))
  if ("key" %in% fetch) {
    wanted <- unique(c(wanted, "chromosome", "base_pair_location",
                       "other_allele", "effect_allele"))
  }
  if ("p_value" %in% fetch && anyNA(picked$p)) wanted <- unique(c(wanted, "z"))
  list(store = store, index = index, picked = picked, threshold = threshold,
       out_columns = out_columns, fetch = fetch, order = order, build = build,
       range = candidates_region_range(region, build), wanted = wanted)
}

# Stage 2: decode keys and values for the candidate rows only, attach p, key
# and exact rank, order and project.  `keys` is an optional shared key slice
# (see candidates_fetch()).
candidates_finish <- function(ctx, threads, keys = NULL) {
  store <- ctx$store; picked <- ctx$picked; rows <- picked$rows
  fetch <- ctx$fetch; order <- ctx$order; wanted <- ctx$wanted
  data <- if (length(rows)) {
    candidates_fetch(store, ctx$index, rows, ctx$range, ctx$build, wanted,
                     threads, keys = keys)
  } else NULL
  if (is.null(data)) {
    data <- pcodec_native_empty_result(unique(c(wanted, "z")))
    data$p_value <- numeric()
  }
  if (nrow(data) && "p_value" %in% fetch) {
    # The decoder's erfc p (the one reconstruction shared with full and
    # selective reads, see pcodec_native_p_from_codes()). Flag-selected
    # central-bin rows (rare: only for flag thresholds above ~4.8e-4) fall back
    # to erfc of the decoded Z, which can differ from the decoder table by an
    # ulp.
    p <- picked$p[match(data$row, rows)]
    if (anyNA(p)) {
      miss <- is.na(p)
      p[miss] <- candidates_native_p(candidates_semantic(store),
                                     exception_z = data$z[miss])
    }
    data$p_value <- p
  }
  if ("key" %in% fetch) {
    data$key <- if (nrow(data)) {
      compressor_variant_key(data$chromosome, data$base_pair_location,
                             data$other_allele, data$effect_allele,
                             build = ctx$build)
    } else character()
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
  keep <- c("row", ctx$out_columns, if (identical(order, "exact")) "exact_rank")
  data <- data[keep]
  row.names(data) <- NULL
  attr(data, "candidate_strategy") <- picked$strategy
  attr(data, "candidate_threshold") <- ctx$threshold
  attr(data, "candidate_order") <- switch(order,
    none = "native_row", reconstructed = "reconstructed_p_approximate",
    exact = "exact_rank")
  data
}

#' Read threshold candidates from several native Pcodec stores
#'
#' Applies [read_candidates()] to each store. Stores that share a variant
#' panel (equal identity signature) decode the position/substitution key
#' blocks once, for the union of their candidate rows, when identity columns
#' (or `key`, or a region) are requested.
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
                                  threads = 1L, bind = FALSE,
                                  strategy = c("auto", "exceptions", "z_stream",
                                               "pvalue_flag", "exact_order")) {
  strategy <- match.arg(strategy)
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
  guard <- function(expr) tryCatch(expr, error = function(e) e)
  # Stage 1 (parallel over stores): select candidate rows.
  ctxs <- pcodec_parallel_lapply(seq_along(stores), function(i) {
    guard(candidates_prepare(stores[[i]], thresholds[[i]], region, columns,
                             order, inner, strategy))
  }, threads = threads)
  # Stage 2: stores sharing a variant panel (equal identity signature) decode
  # the position/substitution blocks once, for the union of their candidate
  # rows (the candidate analogue of read_sumstats_batch() identity reuse).
  ok <- !vapply(ctxs, inherits, logical(1), "error")
  shared <- vector("list", length(stores))
  need_key <- function(ctx) {
    id <- c("global_position", "substitution", "chromosome", "base_pair_location",
            "reference_allele", "alternate_allele", "effect_allele", "other_allele")
    !is.null(ctx$range) || any(id %in% ctx$wanted)
  }
  if (sum(ok) > 1L) {
    sig <- rep(NA_character_, length(stores))
    sig[ok] <- vapply(ctxs[ok], function(ctx) {
      tryCatch(pcodec_identity_signature(ctx$store), error = function(e) NA_character_)
    }, character(1))
    sig[ok & !vapply(seq_along(ctxs), function(i) ok[i] && need_key(ctxs[[i]]),
                     logical(1))] <- NA_character_
    for (g in unique(stats::na.omit(sig))) {
      members <- which(sig %in% g)
      if (length(members) < 2L) next
      union_rows <- sort(unique(unlist(lapply(ctxs[members], function(ctx) {
        ctx$picked$rows
      }), use.names = FALSE)))
      first <- ctxs[[members[1L]]]
      keys <- if (length(union_rows)) guard(candidates_fetch_keys(
        first$store, first$index, union_rows, first$range, threads)) else NULL
      if (inherits(keys, "error")) next
      .pcodec_batch_trace$identity_resolutions <-
        .pcodec_batch_trace$identity_resolutions + 1L
      for (i in members) shared[[i]] <- keys
    }
  }
  pieces <- pcodec_parallel_lapply(seq_along(stores), function(i) {
    if (!ok[[i]]) return(ctxs[[i]])
    guard(candidates_finish(ctxs[[i]], inner, keys = shared[[i]]))
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

candidates_stream_reader <- function(store) {
  list(
    open = function(spec) file(file.path(store$path, spec$file), open = "rb"),
    read_at = function(con, location) {
      seek(con, where = as.numeric(location$offset), origin = "start")
      blob <- readBin(con, raw(), n = as.integer(location$length), endian = "little")
      if (length(blob) != as.integer(location$length)) {
        stop("native Pcodec stream is truncated", call. = FALSE)
      }
      blob
    })
}

# Key (position + substitution) blocks touched by `rows`, restricted to the
# region range.  Returns list(row, position, substitution) sorted by row, or
# NULL when nothing survives.  Only key blocks that contain a requested row (and
# overlap the region) are decoded.
candidates_fetch_keys <- function(store, index, rows, range, threads) {
  key_blocks <- pcodec_native_index_blocks(index, "key")
  num <- function(blocks, field) vapply(blocks, function(b) as.numeric(b[[field]]),
                                        numeric(1))
  key_stops <- num(key_blocks, "row_stop")
  key_starts <- num(key_blocks, "row_start")
  key_first <- num(key_blocks, "first_position")
  key_last <- num(key_blocks, "last_position")
  delta <- identical(index$position_encoding, "delta_u32_within_block")
  key_id <- findInterval(rows, key_stops) + 1L
  todo <- unique(key_id)
  if (!is.null(range)) todo <- todo[key_last[todo] >= range[1] & key_first[todo] <= range[2]]
  io <- candidates_stream_reader(store)
  open_stream <- io$open
  read_at <- io$read_at
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
  keys <- unlist(keys, recursive = FALSE, use.names = FALSE)
  keys <- keys[!vapply(keys, is.null, logical(1))]
  if (!length(keys)) return(NULL)
  sel_row <- unlist(lapply(keys, `[[`, "row"), use.names = FALSE)
  o <- order(sel_row)
  list(row = as.integer(sel_row[o]),
       position = unlist(lapply(keys, `[[`, "position"), use.names = FALSE)[o],
       substitution = unlist(lapply(keys, `[[`, "substitution"), use.names = FALSE)[o])
}

# `keys` (optional): precomputed candidates_fetch_keys() result for a superset
# of `rows` under the same region (same-panel reuse); sliced, not re-decoded.
candidates_fetch <- function(store, index, rows, range, build, wanted,
                             threads, keys = NULL) {
  manifest <- store$manifest
  semantic <- manifest$semantic_codec %||% list()
  value_blocks <- pcodec_native_index_blocks(index, "value")
  num <- function(blocks, field) vapply(blocks, function(b) as.numeric(b[[field]]),
                                        numeric(1))
  identity <- c("global_position", "substitution", "chromosome",
                "base_pair_location", "reference_allele", "alternate_allele",
                "effect_allele", "other_allele")
  need_identity <- !is.null(range) || any(identity %in% wanted)
  io <- candidates_stream_reader(store)
  open_stream <- io$open
  read_at <- io$read_at
  if (need_identity) {
    if (is.null(keys)) {
      keys <- candidates_fetch_keys(store, index, rows, range, threads)
    } else {
      at <- match(rows, keys$row)
      at <- at[!is.na(at)]
      keys <- if (length(at)) lapply(keys, function(v) v[at]) else NULL
    }
    if (is.null(keys)) return(NULL)
    sel_row <- keys$row
    sel_position <- keys$position
    sel_substitution <- keys$substitution
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

#' Feature flags of this CompreSSoR build
#'
#' Lets dependent packages feature-detect capabilities without version
#' parsing. `"candidates_one_pass"` means [read_candidates()] returns values,
#' `key`, `p_value` (bit-identical to [read_sumstats()]) and exact ranks for
#' the candidate rows in a single pass, and [read_candidates_batch()] reuses
#' same-panel identity.
#'
#' @return A character vector of capability names.
#' @export
compressor_capabilities <- function() {
  c("candidates_one_pass", "candidate_key_column", "p_value_shared_reconstruction")
}
