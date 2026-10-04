# Threshold/region candidate extraction for native Pcodec stores.
#
# Candidate membership is defined as the full-read convention
# `p_value <= threshold`, where p is the p-value reconstructed from the lossy
# stored Z (the native decoder's erfc table for central bins and erfc of the
# exact float32 Z for exceptions).  Three strategies reach that set at
# different cost; see read_candidates().

candidates_semantic <- function(store) {
  pcodec_native_semantic_params(store$manifest$semantic_codec)
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
    semantic$se_range[1], semantic$se_range[2], NULL,
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
    candidates_check_flag_count(store, length(rows), "flag stream")
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

# The manifest records the number of flagged rows at write time
# (`domains$pvalue_flag$hit_rows`).  A pvalue_flag selection without a region
# must return exactly that many rows; anything else is partial data, so stop
# rather than return it.  Stores whose manifest lacks the count are not checked.
candidates_check_flag_count <- function(store, n, what) {
  expected <- ((store$manifest$domains %||% list())$pvalue_flag %||% list())$hit_rows
  if (is.null(expected) || length(expected) != 1L || is.na(expected)) return(invisible(TRUE))
  if (!identical(as.numeric(n), as.numeric(expected))) {
    stop("candidate row count mismatch in store '", store$path, "': ", what,
         " gave ", n, " pvalue_flag rows but the manifest records ", expected,
         " flagged rows; refusing to return partial data", call. = FALSE)
  }
  invisible(TRUE)
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
    }, threads = threads, labels = store$path, what = "p-value flag read")
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
  # Candidate reads use only the cached block matrices, never the parsed
  # index; `index` stays NULL in the context (cheap to return from a worker).
  index <- NULL
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
  # Without a region every selected row must come back, in row order.
  if (is.null(ctx$range) && !identical(as.integer(data$row), as.integer(rows))) {
    stop("candidate rows lost in decode for store '", store$path, "': selected ",
         length(rows), " rows, decoded ", nrow(data), call. = FALSE)
  }
  if (identical(picked$strategy, "pvalue_flag") && is.null(ctx$range)) {
    candidates_check_flag_count(store, nrow(data), "the candidate read")
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
#' Applies [read_candidates()] to each store, stores in parallel, and returns
#' exactly the per-store results. With
#' `options(CompreSSoR.candidates_share_keys = TRUE)`, stores that share a
#' variant panel (equal identity signature, each verified against its own
#' manifest) decode the position/substitution key blocks once, for the union
#' of their candidate rows, when identity columns (or `key`, or a region) are
#' requested; by default each store decodes its own candidate keys, which with
#' the native selective reader is faster.
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
  labels <- vapply(seq_along(stores), function(i) {
    candidates_store_label(stores[[i]], names(stores)[i])
  }, character(1))
  need_key <- function(ctx) {
    id <- c("global_position", "substitution", "chromosome", "base_pair_location",
            "reference_allele", "alternate_allele", "effect_allele", "other_allele")
    !is.null(ctx$range) || any(id %in% ctx$wanted)
  }
  share <- length(stores) > 1L && candidates_share_keys()
  if (!share) {
    # One pass per store: select and fetch in the same worker.
    pieces <- pcodec_parallel_lapply(seq_along(stores), function(i) {
      guard(candidates_finish(candidates_prepare(stores[[i]], thresholds[[i]], region,
                                                 columns, order, inner, strategy),
                              inner))
    }, threads = threads, labels = labels, what = "candidate read")
    return(candidates_batch_result(pieces, stores, labels, thresholds, bind))
  }
  # Stage 1 (parallel over stores): select candidate rows; each worker also
  # returns the store's compact block matrices and, when keys may be shared,
  # its identity signature, so the parent parses no index. A store's own
  # error is carried as a condition object and reported below with its
  # label; a dead worker (NULL) is rejected by pcodec_parallel_lapply().
  ctxs <- pcodec_parallel_lapply(seq_along(stores), function(i) {
    guard({
      ctx <- candidates_prepare(stores[[i]], thresholds[[i]], region, columns,
                                order, inner, strategy)
      ctx$mats <- pcodec_native_select_matrices(ctx$store)
      ctx$signature <- if (share && need_key(ctx)) {
        tryCatch(pcodec_identity_signature(ctx$store), error = function(e) NA_character_)
      } else NA_character_
      ctx
    })
  }, threads = threads, labels = labels, what = "candidate selection")
  ctxs <- lapply(ctxs, function(x) {
    if (is.list(x) && !inherits(x, "error") && is.list(x$picked)) return(x)
    if (inherits(x, "error")) return(x)
    simpleError("candidate selection worker returned no context")
  })
  ok <- !vapply(ctxs, inherits, logical(1), "error")
  if (length(stores) > 1L && threads > 1L) {
    for (ctx in ctxs[ok]) pcodec_native_register_matrices(ctx$store, ctx$mats)
  }
  # Stage 2: stores sharing a variant panel (equal identity signature) decode
  # the position/substitution blocks once, for the union of their candidate
  # rows (the candidate analogue of read_sumstats_batch() identity reuse).
  # `shared` must keep one slot per store: assign with `[<-` and list(), never
  # `shared[[i]] <- NULL`, which deletes slot i and shifts every later store's
  # keys onto the wrong store (the cause of silently dropped candidates).
  shared <- vector("list", length(stores))
  if (share && sum(ok) > 1L) {
    sig <- rep(NA_character_, length(stores))
    sig[ok] <- vapply(ctxs[ok], function(ctx) ctx$signature %||% NA_character_,
                      character(1))
    grouped <- which(!is.na(sig) & sig %in% sig[duplicated(sig) & !is.na(sig)])
    # Every store about to receive another store's decoded keys must first
    # match its own manifest checksums (a corrupt store errors here).
    pcodec_verify_identity_files(lapply(ctxs[grouped], `[[`, "store"),
                                 threads = threads)
    groups <- lapply(unique(sig[grouped]), function(g) which(sig %in% g))
    # Union key decode per group; groups run in parallel, each on
    # threads %/% groups native threads.
    group_threads <- max(1L, threads %/% max(1L, length(groups)))
    group_keys <- pcodec_parallel_lapply(groups, function(members) {
      union_rows <- sort(unique(unlist(lapply(ctxs[members], function(ctx) {
        ctx$picked$rows
      }), use.names = FALSE)))
      # No candidate rows in the group (or none in the region): nothing to
      # share; each member's own finish step handles its empty result.
      if (!length(union_rows)) return(PCODEC_PARALLEL_EMPTY)
      first <- ctxs[[members[1L]]]
      keys <- guard(candidates_fetch_keys(
        first$store, first$index, union_rows, first$range, group_threads))
      if (is.null(keys) || inherits(keys, "error")) PCODEC_PARALLEL_EMPTY else keys
    }, threads = if (length(groups) > 1L) threads else 1L,
    what = "shared candidate key decode")
    for (g in seq_along(groups)) {
      if (pcodec_parallel_is_empty(group_keys[[g]])) next
      .pcodec_batch_trace$identity_resolutions <-
        .pcodec_batch_trace$identity_resolutions + 1L
      shared[groups[[g]]] <- list(group_keys[[g]])
    }
  }
  stopifnot(length(shared) == length(stores))
  pieces <- pcodec_parallel_lapply(seq_along(stores), function(i) {
    if (!ok[[i]]) return(ctxs[[i]])
    guard(candidates_finish(ctxs[[i]], inner, keys = shared[[i]]))
  }, threads = threads, labels = labels, what = "candidate decode")
  candidates_batch_result(pieces, stores, labels, thresholds, bind)
}

# Check and assemble read_candidates_batch() results (one per store).
candidates_batch_result <- function(pieces, stores, labels, thresholds, bind) {
  failed <- vapply(pieces, function(x) !is.data.frame(x), logical(1))
  if (any(failed)) {
    first <- which(failed)[1L]
    msg <- if (inherits(pieces[[first]], "error")) conditionMessage(pieces[[first]])
           else "worker failed"
    stop("failed to read candidates from store ", first, " ('", labels[first],
         "'): ", msg, call. = FALSE)
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

# Whether read_candidates_batch() decodes the keys of same-panel stores once
# (option CompreSSoR.candidates_share_keys, default FALSE). With the native
# selective reader each store's candidate keys cost milliseconds, so one
# pass per store (no grouping, verification or second fork round) is faster.
candidates_share_keys <- function() {
  isTRUE(getOption("CompreSSoR.candidates_share_keys", FALSE))
}

# A store's name in a batch (its list name, else its path) for error messages.
candidates_store_label <- function(store, name = NULL) {
  if (!is.null(name) && !is.na(name) && nzchar(name)) return(name)
  if (inherits(store, "compressor_store")) return(store$path)
  if (is.character(store) && length(store) == 1L && !is.na(store)) return(store)
  "<store>"
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
  # An empty overlap is an inverted range, which every caller treats as
  # matching no rows.
  pcodec_native_region_range(region, build) %||% c(1, 0)
}


# One native selective read (mode 0: sorted zero-based row IDs) of `streams`
# for `rows`: only the key/value blocks holding those rows are read and
# decoded, every payload file is opened once, and blocks are decoded on up
# to `threads` native threads.
candidates_native_select <- function(store, index, rows, streams, identity,
                                     threads) {
  mats <- pcodec_native_select_matrices(store, index)
  threads <- pcodec_validate_threads(threads)
  check_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(check_limit) && check_limit != "false") threads <- min(threads, 2L)
  .Call("compressor_read_pcodec_native_select", mats$files, mats$position,
        mats$substitution, mats$z, mats$eaf, mats$se, mats$exceptions,
        as.numeric(store$manifest$n_rows %||% store$manifest$rows),
        as.character(streams), mats$exception_codec,
        as.integer(threads), 0L, as.numeric(rows), isTRUE(identity),
        isTRUE(mats$delta), PACKAGE = "CompreSSoR")
}

# Keys (position + substitution) of `rows` (sorted, zero-based), restricted
# to the region range. Returns list(row, position, substitution) sorted by
# row, or NULL when nothing survives. Only the key blocks that hold a
# requested row are decoded.
candidates_fetch_keys <- function(store, index, rows, range, threads) {
  if (!is.null(range) && length(rows)) {
    # Skip rows whose key block lies outside the region (block anchors from
    # the cached matrices: row_stop, first and last position).
    pos <- pcodec_native_select_matrices(store, index)$position
    block <- findInterval(rows, pos[, 5L]) + 1L
    rows <- rows[pos[block, 7L] >= range[1] & pos[block, 6L] <= range[2]]
  }
  if (!length(rows)) return(NULL)
  res <- candidates_native_select(store, index, rows,
                                  c("position", "substitution"), TRUE, threads)
  keep <- if (is.null(range)) rep.int(TRUE, length(res$rows)) else
    res$position >= range[1] & res$position <= range[2]
  if (!any(keep)) return(NULL)
  list(row = as.integer(res$rows[keep]), position = as.numeric(res$position[keep]),
       substitution = as.integer(res$substitution[keep]))
}

# `keys` (optional): precomputed candidates_fetch_keys() result for a superset
# of `rows` under the same region (same-panel reuse); sliced, not re-decoded.
# Values come from one native selective read of the surviving rows and are
# decoded by the native decoder, so they are bit-identical to a full read.
candidates_fetch <- function(store, index, rows, range, build, wanted,
                             threads, keys = NULL) {
  semantic <- store$manifest$semantic_codec %||% list()
  identity <- c("global_position", "substitution", "chromosome",
                "base_pair_location", "reference_allele", "alternate_allele",
                "effect_allele", "other_allele")
  need_identity <- !is.null(range) || any(identity %in% wanted)
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
  } else {
    sel_row <- as.integer(rows)
  }

  out <- list(row = sel_row)
  if (any(identity %in% wanted)) {
    cols <- pcodec_native_key_columns(keys$position, keys$substitution, build = build)
    cols$global_position <- keys$position
    cols$substitution <- keys$substitution
    for (nm in intersect(identity, wanted)) out[[nm]] <- cols[[nm]]
  }
  need_z <- any(c("z", "beta") %in% wanted)
  need_se <- any(c("standard_error", "beta") %in% wanted)
  need_eaf <- "effect_allele_frequency" %in% wanted || need_se
  needed <- c(if (need_z) "z", if (need_eaf) "eaf", if (need_se) "se")
  if (length(needed) && length(sel_row)) {
    res <- candidates_native_select(store, index, sel_row, needed, FALSE, threads)
    if (!identical(as.integer(res$rows), as.integer(sel_row))) {
      stop("candidate value decode returned other rows than requested in store '",
           store$path, "'", call. = FALSE)
    }
    decoded <- pcodec_native_decode_rows(
      semantic, sel_row, list(z = res$z, se = res$se, eaf = res$eaf),
      list(index = res$exc_index, z = res$exc_z, log2se = res$exc_log2se,
           eaf = res$exc_eaf, flags = res$exc_flags),
      want_beta = "beta" %in% wanted)
    if (need_z) out$z <- decoded$z
    if (need_se) out$standard_error <- decoded$standard_error
    if (need_eaf) out$effect_allele_frequency <- decoded$effect_allele_frequency
    if ("beta" %in% wanted) out$beta <- decoded$beta
  }
  out <- out[c("row", intersect(c(identity, "z", "beta", "standard_error",
                                  "effect_allele_frequency"), names(out)))]
  structure(out, class = "data.frame", row.names = .set_row_names(length(sel_row)))
}

#' Feature flags of this CompreSSoR build
#'
#' Lets dependent packages feature-detect capabilities without version
#' parsing. `"candidates_one_pass"` means [read_candidates()] returns values,
#' `key`, `p_value` (bit-identical to [read_sumstats()]) and exact ranks for
#' the candidate rows in a single pass, and [read_candidates_batch()] can
#' reuse same-panel identity. `"candidates_batch_rows_checked"` means
#' [read_candidates_batch()] returns exactly the per-store [read_candidates()]
#' result for every batch composition (the shared-identity row loss in
#' 0.7.0 is fixed) and both readers stop, rather than return partial data,
#' when a decode loses selected rows or a `pvalue_flag` read disagrees with
#' the manifest's flagged-row count. `"reads_bit_identical"` means every
#' reader (full, row-ID, key, region, candidate and batched reads) decodes a
#' row with the same native decoder, so values are bit-identical to the full
#' read. `"integrity_verified"` means `validate_compressor(full = TRUE)`
#' checks every payload checksum, the native index is checked whenever it is
#' parsed, and batched identity sharing verifies each member's identity
#' streams.
#'
#' @return A character vector of capability names.
#' @export
compressor_capabilities <- function() {
  c("candidates_one_pass", "candidate_key_column", "p_value_shared_reconstruction",
    "candidates_batch_rows_checked", "reads_bit_identical", "integrity_verified")
}
