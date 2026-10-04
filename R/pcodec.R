## Native Pcodec store interface.
##
## The historical Python-backed implementation lives under archive/ and is not
## part of the installed package. New stores are native 0.4 stores only.

pcodec_manifest_checksum_path <- function(path) {
  file.path(dirname(path), "manifest.sha256")
}

pcodec_payload_sha256 <- function(files) {
  if (!length(files)) return(digest::digest("", algo = "sha256", serialize = FALSE))
  file_names <- names(files)
  if (is.null(file_names) || anyNA(file_names) || any(!nzchar(file_names))) {
    stop("Pcodec payload hash inputs must be named", call. = FALSE)
  }
  entries <- vapply(sort(file_names), function(name) {
    item <- files[[name]]
    paste(name, item$bytes, item$sha256, sep = "\t")
  }, character(1))
  digest::digest(paste(entries, collapse = "\n"), algo = "sha256", serialize = FALSE)
}

# SHA-256 hex digests of files (NA for a file that cannot be read), hashed
# natively on up to `threads` threads; equal to digest::digest(path, algo =
# "sha256", file = TRUE).
pcodec_sha256_files <- function(paths, threads = 1L) {
  paths <- as.character(paths)
  if (!length(paths)) return(character())
  if (is.loaded("compressor_sha256_files", PACKAGE = "CompreSSoR")) {
    return(.Call("compressor_sha256_files", paths, as.integer(threads),
                 PACKAGE = "CompreSSoR"))
  }
  vapply(paths, function(path) {
    if (!file.exists(path) || dir.exists(path)) return(NA_character_)
    digest::digest(path, algo = "sha256", file = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

# Compare files of a store with its manifest integrity record (byte count and
# sha256). `names` are relative payload paths (default: every recorded file).
# With `payload = TRUE` the aggregate payload_sha256 is recomputed from the
# per-file record and every file named in manifest$files must have a record.
# Returns a character vector of problems (empty when everything matches).
pcodec_integrity_problems <- function(store, names = NULL, threads = 1L,
                                      payload = is.null(names)) {
  force(payload)  # its default depends on `names`, which is reassigned below
  integrity <- store$manifest$integrity
  files <- integrity$files
  if (!is.list(files) || !length(files)) {
    return("the manifest has no per-file integrity record")
  }
  if (is.null(names)) names <- names(files)
  problems <- character()
  unknown <- setdiff(names, names(files))
  if (length(unknown)) {
    problems <- c(problems, paste0(unknown, ": no integrity record in the manifest"))
  }
  names <- intersect(names, names(files))
  paths <- file.path(store$path, names)
  info <- file.info(paths, extra_cols = FALSE)
  size <- info$size
  missing <- is.na(size) | info$isdir %in% TRUE
  if (any(missing)) problems <- c(problems, paste0(names[missing], ": missing"))
  recorded_bytes <- vapply(files[names], function(item) {
    as.numeric(item$bytes %||% NA_real_)
  }, numeric(1))
  recorded_sha <- vapply(files[names], function(item) {
    tolower(as.character(item$sha256 %||% NA_character_))
  }, character(1))
  wrong_size <- !missing & !is.na(recorded_bytes) & size != recorded_bytes
  if (any(wrong_size)) {
    problems <- c(problems, sprintf(
      "%s: %s bytes on disk but the manifest records %s (truncated or modified)",
      names[wrong_size], format(size[wrong_size], scientific = FALSE),
      format(recorded_bytes[wrong_size], scientific = FALSE)))
  }
  check <- !missing & !wrong_size
  if (any(check)) {
    observed <- pcodec_sha256_files(paths[check], threads = threads)
    bad <- is.na(observed) | is.na(recorded_sha[check]) |
      tolower(observed) != recorded_sha[check]
    if (any(bad)) {
      problems <- c(problems, paste0(names[check][bad],
                                     ": sha256 mismatch (corrupt or modified file)"))
    }
  }
  if (isTRUE(payload)) {
    if (!identical(tolower(as.character(integrity$payload_sha256 %||% "")),
                   pcodec_payload_sha256(files))) {
      problems <- c(problems,
                    "payload_sha256 does not match the per-file integrity record")
    }
    listed <- unname(unlist(store$manifest$files))
    unrecorded <- setdiff(listed, names(files))
    if (length(unrecorded)) {
      problems <- c(problems, paste0(unrecorded,
                                     ": payload file has no integrity record"))
    }
  }
  problems
}

# Batch identity sharing serves one store's decoded keys to every store with
# the same identity signature, and the signature is built from manifest
# hashes. Before a store takes part, its position and substitution streams
# are hashed and compared with its own manifest record, so a corrupt or
# half-copied store fails instead of being masked by another store's keys.
# (The index, which holds the block anchors, is verified whenever it is
# parsed.) A verified result is cached per (path, size, mtime, recorded hash).
.pcodec_identity_verified <- new.env(parent = emptyenv())

pcodec_identity_file_names <- function(store) {
  files <- store$manifest$files
  c(files$position, files$substitution)
}

pcodec_identity_verify_stamp <- function(store) {
  names <- pcodec_identity_file_names(store)
  if (length(names) != 2L) return(NULL)
  paths <- file.path(store$path, names)
  info <- file.info(paths, extra_cols = FALSE)
  recorded <- vapply(names, function(name) {
    as.character(store$manifest$integrity$files[[name]]$sha256 %||% NA_character_)
  }, character(1))
  paste(paths, info$size, format(as.numeric(info$mtime), digits = 17), recorded,
        sep = "|", collapse = ";")
}

# Verify the identity streams of `stores` (a list of open stores). Stores
# whose manifests record the same identity hashes form a set: one member
# (preferably one verified earlier) is hashed and compared with its record,
# and every other member's files are compared byte for byte with that
# verified member's files -- equal bytes then match the same recorded hash,
# at memory-compare rather than SHA-256 cost. Stops on the first store that
# fails, naming it. A verified result is cached per file stamp.
pcodec_verify_identity_files <- function(stores, threads = 1L) {
  if (!length(stores)) return(invisible(TRUE))
  stamps <- lapply(stores, pcodec_identity_verify_stamp)
  for (i in seq_along(stores)) {
    if (is.null(stamps[[i]])) {
      stop("store '", stores[[i]]$path, "' has no identity stream record; ",
           "it cannot take part in shared identity decoding", call. = FALSE)
    }
  }
  verified <- vapply(seq_along(stores), function(i) {
    identical(.pcodec_identity_verified[[stores[[i]]$path]], stamps[[i]])
  }, logical(1))
  if (all(verified)) return(invisible(TRUE))
  names_of <- lapply(stores, pcodec_identity_file_names)
  recorded <- vapply(seq_along(stores), function(i) {
    paste(vapply(names_of[[i]], function(name) {
      tolower(as.character(stores[[i]]$manifest$integrity$files[[name]]$sha256 %||% NA_character_))
    }, character(1)), collapse = "|")
  }, character(1))
  fail <- function(store, files) {
    stop("store '", store$path, "' failed identity verification: ",
         paste(files, collapse = ", "),
         " do not match the store's manifest checksums (corrupt, truncated ",
         "or half-copied store); refusing to share decoded identity with it",
         call. = FALSE)
  }
  paths_of <- function(i) file.path(stores[[i]]$path, names_of[[i]])
  # One reference per recorded-hash set: a verified member, else the first,
  # hashed now (all such references in one threaded call).
  sets <- split(seq_along(stores), recorded)
  reference <- vapply(sets, function(members) {
    hit <- members[verified[members]]
    if (length(hit)) hit[1L] else members[1L]
  }, integer(1))
  to_hash <- reference[!verified[reference]]
  if (length(to_hash)) {
    digests <- pcodec_sha256_files(unlist(lapply(to_hash, paths_of)), threads = threads)
    at <- 0L
    for (i in to_hash) {
      n <- names_of[[i]]
      observed <- tolower(digests[at + seq_along(n)])
      at <- at + length(n)
      expected <- strsplit(recorded[[i]], "|", fixed = TRUE)[[1L]]
      bad <- is.na(observed) | observed != expected
      if (any(bad)) fail(stores[[i]], n[bad])
      .pcodec_identity_verified[[stores[[i]]$path]] <- stamps[[i]]
      verified[i] <- TRUE
    }
  }
  # Every other unverified member: byte-compare with its set's reference.
  todo <- which(!verified)
  if (length(todo)) {
    ref_of <- reference[match(recorded[todo], names(sets))]
    a <- unlist(lapply(todo, paths_of))
    b <- unlist(lapply(ref_of, paths_of))
    equal <- pcodec_files_equal(a, b, threads = threads)
    at <- 0L
    for (t in seq_along(todo)) {
      i <- todo[t]
      n <- names_of[[i]]
      same <- equal[at + seq_along(n)]
      at <- at + length(n)
      if (!all(same)) fail(stores[[i]], n[!same])
      .pcodec_identity_verified[[stores[[i]]$path]] <- stamps[[i]]
    }
  }
  invisible(TRUE)
}

# Pairwise byte equality of files (natively, `threads` pairs at a time).
pcodec_files_equal <- function(a, b, threads = 1L) {
  if (!length(a)) return(logical())
  if (is.loaded("compressor_files_equal", PACKAGE = "CompreSSoR")) {
    return(.Call("compressor_files_equal", as.character(a), as.character(b),
                 as.integer(threads), PACKAGE = "CompreSSoR"))
  }
  mapply(function(x, y) {
    sx <- file.size(x); sy <- file.size(y)
    !is.na(sx) && !is.na(sy) && sx == sy &&
      identical(readBin(x, raw(), sx), readBin(y, raw(), sy))
  }, a, b, USE.NAMES = FALSE)
}

pcodec_canonical_manifest_sha256 <- function(manifest) {
  strip_observational <- function(value) {
    if (!is.list(value)) return(value)
    names_value <- names(value)
    if (!is.null(names_value)) {
      value <- value[!names_value %in% c("timings", "elapsed_seconds",
                                         "read_elapsed_seconds",
                                         "projection_elapsed_seconds",
                                         "build_info",
                                         # thread counts: observational, so
                                         # the hash is thread-independent
                                         "threads", "threads_requested",
                                         "writer", "requested_workers",
                                         "effective_workers")]
    }
    lapply(value, strip_observational)
  }
  canonical <- strip_observational(manifest)
  canonical$created_utc <- NULL
  # Wall-clock timings and read-duration observations are metadata, not part
  # of the deterministic store contract or its canonical manifest identity.
  if (!is.null(canonical$integrity)) {
    canonical$integrity$canonical_sha256 <- NULL
  }
  payload <- jsonlite::toJSON(canonical, auto_unbox = TRUE, pretty = FALSE,
                              null = "null", digits = 17)
  digest::digest(payload, algo = "sha256", serialize = FALSE)
}

# The form a manifest takes after write_manifest() and read_manifest(), as far
# as serialization is concerned: jsonlite pretty-prints an atomic vector of
# length other than one on a single line but the list it reads back one
# element per line, so atomic vectors become unnamed lists of scalars. Length-
# one values and lists of scalars serialize the same either way. Anything
# else (data frames, factors, classed or dimensioned values, vectors holding
# NA/NaN/Inf) takes a real JSON round trip, so write_manifest() of the result
# produces exactly the bytes of the former write/read_json/re-write sequence.
manifest_json_normalise <- function(x) {
  if (is.null(x)) return(x)
  attribute_names <- names(attributes(x))
  plain <- !length(setdiff(attribute_names, "names"))
  if (is.list(x)) {
    if (!plain) return(manifest_json_round_trip(x))
    return(lapply(x, manifest_json_normalise))
  }
  if (is.atomic(x) && plain) {
    if (length(x) == 1L) return(x)
    if (anyNA(x) || (is.double(x) && any(is.infinite(x)))) {
      return(manifest_json_round_trip(x))
    }
    return(as.list(unname(x)))
  }
  manifest_json_round_trip(x)
}

manifest_json_round_trip <- function(x) {
  jsonlite::parse_json(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = 17),
                       simplifyVector = FALSE)
}

# The sealing step of seal_pcodec_manifest() on an in-memory manifest.
pcodec_seal_manifest_value <- function(manifest) {
  if (!is.null(manifest$integrity$files) &&
      !is.null(manifest$integrity$payload_sha256)) {
    manifest$integrity$canonical_sha256 <- pcodec_canonical_manifest_sha256(manifest)
  }
  manifest
}

# Write an already-sealed manifest and its manifest.sha256 record. With
# `atomic = TRUE` (a manifest inside a committed store) both are written to
# temporary files beside their targets and renamed into place, manifest first
# and checksum immediately after, so neither file is ever seen truncated.
write_pcodec_manifest <- function(manifest, path, atomic = FALSE) {
  checksum_path <- pcodec_manifest_checksum_path(path)
  if (!isTRUE(atomic)) {
    write_manifest(manifest, path)
    checksum <- digest::digest(path, algo = "sha256", file = TRUE)
    writeLines(checksum, checksum_path, useBytes = TRUE)
    return(invisible(checksum))
  }
  manifest_tmp <- tempfile(".manifest-", tmpdir = dirname(path), fileext = ".json")
  checksum_tmp <- tempfile(".manifest-", tmpdir = dirname(path), fileext = ".sha256")
  on.exit(unlink(c(manifest_tmp, checksum_tmp), force = TRUE), add = TRUE)
  write_manifest(manifest, manifest_tmp)
  checksum <- digest::digest(manifest_tmp, algo = "sha256", file = TRUE)
  writeLines(checksum, checksum_tmp, useBytes = TRUE)
  if (!file.rename(manifest_tmp, path) || !file.rename(checksum_tmp, checksum_path)) {
    stop("could not replace the store manifest at ", path, call. = FALSE)
  }
  invisible(checksum)
}

seal_pcodec_manifest <- function(path) {
  manifest <- read_manifest(path)
  if (!is.null(manifest$integrity$files) &&
      !is.null(manifest$integrity$payload_sha256)) {
    manifest$integrity$canonical_sha256 <- pcodec_canonical_manifest_sha256(manifest)
    write_manifest(manifest, path)
  }
  checksum <- digest::digest(path, algo = "sha256", file = TRUE)
  writeLines(checksum, pcodec_manifest_checksum_path(path), useBytes = TRUE)
  invisible(checksum)
}

# `expected` and `observed` may be supplied by a caller that already read the
# checksum record and hashed the manifest (open_compressor does, so the file is
# hashed once per open); standalone callers get the original behaviour.
verify_pcodec_manifest <- function(path, expected = NULL, observed = NULL) {
  if (is.null(expected)) {
    checksum_path <- pcodec_manifest_checksum_path(path)
    if (!file.exists(checksum_path)) {
      stop("Pcodec store is missing manifest.sha256", call. = FALSE)
    }
    expected <- readLines(checksum_path, warn = FALSE, n = 1L)
  }
  expected <- trimws(expected)
  if (length(expected) != 1L || !grepl("^[0-9a-fA-F]{64}$", expected)) {
    stop("Pcodec manifest checksum record is malformed", call. = FALSE)
  }
  if (is.null(observed)) observed <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!identical(tolower(expected), tolower(observed))) {
    stop("Pcodec manifest checksum mismatch", call. = FALSE)
  }
  invisible(TRUE)
}

pcodec_open_store_cached <- function(store) {
  if (inherits(store, "compressor_store")) return(store)
  open_compressor(store)
}

pcodec_write_store <- function(data, output, metadata = list(), eaf_coverage = NULL,
                               finalize = TRUE) {
  required <- c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error",
                "effect_allele_frequency", "z")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Pcodec input is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!nrow(data)) stop("cannot write an empty Pcodec store", call. = FALSE)
  if (!pcodec_native_available()) {
    stop("CompreSSoR requires its native Pcodec backend; install Rust/Cargo and reinstall the package",
         call. = FALSE)
  }
  pcodec_native_write_store(data, output, metadata = metadata,
                            eaf_coverage = eaf_coverage, finalize = finalize)
}

pcodec_native_projection <- function(out, columns = NULL) {
  source_bytes_read <- attr(out, "source_bytes_read", exact = TRUE)
  if ("row" %in% names(out)) out$row <- NULL
  if (is.null(columns)) {
    attr(out, "source_bytes_read") <- source_bytes_read
    return(out)
  }
  missing <- setdiff(columns, names(out))
  if (length(missing)) {
    stop("requested columns are not present: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  out <- out[columns]
  attr(out, "source_bytes_read") <- source_bytes_read
  out
}

# Zero-based row IDs as integers. Non-integer, missing or non-finite values
# are an error (as.integer() would silently truncate 1.5 to row 1).
pcodec_row_ids <- function(x) {
  if (!is.numeric(x)) stop("row IDs must be numeric", call. = FALSE)
  if (is.integer(x)) {
    if (anyNA(x)) stop("row IDs must not be missing", call. = FALSE)
    return(x)
  }
  if (anyNA(x)) stop("row IDs must not be missing", call. = FALSE)
  if (any(!is.finite(x) | x != floor(x) | abs(x) > .Machine$integer.max)) {
    stop("row IDs must be whole numbers (zero-based integer row IDs)", call. = FALSE)
  }
  as.integer(x)
}

pcodec_validate_threads <- function(threads, label = "threads") {
  if (length(threads) != 1L || !is.numeric(threads) || is.na(threads) ||
      !is.finite(threads) || threads < 1 || threads != floor(threads)) {
    stop(label, " must be one positive integer", call. = FALSE)
  }
  as.integer(threads)
}

pcodec_native_default_threads <- function(region = NULL, variants = NULL,
                                           threads = NULL) {
  if (!is.null(threads)) return(pcodec_validate_threads(threads))
  configured <- getOption("CompreSSoR.pcodec.threads", NULL)
  if (!is.null(configured)) return(pcodec_validate_threads(configured))
  # Whole-file reads benefit from independent stream decoding in parallel.
  # Regional, key and row reads start from one thread here; the native
  # selective reader then raises it with the number of key blocks touched
  # (pcodec_native_select_auto_threads()).
  if (!is.null(region) || !is.null(variants)) 1L else 4L
}

# lapply() over X, forked over up to `threads` workers on Unix-like systems.
# Every result is checked centrally (pcodec_parallel_check()): a worker error
# is re-raised, and a missing result -- mclapply() returns NULL for a forked
# worker that died (OOM killer, signal) -- is an error, never data. FUN must
# therefore not return NULL; return PCODEC_PARALLEL_EMPTY (or another
# non-NULL value) for "nothing". `labels` (one per element of X, e.g. store
# paths) attributes a failure to the element it came from.
PCODEC_PARALLEL_EMPTY <- structure(list(), class = "compressor_parallel_empty")

pcodec_parallel_is_empty <- function(x) inherits(x, "compressor_parallel_empty")

pcodec_parallel_lapply <- function(X, FUN, threads = 1L, labels = NULL,
                                   what = "parallel worker") {
  threads <- pcodec_validate_threads(threads)
  if (length(labels) == 1L && length(X) != 1L) labels <- rep(labels, length(X))
  if (!is.null(labels) && length(labels) != length(X)) {
    stop("internal error: parallel labels must match the work items", call. = FALSE)
  }
  # R CMD check sets this guard to prevent packages from spawning an
  # uncontrolled number of workers. Respect it while retaining the native
  # four-thread default for ordinary whole-file reads.
  check_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(check_limit) && check_limit != "false") threads <- min(threads, 2L)
  serial <- length(X) <= 1L || threads <= 1L ||
    # Forked workers are safe here because each worker opens independent
    # files and calls the standalone Pcodec decoder on private R objects.
    # Windows has no fork backend; retain deterministic serial behaviour.
    .Platform$OS.type == "windows"
  if (serial) {
    out <- if (is.null(labels)) lapply(X, FUN) else {
      lapply(seq_along(X), function(i) {
        tryCatch(FUN(X[[i]]), error = function(e) {
          pcodec_parallel_raise(e, labels[[i]], what)
        })
      })
    }
    return(pcodec_parallel_check(out, length(X), labels, what))
  }
  # Each element's own error is caught in the worker and returned as a
  # marked condition, so it is attributed to its element (mclapply() would
  # mark every element of a failed prescheduled job). Failures are re-raised
  # by pcodec_parallel_check(); drop mclapply()'s generic warnings about them.
  guarded <- function(x) {
    tryCatch(FUN(x), error = function(e) {
      structure(list(condition = e), class = "compressor_parallel_error")
    })
  }
  out <- withCallingHandlers(
    parallel::mclapply(X, guarded, mc.cores = min(threads, length(X)),
                       mc.preschedule = TRUE),
    warning = function(w) {
      if (grepl("encountered errors? in user code|did not deliver",
                conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
  pcodec_parallel_check(out, length(X), labels, what)
}

pcodec_parallel_raise <- function(condition, label, what) {
  is_condition <- inherits(condition, "condition")
  if (is.null(label) && is_condition) stop(condition)
  message <- if (is_condition) conditionMessage(condition) else
    sub("^Error[^:]*: ", "", trimws(paste(as.character(condition), collapse = " ")))
  if (is.null(label)) stop(message, call. = FALSE)
  stop(what, " failed for '", label, "': ", message, call. = FALSE)
}

# Central result check for pcodec_parallel_lapply(): exactly one result per
# work item; a try-error re-raises the worker's condition; NULL (a dead
# forked worker) is an error.
pcodec_parallel_check <- function(out, expected, labels = NULL,
                                  what = "parallel worker") {
  if (length(out) != expected) {
    stop(what, " returned ", length(out), " results for ", expected,
         " work items", call. = FALSE)
  }
  if (!expected) return(out)
  marked <- vapply(out, inherits, logical(1), "compressor_parallel_error")
  if (any(marked)) {
    at <- which(marked)[1L]
    pcodec_parallel_raise(out[[at]]$condition, labels[at], what)
  }
  failed <- vapply(out, inherits, logical(1), "try-error")
  if (any(failed)) {
    at <- which(failed)[1L]
    bad <- out[[at]]
    condition <- attr(bad, "condition")
    pcodec_parallel_raise(if (inherits(condition, "condition")) condition else bad,
                          labels[at], what)
  }
  missing <- vapply(out, is.null, logical(1))
  if (any(missing)) {
    at <- which(missing)[1L]
    stop(what, " for ", if (is.null(labels)) paste("item", at) else
      paste0("'", labels[[at]], "'"),
      " returned no result (the worker process died, for example killed ",
      "for memory); refusing to return partial data", call. = FALSE)
  }
  out
}

pcodec_read_store <- function(store, region = NULL, variants = NULL,
                               columns = NULL, threads = NULL) {
  store <- pcodec_open_store_cached(store)
  if (!store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS) {
    stop("this CompreSSoR build reads native 0.4 stores only; the historical Python-backed store is archived",
         call. = FALSE)
  }
  pcodec_native_projection(
    pcodec_native_read_store(store, region = region, variants = variants,
                             columns = columns, threads = threads),
    columns = columns
  )
}

# Internal trace counters: number of identity resolutions performed by batched
# reads. Tests use this to assert identity work happens once per panel group.
.pcodec_batch_trace <- new.env(parent = emptyenv())
.pcodec_batch_trace$identity_resolutions <- 0L
.pcodec_batch_trace$groups <- 0L

PCODEC_IDENTITY_COLUMNS <- c("global_position", "substitution", "chromosome",
                             "base_pair_location", "reference_allele",
                             "alternate_allele", "effect_allele", "other_allele")

# Identity signature of a store: SHA-256 of the position and substitution
# streams (from the manifest integrity record) plus row count and build. Two
# stores with equal signatures have identical row-id -> variant mappings.
# Stores without an integrity record get a unique signature (no sharing).
pcodec_identity_signature <- function(store) {
  m <- store$manifest
  files <- m$integrity$files
  pos <- m$files$position
  sub <- m$files$substitution
  hp <- if (!is.null(pos)) files[[pos]]$sha256 else NULL
  hs <- if (!is.null(sub)) files[[sub]]$sha256 else NULL
  if (is.null(hp) || is.null(hs)) return(paste0("unique:", normalizePath(store$path)))
  # Positions are delta-coded within key blocks, so the stream bytes alone do
  # not fix the identity: the block anchors (row ranges, first and last
  # positions) and the position encoding, from the index, must match too.
  # They are taken from the cached block matrices (no parsed index needed).
  mats <- pcodec_native_select_matrices(store)
  anchors <- unname(mats$position[, 4:7, drop = FALSE])
  paste(hp, hs, as.character(m$n_rows %||% m$rows),
        as.character(m$genome_build %||% "GRCh38"),
        if (isTRUE(mats$delta)) "delta_u32_within_block" else "",
        digest::digest(anchors, algo = "sha1"), sep = "|")
}

# Batched reads. Two strategies, both returning what read_sumstats() returns
# for each store:
# * one pass per store (the default): each forked worker opens its store and
#   reads keys -> rows -> values in one native call; a shared key list is
#   normalised once and parsed once per genome build per worker.
# * panel sharing (opt-in, options(CompreSSoR.batch_share_panels = TRUE)):
#   stores with the same variant panel (equal identity signature) resolve
#   keys, row IDs and regions to rows once, then decode values only, stores
#   in parallel with threads %/% length(stores) decoder threads each; stores
#   alone in their group are read in one pass each. Members are verified
#   against their own manifests first (pcodec_verify_identity_files()).
# Sharing was the default while selective reads decoded keys in R. With the
# native selective reader a per-store key resolution costs milliseconds, and
# sharing needs a parent-side grouping pass plus identity verification: on
# BluePebble (8 threads) one pass per store was as fast on 20 shared-panel
# simulated stores (1k keys) and faster on 300 UKB-PPP stores (6k keys).
pcodec_batch_share_panels <- function(k, threads) {
  forced <- getOption("CompreSSoR.batch_share_panels", NULL)
  if (!is.null(forced)) return(k > 1L && isTRUE(forced))
  FALSE
}

pcodec_read_stores <- function(stores, variants = NULL, columns, threads = 1L,
                               region = NULL, annotate = FALSE) {
  if (!length(stores)) stop("stores must be non-empty", call. = FALSE)
  k <- length(stores)
  if (!is.list(variants) || is.data.frame(variants)) variants <- rep(list(variants), k)
  if (is.null(variants) || length(variants) != k) {
    stop("stores and variants must have the same non-zero length", call. = FALSE)
  }
  if (!is.list(region)) region <- rep(list(region), k)
  if (length(region) != k) stop("region must be one value or one per store", call. = FALSE)
  if (!length(columns) || anyNA(columns) || any(!nzchar(columns))) {
    stop("columns must contain at least one non-empty column name", call. = FALSE)
  }
  columns <- unique(as.character(columns))
  threads <- pcodec_validate_threads(threads)
  share <- pcodec_batch_share_panels(k, threads)
  store_names <- names(stores)
  normalise <- function(keys) {
    if (is.null(keys)) return(NULL)
    if (is.character(keys)) {
      if (anyNA(keys) || any(!nzchar(trimws(keys)))) {
        stop("each variants element must contain canonical variant keys", call. = FALSE)
      }
      return(unique(trimws(keys)))
    }
    if (is.numeric(keys)) return(unique(pcodec_row_ids(keys)))
    stop("each variants element must contain canonical variant keys or row IDs",
         call. = FALSE)
  }
  # Normalise each distinct request once. A request shared by every store (the
  # usual case) is one object, so identical() is a pointer comparison.
  distinct <- list()
  slot <- integer(k)
  for (i in seq_len(k)) {
    j <- 0L
    for (u in seq_along(distinct)) {
      if (identical(distinct[[u]]$raw, variants[[i]])) { j <- u; break }
    }
    if (!j) {
      distinct[[length(distinct) + 1L]] <- list(raw = variants[[i]],
                                                norm = normalise(variants[[i]]))
      j <- length(distinct)
    }
    slot[i] <- j
  }
  variants <- lapply(slot, function(j) distinct[[j]]$norm)
  id_cols <- intersect(columns, PCODEC_IDENTITY_COLUMNS)
  value_cols <- setdiff(columns, id_cols)
  rows_only <- vapply(seq_len(k), function(i) {
    is.numeric(variants[[i]]) && is.null(region[[i]]) && !length(id_cols)
  }, logical(1))
  full <- vapply(seq_len(k), function(i) {
    is.null(variants[[i]]) && is.null(region[[i]])
  }, logical(1))

  if (!share) {
    # One pass per store, in parallel: each worker opens its store and reads
    # keys -> rows -> values in one native call (one fork round, no parent
    # work per store). A shared key list is parsed once per genome build per
    # worker process.
    resolutions <- sum(!rows_only & !full)
    .pcodec_batch_trace$identity_resolutions <-
      .pcodec_batch_trace$identity_resolutions + resolutions
    .pcodec_batch_trace$groups <- .pcodec_batch_trace$groups + resolutions
    labels <- vapply(seq_len(k), function(i) {
      candidates_store_label(stores[[i]], store_names[i])
    }, character(1))
    parse_cache <- new.env(parent = emptyenv())
    inner <- max(1L, threads %/% k)
    decoded <- pcodec_parallel_lapply(seq_len(k), function(i) {
      store <- pcodec_open_store_cached(stores[[i]])
      if (!store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS) {
        stop("this CompreSSoR build reads native 0.4 stores only; the historical ",
             "Python-backed store is archived", call. = FALSE)
      }
      targets <- NULL
      if (is.character(variants[[i]])) {
        build <- compressor_normalize_build(store$manifest$genome_build %||% "GRCh38")
        ck <- paste(slot[i], build)
        targets <- parse_cache[[ck]]
        if (is.null(targets)) {
          targets <- pcodec_native_target_keys(variants[[i]], build = build, trim = TRUE)
          parse_cache[[ck]] <- targets
        }
      }
      out <- pcodec_native_projection(
        pcodec_native_read_store(store, region = region[[i]], variants = variants[[i]],
                                 columns = columns, threads = inner,
                                 key_targets = targets),
        columns = columns)
      if (isTRUE(annotate)) attr(out, "compressor_store_meta") <- pcodec_store_meta(store)
      out
    }, threads = threads, labels = labels, what = "batched read")
    pcodec_batch_check_workers(decoded, k, labels)
    names(decoded) <- store_names
    return(decoded)
  }

  opened <- pcodec_batch_open(stores, threads, signatures = TRUE)
  stores <- opened$stores
  store_labels <- vapply(stores, function(store) store$path, character(1))

  # One pass per store for the stores in `idx`, in parallel; a shared key list
  # is parsed once per genome build. Worker errors are re-raised.
  parsed <- list()
  one_pass <- function(idx) {
    if (!length(idx)) return(list())
    targets <- vector("list", length(idx))
    for (t in seq_along(idx)) {
      i <- idx[t]
      if (!is.character(variants[[i]])) next
      build <- compressor_normalize_build(stores[[i]]$manifest$genome_build %||% "GRCh38")
      ck <- paste(slot[i], build)
      if (is.null(parsed[[ck]])) {
        parsed[[ck]] <<- pcodec_native_target_keys(variants[[i]], build = build,
                                                   trim = TRUE)
      }
      targets[t] <- list(parsed[[ck]])
    }
    inner <- max(1L, threads %/% length(idx))
    out <- pcodec_parallel_lapply(seq_along(idx), function(t) {
      i <- idx[t]
      pcodec_native_projection(
        pcodec_native_read_store(stores[[i]], region = region[[i]],
                                 variants = variants[[i]], columns = columns,
                                 threads = inner, key_targets = targets[[t]]),
        columns = columns)
    }, threads = threads, labels = store_labels[idx], what = "batched read")
    pcodec_batch_check_workers(out, length(idx), store_labels[idx])
  }

  signature <- opened$signature
  # Resolve each distinct (identity group, request) once. Pure value reads
  # of an explicit row-id request need no identity work at all. Each distinct
  # request is hashed once.
  request_hash <- list()
  request_key <- vapply(seq_len(k), function(i) {
    ck <- paste(slot[i], paste(region[[i]], collapse = "\r"), sep = "\n")
    if (is.null(request_hash[[ck]])) {
      request_hash[[ck]] <<- digest::digest(list(variants[[i]], region[[i]]),
                                            algo = "sha1")
    }
    request_hash[[ck]]
  }, character(1))
  resolve_key <- paste(signature, request_key, sep = "#")
  group_size <- table(resolve_key)
  # A store alone in its group needs no shared resolution: those stores are
  # read in one pass each, in parallel.
  single <- !rows_only & !full & as.integer(group_size[resolve_key]) == 1L
  # Stores whose rows (and identity columns) come from another member of
  # their group: check their identity streams against their own manifests
  # first (see pcodec_verify_identity_files()).
  pcodec_verify_identity_files(stores[which(!single & !rows_only & !full)],
                               threads = threads)
  resolved <- list()
  for (i in seq_len(k)) {
    key <- resolve_key[[i]]
    if (single[i] || !is.null(resolved[[key]])) next
    if (rows_only[i]) {
      resolved[[key]] <- list(rows = sort(variants[[i]]), identity = NULL, bytes = 0,
                              direct = TRUE)
      next
    }
    if (full[i]) {
      resolved[[key]] <- list(full = TRUE)
      next
    }
    .pcodec_batch_trace$identity_resolutions <- .pcodec_batch_trace$identity_resolutions + 1L
    ident <- pcodec_native_read_store(
      stores[[i]], region = region[[i]], variants = variants[[i]],
      columns = if (length(id_cols)) id_cols else "base_pair_location",
      threads = threads
    )
    resolved[[key]] <- list(rows = as.integer(ident$row), identity = ident,
                            bytes = attr(ident, "source_bytes_read", exact = TRUE) %||% 0)
  }
  .pcodec_batch_trace$identity_resolutions <-
    .pcodec_batch_trace$identity_resolutions + sum(single)
  .pcodec_batch_trace$groups <- .pcodec_batch_trace$groups + length(unique(signature))

  decoded <- vector("list", k)
  decoded[which(single)] <- one_pass(which(single))
  need <- which(!single)
  inner <- max(1L, threads %/% max(1L, length(need)))
  decoded[need] <- pcodec_parallel_lapply(need, function(i) {
    res <- resolved[[resolve_key[[i]]]]
    if (isTRUE(res$full)) {
      return(pcodec_read_store(stores[[i]], columns = columns, threads = inner))
    }
    if (!length(res$rows)) {
      return(pcodec_native_projection(pcodec_native_empty_result(columns), columns))
    }
    out <- NULL
    bytes <- res$bytes
    if (length(value_cols)) {
      vals <- pcodec_native_read_store(stores[[i]], variants = res$rows,
                                       columns = value_cols, threads = inner)
      bytes <- bytes + (attr(vals, "source_bytes_read", exact = TRUE) %||% 0)
      out <- if (is.null(res$identity)) vals else {
        cbind(res$identity[setdiff(names(res$identity), "row")], vals[value_cols])
      }
      if (!"row" %in% names(out)) out <- cbind(row = vals$row, out)
    } else {
      out <- res$identity
    }
    row.names(out) <- NULL
    attr(out, "source_bytes_read") <- bytes
    pcodec_native_projection(out, columns)
  }, threads = threads, labels = store_labels[need], what = "batched read")
  # Shared-panel workers must not hand back anything but a store's data.
  pcodec_batch_check_workers(decoded, k, store_labels)
  if (isTRUE(annotate)) {
    for (i in seq_len(k)) {
      attr(decoded[[i]], "compressor_store_meta") <- pcodec_store_meta(stores[[i]])
    }
  }
  names(decoded) <- store_names
  decoded
}

# Build and row count of an opened store, for request_index.
pcodec_store_meta <- function(store) {
  m <- store$manifest
  list(build = compressor_normalize_build(m$genome_build %||% "GRCh38"),
       n_rows = as.numeric(m$n_rows %||% m$rows))
}

# Open the stores of a batch (and load their cached block matrices, plus the
# identity signature when panels are shared). With several stores and
# threads this runs in forked workers, which parses manifests and indexes in
# parallel; the parent then keeps each store's compact block matrices in its
# cache (never the parsed index), so later stages and forks reuse them.
pcodec_batch_open <- function(stores, threads, signatures = FALSE) {
  labels <- vapply(seq_along(stores), function(i) {
    candidates_store_label(stores[[i]], names(stores)[i])
  }, character(1))
  open_one <- function(x) {
    store <- pcodec_open_store_cached(x)
    if (!store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS) {
      stop("this CompreSSoR build reads native 0.4 stores only; the historical ",
           "Python-backed store is archived", call. = FALSE)
    }
    mats <- pcodec_native_select_matrices(store)
    list(store = store, mats = mats,
         signature = if (signatures) pcodec_identity_signature(store) else NA_character_)
  }
  parallel <- length(stores) >= PCODEC_BATCH_PARALLEL_OPEN && threads > 1L
  out <- pcodec_parallel_lapply(stores, open_one,
                                threads = if (parallel) threads else 1L,
                                labels = labels, what = "store open")
  if (parallel) {
    for (item in out) pcodec_native_register_matrices(item$store, item$mats)
  }
  list(stores = lapply(out, `[[`, "store"),
       signature = vapply(out, `[[`, character(1), "signature"))
}

# Batches of at least this many stores are opened in parallel workers.
PCODEC_BATCH_PARALLEL_OPEN <- 8L

# Every batched result must be a data frame (pcodec_parallel_lapply() has
# already re-raised worker errors and rejected dead workers).
pcodec_batch_check_workers <- function(out, expected, labels = NULL) {
  if (length(out) != expected) {
    stop("batched read returned ", length(out), " results for ", expected,
         " stores", call. = FALSE)
  }
  missing <- !vapply(out, is.data.frame, logical(1))
  if (any(missing)) {
    at <- which(missing)[1L]
    stop("batched read worker for store ",
         if (is.null(labels)) at else paste0("'", labels[[at]], "'"),
         " returned no data (worker failed)", call. = FALSE)
  }
  out
}

pcodec_validate_store <- function(store, full = FALSE) {
  store <- pcodec_open_store_cached(store)
  if (!store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS) {
    stop("this CompreSSoR build validates native 0.4 stores only; the historical Python-backed store is archived",
         call. = FALSE)
  }
  pcodec_native_validate_store(store, full = full)
}
