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

pcodec_canonical_manifest_sha256 <- function(manifest) {
  strip_observational <- function(value) {
    if (!is.list(value)) return(value)
    names_value <- names(value)
    if (!is.null(names_value)) {
      value <- value[!names_value %in% c("timings", "elapsed_seconds",
                                         "read_elapsed_seconds")]
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

verify_pcodec_manifest <- function(path) {
  checksum_path <- pcodec_manifest_checksum_path(path)
  if (!file.exists(checksum_path)) {
    stop("Pcodec store is missing manifest.sha256", call. = FALSE)
  }
  expected <- trimws(readLines(checksum_path, warn = FALSE, n = 1L))
  if (length(expected) != 1L || !grepl("^[0-9a-fA-F]{64}$", expected)) {
    stop("Pcodec manifest checksum record is malformed", call. = FALSE)
  }
  observed <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!identical(tolower(expected), tolower(observed))) {
    stop("Pcodec manifest checksum mismatch", call. = FALSE)
  }
  invisible(TRUE)
}

pcodec_open_store_cached <- function(store) {
  if (inherits(store, "compressor_store")) return(store)
  open_compressor(store)
}

pcodec_write_store <- function(data, output, metadata = list()) {
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
  pcodec_native_write_store(data, output, metadata = metadata)
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
  # Whole-file reads benefit from independent stream decoding in parallel;
  # regional and canonical-key reads are block-selective and avoid forking.
  if (!is.null(region) || !is.null(variants)) 1L else 4L
}

pcodec_parallel_lapply <- function(X, FUN, threads = 1L) {
  threads <- pcodec_validate_threads(threads)
  # R CMD check sets this guard to prevent packages from spawning an
  # uncontrolled number of workers. Respect it while retaining the native
  # four-thread default for ordinary whole-file reads.
  check_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", ""))
  if (nzchar(check_limit) && check_limit != "false") threads <- min(threads, 2L)
  if (length(X) <= 1L || threads <= 1L) return(lapply(X, FUN))
  # Forked workers are safe here because each worker opens independent files
  # and calls the standalone Pcodec decoder on private R objects. Windows has
  # no fork backend; retain deterministic serial behaviour there.
  if (.Platform$OS.type == "windows") return(lapply(X, FUN))
  parallel::mclapply(X, FUN, mc.cores = min(threads, length(X)),
                     mc.preschedule = TRUE)
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

pcodec_read_stores <- function(stores, variants, columns, threads = 1L) {
  if (!length(stores) || length(stores) != length(variants)) {
    stop("stores and variants must have the same non-zero length", call. = FALSE)
  }
  if (!length(columns) || anyNA(columns) || any(!nzchar(columns))) {
    stop("columns must contain at least one non-empty column name", call. = FALSE)
  }
  threads <- pcodec_validate_threads(threads)
  stores <- lapply(stores, pcodec_open_store_cached)
  if (any(!vapply(stores, function(store) {
    store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS
  }, logical(1)))) {
    stop("this CompreSSoR build reads native 0.4 stores only; the historical Python-backed store is archived",
         call. = FALSE)
  }
  variants <- lapply(variants, function(keys) {
    if (!is.character(keys) || anyNA(keys) || any(!nzchar(trimws(keys)))) {
      stop("each variants element must contain canonical variant keys", call. = FALSE)
    }
    unique(trimws(keys))
  })
  decoded <- pcodec_parallel_lapply(seq_along(stores), function(i) {
    pcodec_read_store(stores[[i]], variants = variants[[i]], columns = columns,
                      threads = 1L)
  }, threads = threads)
  names(decoded) <- names(stores)
  decoded
}

pcodec_validate_store <- function(store, full = FALSE) {
  store <- pcodec_open_store_cached(store)
  if (!store$manifest$format_version %in% PCODEC_NATIVE_SUPPORTED_FORMATS) {
    stop("this CompreSSoR build validates native 0.4 stores only; the historical Python-backed store is archived",
         call. = FALSE)
  }
  pcodec_native_validate_store(store, full = full)
}
