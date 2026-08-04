## Native Pcodec store interface.
##
## The historical Python-backed implementation lives under archive/ and is not
## part of the installed package. New stores are native 0.4 stores only.

pcodec_manifest_checksum_path <- function(path) {
  file.path(dirname(path), "manifest.sha256")
}

seal_pcodec_manifest <- function(path) {
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

pcodec_validate_threads <- function(threads) {
  if (length(threads) != 1L || !is.numeric(threads) || is.na(threads) ||
      !is.finite(threads) || threads < 1 || threads != floor(threads)) {
    stop("threads must be one positive integer", call. = FALSE)
  }
  as.integer(threads)
}

pcodec_parallel_lapply <- function(X, FUN, threads = 1L) {
  threads <- pcodec_validate_threads(threads)
  if (length(X) <= 1L || threads <= 1L) return(lapply(X, FUN))
  # Forked workers are safe here because each worker opens independent files
  # and calls the standalone Pcodec decoder on private R objects. Windows has
  # no fork backend; retain deterministic serial behaviour there.
  if (.Platform$OS.type == "windows") return(lapply(X, FUN))
  parallel::mclapply(X, FUN, mc.cores = min(threads, length(X)),
                     mc.preschedule = TRUE)
}

pcodec_read_store <- function(store, region = NULL, variants = NULL,
                               columns = NULL, threads = 1L) {
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
