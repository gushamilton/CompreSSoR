#' Return the CompreSSoR reference cache directory
#'
#' @param cache_dir Optional cache directory. If omitted, the
#'   COMPRESSOR_REFERENCE_CACHE environment variable is used, followed by
#'   the standard per-user R data directory.
#' @return A directory path.
#' @export
reference_cache_dir <- function(cache_dir = NULL) {
  cache_dir <- cache_dir %||% Sys.getenv("COMPRESSOR_REFERENCE_CACHE", unset = "")
  if (!nzchar(cache_dir)) cache_dir <- tools::R_user_dir("CompreSSoR", which = "cache")
  normalizePath(path.expand(cache_dir), mustWork = FALSE)
}

#' Describe the default GRCh38 reference
#'
#' The default is an external, immutable EBI/Ensembl release-95 GRCh38
#' reference, like a CRAM reference. Set `COMPRESSOR_CANONICAL_REFERENCE` to
#' the partitioned dictionary directory produced by `build_ebi_reference()`.
#'
#' @param cache_dir Optional cache directory to use when resolving the
#'   descriptor.
#' @return A reference descriptor.
#' @export
grch38_reference <- function(cache_dir = NULL) {
  local <- Sys.getenv("COMPRESSOR_CANONICAL_REFERENCE", unset = "")
  variants <- if (nzchar(local)) {
    path.expand(local)
  } else {
    NULL
  }
  list(
    id = "ebi_ensembl95_grch38_all_v1",
    build = "GRCh38",
    source = "EBI/Ensembl GRCh38 all-variant reference; build with build_ebi_reference()",
    cache_dir = cache_dir,
    variants = variants
  )
}

reference_file_digest <- function(path, algo = "sha256") {
  digest::digest(path, algo = algo, file = TRUE)
}

.reference_integrity_state <- new.env(parent = emptyenv())

verify_partitioned_reference <- function(path, manifest) {
  relative <- manifest$variants_file %||% "variants.parquet"
  master <- file.path(path, relative)
  expected <- manifest$sha256 %||% ""
  if (!file.exists(master)) {
    stop("partitioned reference is missing its master variant file: ", master,
         call. = FALSE)
  }
  if (!is.character(expected) || length(expected) != 1L ||
      !grepl("^[0-9a-fA-F]{64}$", expected)) {
    stop("partitioned reference manifest has no valid SHA-256", call. = FALSE)
  }
  info <- file.info(master)
  cache_key <- paste(normalizePath(master, mustWork = TRUE), info$size,
                     as.numeric(info$mtime), tolower(expected), sep = "\r")
  if (isTRUE(.reference_integrity_state[[cache_key]])) return(invisible(TRUE))
  observed <- reference_file_digest(master, "sha256")
  if (!identical(tolower(observed), tolower(expected))) {
    stop("partitioned reference master file failed its declared SHA-256: ", master,
         call. = FALSE)
  }
  .reference_integrity_state[[cache_key]] <- TRUE
  invisible(TRUE)
}

reference_asset_filename <- function(asset) {
  filename <- asset$filename %||% basename(sub("[?].*$", "", asset$url %||% ""))
  filename <- basename(filename)
  if (!nzchar(filename) || identical(filename, ".")) {
    stop("reference download asset needs a filename", call. = FALSE)
  }
  filename
}

reference_asset_path <- function(asset, cache_dir) {
  if (!is.null(asset$path)) return(normalizePath(path.expand(asset$path), mustWork = TRUE))
  file.path(cache_dir, reference_asset_filename(asset))
}

verify_reference_asset <- function(path, asset) {
  if (!file.exists(path)) return(FALSE)
  if (!is.null(asset$md5) && nzchar(asset$md5)) {
    if (!identical(tolower(reference_file_digest(path, "md5")), tolower(asset$md5))) return(FALSE)
  }
  if (!is.null(asset$sha256) && nzchar(asset$sha256)) {
    if (!identical(tolower(reference_file_digest(path, "sha256")), tolower(asset$sha256))) return(FALSE)
  }
  TRUE
}

download_reference_asset <- function(asset, cache_dir, overwrite = FALSE) {
  if (!is.null(asset$path)) {
    path <- reference_asset_path(asset, cache_dir)
    if (!verify_reference_asset(path, asset)) {
      stop("local reference asset failed its checksum: ", path, call. = FALSE)
    }
    return(list(path = path, downloaded = FALSE))
  }
  if (is.null(asset$url) || !nzchar(asset$url)) {
    stop("reference asset needs either path or url", call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  target <- reference_asset_path(asset, cache_dir)
  if (file.exists(target) && !isTRUE(overwrite)) {
    if (verify_reference_asset(target, asset)) {
      return(list(path = target, downloaded = FALSE))
    }
    stop("cached reference asset failed its checksum; use overwrite = TRUE: ", target,
         call. = FALSE)
  }
  temporary <- tempfile("compressor-reference-", tmpdir = cache_dir)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout %||% 60))
  on.exit(options(timeout = old_timeout), add = TRUE)
  status <- tryCatch(
    utils::download.file(asset$url, temporary, mode = "wb", quiet = TRUE,
                         method = "libcurl"),
    error = function(e) stop("could not download reference asset: ", conditionMessage(e), call. = FALSE)
  )
  if (!identical(status, 0L) || !file.exists(temporary)) {
    stop("could not download reference asset: ", asset$url, call. = FALSE)
  }
  if (!verify_reference_asset(temporary, asset)) {
    stop("downloaded reference asset failed its checksum: ", asset$url, call. = FALSE)
  }
  if (file.exists(target)) unlink(target, force = TRUE)
  if (!file.rename(temporary, target)) stop("could not install reference asset: ", target, call. = FALSE)
  list(path = target, downloaded = TRUE)
}

reference_descriptor <- function(reference) {
  if (is.character(reference) && length(reference) == 1L &&
      toupper(reference) %in% c("GRCH38", "HG38")) {
    return(grch38_reference())
  }
  if (is.data.frame(reference)) {
    return(list(id = "in_memory", build = "GRCh38", variants = reference,
                metadata = list(id = "in_memory", build = "GRCh38")))
  }
  if (is.character(reference) && length(reference) == 1L && file.exists(reference)) {
    return(list(id = "local", build = "GRCh38", variants = reference,
                metadata = list(id = "local", build = "GRCh38", local_path = normalizePath(reference))))
  }
  if (is.list(reference)) {
    if (is.null(reference$build)) reference$build <- "GRCh38"
    if (is.null(reference$variants) &&
        is.null(reference$assets) &&
        !identical(reference$id %||% "", "ebi_ensembl95_grch38_all_v1") &&
        !identical(reference$id %||% "", "canonical_grch38_snp_v1")) {
      stop("reference descriptor needs a variants asset", call. = FALSE)
    }
    return(reference)
  }
  stop("reference must be 'GRCh38', a GRCh38 data.frame/path, or a reference descriptor", call. = FALSE)
}

#' Resolve or download a configured reference into the CompreSSoR cache
#'
#' @param reference "GRCh38", a reference descriptor, or a local reference
#'   path. The built-in `"GRCh38"` descriptor requires
#'   `COMPRESSOR_CANONICAL_REFERENCE` to point at a dictionary built by
#'   [build_ebi_reference()]; it has no implicit remote download URL.
#' @param cache_dir Optional cache directory.
#' @param overwrite Re-download assets that are already cached.
#' @return A resolved reference descriptor containing local asset paths and
#'   checksums.
#' @export
download_reference <- function(reference = "GRCh38", cache_dir = NULL, overwrite = FALSE) {
  descriptor <- reference_descriptor(reference)
  build <- toupper(gsub("[.]", "", descriptor$build %||% "GRCh38"))
  if (!identical(build, "GRCH38")) {
    stop("CompreSSoR currently requires a GRCh38 reference", call. = FALSE)
  }
  cache_dir <- cache_dir %||% descriptor$cache_dir %||% reference_cache_dir()
  descriptor$cache_dir <- cache_dir
  assets <- descriptor$variants
  if (is.null(assets) && !is.null(descriptor$assets)) {
    stop("EBI reference has remote chromosome assets; call build_ebi_reference() to materialise the dictionary",
         call. = FALSE)
  }
  if (is.data.frame(assets)) {
    descriptor$metadata <- descriptor$metadata %||% list(id = descriptor$id %||% "in_memory", build = "GRCh38")
    descriptor$metadata$local_path <- NULL
    descriptor$metadata$rows <- nrow(assets)
    descriptor$variants <- assets
    return(descriptor)
  }
  if (is.character(assets) && length(assets) == 1L && dir.exists(assets)) {
    descriptor$variants <- normalizePath(assets, mustWork = TRUE)
    manifest_path <- file.path(descriptor$variants, "manifest.json")
    if (!file.exists(manifest_path)) {
      stop("partitioned reference is missing manifest.json: ", descriptor$variants, call. = FALSE)
    }
    manifest <- read_manifest(manifest_path)
    verify_partitioned_reference(descriptor$variants, manifest)
    descriptor$id <- manifest$id %||% descriptor$id
    descriptor$build <- manifest$build %||% descriptor$build
    descriptor$metadata <- utils::modifyList(descriptor$metadata %||% list(), manifest)
    descriptor$metadata$local_path <- descriptor$variants
    descriptor$metadata$filename <- basename(descriptor$variants)
    return(descriptor)
  } else if (is.character(assets) && length(assets) == 1L && file.exists(assets)) {
    asset <- list(path = assets)
  } else if (is.character(assets) && length(assets) == 1L) {
    asset <- list(url = assets)
  } else if (is.list(assets) && (!is.null(assets$url) || !is.null(assets$path))) {
    asset <- assets
  } else {
    stop("no canonical GRCh38 reference is configured; set COMPRESSOR_CANONICAL_REFERENCE or build the EBI dictionary with build_ebi_reference()",
         call. = FALSE)
  }
  resolved <- download_reference_asset(asset, cache_dir, overwrite = overwrite)
  descriptor$variants <- resolved$path
  metadata <- descriptor$metadata %||% list()
  manifest_path <- if (grepl("[.]parquet$", resolved$path, ignore.case = TRUE)) {
    sub("[.]parquet$", ".manifest.json", resolved$path, ignore.case = TRUE)
  } else NULL
  if (!is.null(manifest_path) && file.exists(manifest_path)) {
    manifest <- read_manifest(manifest_path)
    metadata <- utils::modifyList(metadata, manifest)
    descriptor$id <- manifest$id %||% descriptor$id
    descriptor$build <- manifest$build %||% descriptor$build
  }
  metadata$id <- metadata$id %||% descriptor$id %||% "custom"
  metadata$build <- descriptor$build %||% "GRCh38"
  metadata$source <- descriptor$source %||% NULL
  metadata$source_url <- asset$url %||% descriptor$source_url %||% NULL
  metadata$local_path <- normalizePath(resolved$path, mustWork = FALSE)
  metadata$md5 <- reference_file_digest(resolved$path, "md5")
  metadata$sha256 <- reference_file_digest(resolved$path, "sha256")
  metadata$downloaded <- isTRUE(resolved$downloaded)
  metadata$filename <- basename(resolved$path)
  source_format <- tryCatch(reference_source_format(resolved$path), error = function(e) NULL)
  source_metadata <- tryCatch(read_reference_metadata(resolved$path, source_format),
                              error = function(e) NULL)
  if (!is.null(source_format) && !identical(source_format, "delimited")) {
    metadata$source_format <- source_format
  }
  if (!is.null(source_metadata)) metadata$source_metadata <- source_metadata
  descriptor$metadata <- metadata
  descriptor
}

#' Resolve a configured reference, downloading an explicit remote asset when necessary
#'
#' @param reference "GRCh38", a reference descriptor, or a local reference
#'   path. `"GRCh38"` requires a preconfigured local canonical dictionary;
#'   only descriptors with an explicit URL are downloaded.
#' @param cache_dir Optional cache directory.
#' @return A resolved reference descriptor.
#' @export
resolve_reference <- function(reference = "GRCh38", cache_dir = NULL) {
  descriptor <- reference_descriptor(reference)
  if (is.data.frame(descriptor$variants)) return(descriptor)
  download_reference(descriptor, cache_dir = cache_dir)
}
