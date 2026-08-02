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
#' The default is the public GWASLab 1KG dbSNP151 GRCh38 autosomal variant
#' table. It is deliberately kept outside the package itself, like a CRAM
#' reference: the first operation that needs it downloads it into the local
#' CompreSSoR cache and records the URL and checksum in the store manifest.
#'
#' @param cache_dir Optional cache directory to use when resolving the
#'   descriptor.
#' @return A reference descriptor.
#' @export
grch38_reference <- function(cache_dir = NULL) {
  list(
    id = "1kg_dbsnp151_hg38_auto",
    build = "GRCh38",
    source = "GWASLab reference catalogue",
    source_url = "https://github.com/Cloufield/gwaslab/blob/main/src/gwaslab/data/reference.json",
    cache_dir = cache_dir,
    variants = list(
      url = "https://www.dropbox.com/s/ouf60n7gdz6cm0g/1kg_dbsnp151_hg38_auto.txt.gz?dl=1",
      filename = "1kg_dbsnp151_hg38_auto.txt.gz",
      md5 = "4c7ef2d2415c18c286219e970fdda972",
      description = "1KG dbSNP151 GRCh38 autosomal variant table"
    )
  )
}

reference_file_digest <- function(path, algo = "sha256") {
  digest::digest(path, algo = algo, file = TRUE)
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
  status <- tryCatch(
    utils::download.file(asset$url, temporary, mode = "wb", quiet = TRUE),
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
    if (is.null(reference$variants)) stop("reference descriptor needs a variants asset", call. = FALSE)
    return(reference)
  }
  stop("reference must be 'GRCh38', a GRCh38 data.frame/path, or a reference descriptor", call. = FALSE)
}

#' Download a reference into the CompreSSoR cache
#'
#' @param reference "GRCh38", a reference descriptor, or a local reference
#'   path.
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
  if (is.data.frame(assets)) {
    descriptor$metadata <- descriptor$metadata %||% list(id = descriptor$id %||% "in_memory", build = "GRCh38")
    descriptor$metadata$local_path <- NULL
    descriptor$metadata$rows <- nrow(assets)
    descriptor$variants <- assets
    return(descriptor)
  }
  if (is.character(assets) && length(assets) == 1L && file.exists(assets)) {
    asset <- list(path = assets)
  } else if (is.character(assets) && length(assets) == 1L) {
    asset <- list(url = assets)
  } else if (is.list(assets) && (!is.null(assets$url) || !is.null(assets$path))) {
    asset <- assets
  } else {
    stop("reference variants must be a data.frame, local path, URL, or asset descriptor", call. = FALSE)
  }
  resolved <- download_reference_asset(asset, cache_dir, overwrite = overwrite)
  descriptor$variants <- resolved$path
  metadata <- descriptor$metadata %||% list()
  metadata$id <- descriptor$id %||% metadata$id %||% "custom"
  metadata$build <- descriptor$build %||% "GRCh38"
  metadata$source <- descriptor$source %||% NULL
  metadata$source_url <- asset$url %||% descriptor$source_url %||% NULL
  metadata$local_path <- normalizePath(resolved$path, mustWork = FALSE)
  metadata$md5 <- reference_file_digest(resolved$path, "md5")
  metadata$sha256 <- reference_file_digest(resolved$path, "sha256")
  metadata$downloaded <- isTRUE(resolved$downloaded)
  metadata$filename <- basename(resolved$path)
  descriptor$metadata <- metadata
  descriptor
}

#' Resolve a reference, downloading it when necessary
#'
#' @param reference "GRCh38", a reference descriptor, or a local reference
#'   path.
#' @param cache_dir Optional cache directory.
#' @return A resolved reference descriptor.
#' @export
resolve_reference <- function(reference = "GRCh38", cache_dir = NULL) {
  descriptor <- reference_descriptor(reference)
  if (is.data.frame(descriptor$variants)) return(descriptor)
  download_reference(descriptor, cache_dir = cache_dir)
}
