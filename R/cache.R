#' Build an optional q8 framed cache
#'
#' The cache uses independently compressed gzip frames so it is portable with
#' base R. It is a serving layer for repeated regional/MR access, not the
#' archival source of truth.
#'
#' @param store A CompreSSoR store object or path.
#' @param output Cache directory. Defaults to `cache.q8` inside the store.
#' @param block_rows Rows per independent frame.
#' @param overwrite Whether an existing cache may be replaced.
#' @return The cache directory path, invisibly.
#' @export
build_cache <- function(store, output = NULL, block_rows = 65536L, overwrite = FALSE) {
  store <- if (inherits(store, "compressor_store")) store else open_compressor(store)
  if (identical(store$manifest$backend, "pcodec")) {
    stop("Pcodec stores are already independently paged and do not use q8 caches",
         call. = FALSE)
  }
  output <- output %||% file.path(store$path, "cache.q8")
  if (dir.exists(output)) {
    if (!overwrite) stop("cache already exists; use overwrite = TRUE", call. = FALSE)
    unlink(output, recursive = TRUE, force = TRUE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  data <- read_sumstats(store, use_cache = FALSE)
  n <- nrow(data)
  cache_meta <- q_encode(data$beta[seq_len(min(n, 1L))], data$standard_error[seq_len(min(n, 1L))],
                         data$effect_allele_frequency[seq_len(min(n, 1L))],
                         z = data$z[seq_len(min(n, 1L))],
                         z_bits = 8L, se_bits = 8L, eaf_bits = 8L,
                         block_rows = block_rows)$metadata
  frame_path <- file.path(output, "frames.bin")
  con <- file(frame_path, open = "wb")
  on.exit(close(con), add = TRUE)
  index <- vector("list", ceiling(n / block_rows))
  for (i in seq_along(index)) {
    start <- (i - 1L) * block_rows + 1L
    stop <- min(n, i * block_rows)
    encoded <- q_encode(data$beta[start:stop], data$standard_error[start:stop],
                        data$effect_allele_frequency[start:stop],
                        z = data$z[start:stop], z_bits = 8L, se_bits = 8L,
                        eaf_bits = 8L, block_rows = block_rows)
    payload <- list(
      row_start = start - 1L,
      main = encoded$main,
      exceptions = encoded$exceptions,
      metadata = encoded$metadata
    )
    compressed <- memCompress(serialize(payload, NULL, version = 3), type = "gzip")
    length_offset <- seek(con, where = 0, origin = "current")
    writeBin(as.integer(length(compressed)), con, size = 4L, endian = "little")
    frame_offset <- seek(con, where = 0, origin = "current")
    writeBin(compressed, con)
    index[[i]] <- list(row_start = start - 1L, row_end = stop - 1L,
                       offset = frame_offset, length = length(compressed),
                       length_offset = length_offset)
  }
  close(con)
  on.exit(NULL, add = FALSE)
  saveRDS(index, file.path(output, "index.rds"))
  manifest <- list(
    format = "CompreSSoR-q8-cache", format_version = "0.1.0",
    profile = "discovery", backend = "framed-gzip", block_rows = block_rows,
    n_rows = n, files = list(frames = "frames.bin", index = "index.rds"),
    codec = cache_meta, source = basename(store$path), created_utc = now_utc()
  )
  write_manifest(manifest, file.path(output, "manifest.json"))
  invisible(output)
}

read_q8_cache <- function(store, region = NULL) {
  cache_path <- file.path(store$path, "cache.q8")
  if (!dir.exists(cache_path)) stop("q8 cache does not exist; call build_cache() first", call. = FALSE)
  manifest <- read_manifest(file.path(cache_path, "manifest.json"))
  index <- readRDS(file.path(cache_path, manifest$files$index))
  bounds <- read_region_bounds(region)
  variants <- read_variant_table(store, bounds = bounds)
  if (!nrow(variants)) return(data.frame())
  wanted <- variants$row
  selected <- vapply(index, function(x) any(wanted >= x$row_start & wanted <= x$row_end), logical(1))
  con <- file(file.path(cache_path, manifest$files$frames), open = "rb")
  on.exit(close(con), add = TRUE)
  decoded <- vector("list", sum(selected))
  k <- 0L
  for (i in which(selected)) {
    frame <- index[[i]]
    seek(con, where = frame$offset, origin = "start")
    compressed <- readBin(con, what = raw(), n = frame$length)
    payload <- unserialize(memDecompress(compressed, type = "gzip"))
    values <- q_decode(payload$main, payload$exceptions, payload$metadata,
                       include_p = TRUE)
    decoded[[k <- k + 1L]] <- data.frame(
      row = payload$row_start + seq_len(nrow(payload$main)) - 1L,
      beta = values$beta,
      z = values$z,
      standard_error = values$standard_error,
      effect_allele_frequency = values$effect_allele_frequency,
      p_value = values$p_value
    )
  }
  values <- do.call(rbind, decoded)
  out <- merge(variants, values, by = "row", sort = FALSE)
  out <- out[order(out$row), , drop = FALSE]
  out$row <- NULL
  out
}
