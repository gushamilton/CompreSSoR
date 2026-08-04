## Native Pcodec backend.
##
## The native path uses the small standalone Pcodec ABI and a block index. It
## deliberately has a new format version: the upstream C ABI does not expose
## the wrapped FileCompressor API used by the older Python-backed stores.

PCODEC_NATIVE_FORMAT <- "0.4.4-pcodec-native"
PCODEC_NATIVE_SUPPORTED_FORMATS <- c("0.4.0-pcodec-native", "0.4.1-pcodec-native",
                                     "0.4.2-pcodec-native", "0.4.3-pcodec-native",
                                     PCODEC_NATIVE_FORMAT)
PCODEC_NATIVE_BLOCK_ROWS <- 65536L
PCODEC_NATIVE_KEY_BLOCK_ROWS <- 131072L
PCODEC_NATIVE_SE_CENTER_ROWS <- 65536L
PCODEC_NATIVE_PAGE_ROWS <- 131072L
PCODEC_NATIVE_LEVEL <- 8L
PCODEC_NATIVE_SE_BITS <- 8L
PCODEC_NATIVE_SE_COUNT <- 254L
PCODEC_NATIVE_SE_RESIDUAL_RANGE <- c(-4, 4)

pcodec_native_available <- function() {
  is.loaded("compressor_pcodec_native_available", PACKAGE = "CompreSSoR") &&
    isTRUE(.Call("compressor_pcodec_native_available", PACKAGE = "CompreSSoR"))
}

pcodec_native_enabled <- function() {
  pcodec_native_available()
}

pcodec_native_compress <- function(values, dtype) {
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  values <- if (dtype == "u32") as.numeric(values) else as.integer(values)
  switch(dtype,
    u8 = .Call("compressor_pcodec_compress_u8", values,
               PCODEC_NATIVE_LEVEL, PCODEC_NATIVE_PAGE_ROWS,
               PACKAGE = "CompreSSoR"),
    u16 = .Call("compressor_pcodec_compress_u16", values,
                PCODEC_NATIVE_LEVEL, PCODEC_NATIVE_PAGE_ROWS,
                PACKAGE = "CompreSSoR"),
    u32 = .Call("compressor_pcodec_compress_u32", values,
                PCODEC_NATIVE_LEVEL, PCODEC_NATIVE_PAGE_ROWS,
                PACKAGE = "CompreSSoR"),
    stop("unsupported native Pcodec dtype: ", dtype, call. = FALSE)
  )
}

pcodec_native_decompress <- function(blob, n, dtype) {
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  switch(dtype,
    u8 = .Call("compressor_pcodec_decompress_u8", blob, as.integer(n),
               PACKAGE = "CompreSSoR"),
    u16 = .Call("compressor_pcodec_decompress_u16", blob, as.integer(n),
                PACKAGE = "CompreSSoR"),
    u32 = .Call("compressor_pcodec_decompress_u32", blob, as.integer(n),
                PACKAGE = "CompreSSoR"),
    stop("unsupported native Pcodec dtype: ", dtype, call. = FALSE)
  )
}

pcodec_native_zstd_compress <- function(blob, level = 19L) {
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  .Call("compressor_zstd_compress", blob, as.integer(level),
        PACKAGE = "CompreSSoR")
}

pcodec_native_zstd_decompress <- function(blob, n) {
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  .Call("compressor_zstd_decompress", blob, as.integer(n),
        PACKAGE = "CompreSSoR")
}

pcodec_native_offsets <- function() {
  lengths <- compressor_grch38_chromosome_lengths
  c(0, cumsum(as.numeric(lengths)))[seq_along(lengths)]
}

pcodec_native_identity <- function(data) {
  chromosome <- toupper(sub("^CHR", "", as.character(data$chromosome),
                             ignore.case = TRUE))
  position <- as.numeric(data$base_pair_location)
  reference <- toupper(as.character(data$other_allele))
  alternate <- toupper(as.character(data$effect_allele))
  offsets <- pcodec_native_offsets()
  names(offsets) <- names(compressor_grch38_chromosome_lengths)
  base_code <- c(A = 0, C = 1, G = 2, T = 3)
  global_position <- unname(offsets[chromosome]) + position - 1
  substitution <- as.integer(base_code[reference]) * 4L +
    as.integer(base_code[alternate])
  if (any(!is.finite(global_position) | global_position < 0 | global_position > 4294967295)) {
    stop("native Pcodec identity position is outside uint32", call. = FALSE)
  }
  if (anyNA(substitution)) stop("native Pcodec identity contains an invalid allele", call. = FALSE)
  list(global_position = global_position, substitution = substitution)
}

pcodec_native_quantise <- function(data, block_rows = PCODEC_NATIVE_SE_CENTER_ROWS) {
  n <- nrow(data)
  z <- as.numeric(data$z)
  se <- as.numeric(data$standard_error)
  eaf <- as.numeric(data$effect_allele_frequency)
  valid_eaf <- is.finite(eaf) & eaf >= 0 & eaf <= 1
  safe_eaf <- pmin(1 - 1e-12, pmax(1e-12, ifelse(valid_eaf, eaf, 0.5)))
  eaf_codes <- as.integer(round(255 * (2 / pi) * asin(sqrt(safe_eaf))))
  eaf_decoded <- sin(pi * eaf_codes / (2 * 255))^2

  z_min <- -3.5
  z_max <- 3.5
  z_count <- 510L
  z_missing <- 510L
  z_exception <- 511L
  z_step <- (z_max - z_min) / z_count
  z_valid <- is.finite(z)
  z_central <- z_valid & z >= z_min & z < z_max
  z_codes <- rep.int(z_missing, n)
  z_codes[z_central] <- as.integer(floor((z[z_central] - z_min) / z_step))
  z_codes[!z_central & z_valid] <- z_exception

  safe_eaf_for_se <- pmin(1 - 1e-12, pmax(1e-12, eaf_decoded))
  valid_se <- is.finite(se) & se > 0
  residual <- rep(NA_real_, n)
  residual_ready <- valid_se & valid_eaf
  residual[residual_ready] <- log2(se[residual_ready]) +
    0.5 * log2(2 * safe_eaf_for_se[residual_ready] *
                 (1 - safe_eaf_for_se[residual_ready]))
  # The physical stream is already uint8. The previous native format used
  # only six semantic bits and therefore turned ordinary, non-template SEs
  # into exceptions. Use the full byte domain: 254 central bins plus one
  # missing and one exact-exception sentinel.
  se_count <- PCODEC_NATIVE_SE_COUNT
  se_missing <- se_count
  se_exception <- se_count + 1L
  se_min <- PCODEC_NATIVE_SE_RESIDUAL_RANGE[1]
  se_max <- PCODEC_NATIVE_SE_RESIDUAL_RANGE[2]
  se_step <- (se_max - se_min) / se_count
  se_codes <- rep.int(se_missing, n)
  centres <- numeric(if (n) ceiling(n / block_rows) else 0L)
  se_central <- is.finite(residual) & valid_se & valid_eaf
  if (n) {
    for (block in seq_along(centres)) {
      start <- (block - 1L) * block_rows + 1L
      stop <- min(n, block * block_rows)
      inside <- start:stop
      good <- inside[se_central[inside]]
      centres[block] <- if (length(good)) stats::median(residual[good], na.rm = TRUE) else 0
      delta <- residual[inside] - centres[block]
      in_range <- se_central[inside] & delta >= se_min & delta < se_max
      local <- rep.int(se_missing, length(inside))
      local[in_range] <- as.integer(floor((delta[in_range] - se_min) / se_step))
      local[se_central[inside] & !in_range] <- se_exception
      local[valid_se[inside] & !valid_eaf[inside]] <- se_exception
      se_codes[inside] <- local
    }
  }

  z_exception_mask <- !z_central & z_valid
  se_exception_mask <- se_codes == se_exception
  exception_mask <- z_exception_mask | se_exception_mask | !valid_eaf
  rows <- which(exception_mask) - 1L
  exceptions <- data.frame(
    row = as.integer(rows),
    z = as.numeric(z[rows + 1L]),
    log2se = as.numeric(log2(ifelse(is.finite(se[rows + 1L]) & se[rows + 1L] > 0,
                                    se[rows + 1L], 1))),
    eaf = as.numeric(eaf[rows + 1L]),
    flags = as.integer(z_exception_mask[rows + 1L]) +
      2L * as.integer(se_exception_mask[rows + 1L]) +
      4L * as.integer(!valid_eaf[rows + 1L])
  )
  list(
    z = as.integer(z_codes), eaf = as.integer(eaf_codes),
    se = as.integer(se_codes), centres = centres, exceptions = exceptions
  )
}

pcodec_native_append_stream <- function(values, path, dtype,
                                         block_rows = PCODEC_NATIVE_BLOCK_ROWS) {
  n <- length(values)
  blocks <- vector("list", if (n) ceiling(n / block_rows) else 0L)
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  offset <- 0
  if (n) {
    for (block in seq_along(blocks)) {
      start <- (block - 1L) * block_rows + 1L
      stop <- min(n, block * block_rows)
      blob <- pcodec_native_compress(values[start:stop], dtype)
      writeBin(blob, connection, useBytes = TRUE)
      blocks[[block]] <- list(
        row_start = start - 1L, row_stop = stop,
        offset = offset, length = length(blob), values = stop - start + 1L
      )
      offset <- offset + length(blob)
    }
  }
  list(file = basename(path), bytes = offset, blocks = blocks)
}

pcodec_native_position_gaps <- function(position, block_rows = PCODEC_NATIVE_BLOCK_ROWS) {
  n <- length(position)
  gaps <- numeric(n)
  if (n) {
    blocks <- ceiling(n / as.integer(block_rows))
    for (block in seq_len(blocks)) {
      start <- (block - 1L) * as.integer(block_rows) + 1L
      stop <- min(n, block * as.integer(block_rows))
      gaps[start] <- 0
      if (stop > start) gaps[(start + 1L):stop] <- diff(position[start:stop])
    }
  }
  gaps
}

pcodec_native_block_template <- function(position, block_rows) {
  n <- length(position)
  block_count <- if (n) ceiling(n / block_rows) else 0L
  lapply(seq_len(block_count), function(block) {
    start <- (block - 1L) * block_rows
    stop <- min(n, block * block_rows)
    list(row_start = start, row_stop = stop,
         first_position = position[start + 1L],
         last_position = position[stop])
  })
}

pcodec_native_validate_block_rows <- function(value, label) {
  value <- as.integer(value)
  if (length(value) != 1L || is.na(value) || value < 1024L ||
      value != 2^round(log2(value))) {
    stop(label, " must be a power of two >= 1024", call. = FALSE)
  }
  value
}

pcodec_native_index_blocks <- function(index, kind = c("value", "key")) {
  kind <- match.arg(kind)
  if (identical(kind, "key") && !is.null(index$key_blocks)) {
    return(index$key_blocks)
  }
  if (identical(kind, "value") && !is.null(index$value_blocks)) {
    return(index$value_blocks)
  }
  index$blocks
}

pcodec_native_exception_bytes <- function(exceptions) {
  if (!nrow(exceptions)) return(raw())
  connection <- rawConnection(raw(0), open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.integer(exceptions$row), connection, size = 4L, endian = "little")
  writeBin(as.numeric(exceptions$z), connection, size = 4L, endian = "little")
  writeBin(as.numeric(exceptions$log2se), connection, size = 4L, endian = "little")
  writeBin(as.numeric(exceptions$eaf), connection, size = 4L, endian = "little")
  writeBin(as.integer(exceptions$flags), connection, size = 1L, endian = "little")
  rawConnectionValue(connection)
}

pcodec_native_read_exception_bytes <- function(blob, count) {
  if (!count) {
    return(data.frame(row = integer(), z = numeric(), log2se = numeric(),
                      eaf = numeric(), flags = integer()))
  }
  if (length(blob) != count * 17L) stop("native Pcodec exception frame is truncated", call. = FALSE)
  connection <- rawConnection(blob, open = "rb")
  on.exit(close(connection), add = TRUE)
  data.frame(
    row = readBin(connection, integer(), n = count, size = 4L, endian = "little"),
    z = readBin(connection, numeric(), n = count, size = 4L, endian = "little"),
    log2se = readBin(connection, numeric(), n = count, size = 4L, endian = "little"),
    eaf = readBin(connection, numeric(), n = count, size = 4L, endian = "little"),
    flags = readBin(connection, integer(), n = count, size = 1L, endian = "little")
  )
}

pcodec_native_write_exceptions <- function(exceptions, output, blocks) {
  path <- file.path(output, "exceptions.bin")
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  offset <- 0
  locations <- vector("list", length(blocks))
  for (block in seq_along(blocks)) {
    inside <- exceptions$row >= blocks[[block]]$row_start &
      exceptions$row < blocks[[block]]$row_stop
    raw_blob <- pcodec_native_exception_bytes(exceptions[inside, , drop = FALSE])
    blob <- if (length(raw_blob)) pcodec_native_zstd_compress(raw_blob, level = 19L) else raw()
    if (length(blob)) writeBin(blob, connection, useBytes = TRUE)
    locations[[block]] <- list(
      offset = offset, length = length(blob), raw_length = length(raw_blob),
      count = sum(inside)
    )
    offset <- offset + length(blob)
  }
  list(file = basename(path), bytes = offset, codec = "zstd", record_bytes = 17L,
       blocks = locations)
}

pcodec_native_write_store <- function(data, output, metadata = list()) {
  if (!pcodec_native_enabled()) {
    stop("native Pcodec is not enabled", call. = FALSE)
  }
  required <- c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error",
                "effect_allele_frequency", "z")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Pcodec input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  identity <- pcodec_native_identity(data)
  order <- order(identity$global_position, identity$substitution, method = "radix")
  ordered_position <- identity$global_position[order]
  ordered_substitution <- identity$substitution[order]
  if (anyDuplicated(paste(ordered_position, ordered_substitution, sep = ":"))) {
    stop("duplicate full REF/ALT identity keys", call. = FALSE)
  }
  ordered <- data[order, , drop = FALSE]
  values <- pcodec_native_quantise(ordered)
  n <- nrow(ordered)
  block_rows <- pcodec_native_validate_block_rows(
    metadata$block_rows %||% PCODEC_NATIVE_BLOCK_ROWS, "native Pcodec block_rows")
  key_block_rows <- pcodec_native_validate_block_rows(
    metadata$key_block_rows %||% PCODEC_NATIVE_KEY_BLOCK_ROWS,
    "native Pcodec key_block_rows")
  block_template <- pcodec_native_block_template(ordered_position, block_rows)
  key_block_template <- pcodec_native_block_template(ordered_position, key_block_rows)
  streams <- list(
    position = pcodec_native_append_stream(
      pcodec_native_position_gaps(ordered_position, key_block_rows),
      file.path(output, "position.pco"), "u32", key_block_rows),
    substitution = pcodec_native_append_stream(
      ordered_substitution, file.path(output, "substitution.pco"), "u8", key_block_rows),
    z = pcodec_native_append_stream(
      values$z, file.path(output, "z.pco"), "u16", block_rows),
    eaf = pcodec_native_append_stream(
      values$eaf, file.path(output, "eaf.pco"), "u8", block_rows),
    se = pcodec_native_append_stream(
      values$se, file.path(output, "se.pco"), "u8", block_rows)
  )
  exception_stream <- pcodec_native_write_exceptions(values$exceptions, output, block_template)
  index <- list(
    format = "CompreSSoR-native-index", version = 3L,
    position_encoding = "delta_u32_within_block",
    rows = n, block_rows = block_rows,
    key_block_rows = key_block_rows, value_block_rows = block_rows,
    blocks = block_template, key_blocks = key_block_template,
    value_blocks = block_template,
    streams = streams, exceptions = exception_stream
  )
  jsonlite::write_json(index, file.path(output, "native.index.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 17)
  files <- list(
    position = "position.pco", substitution = "substitution.pco",
    z = "z.pco", eaf = "eaf.pco", se = "se.pco",
    exceptions = "exceptions.bin", index = "native.index.json"
  )
  chromosomes <- compressor_grch38_chromosome_lengths
  offsets <- pcodec_native_offsets()
  manifest <- list(
    format = "CompreSSoR", format_version = PCODEC_NATIVE_FORMAT,
    backend = "pcodec", profile = "standard", rows = n, n_rows = n,
    block_rows = block_rows, key_block_rows = key_block_rows,
    value_block_rows = block_rows,
    identity = list(
      encoding = "native_global_position_plus_full_ref_alt_code",
      position_storage = "within-block delta-coded uint32; block first positions are in the index",
      external_reference_required = FALSE,
      chromosome_lengths = as.list(chromosomes),
      chromosome_offsets = as.list(offsets),
      effect_allele_is_alt = TRUE, other_allele_is_ref = TRUE
    ),
    semantic_codec = list(
      name = "z9/eaf8/se8", z_bits = 9L, eaf_bits = 8L,
      se_bits = PCODEC_NATIVE_SE_BITS, z_range = c(-3.5, 3.5),
      se_count = PCODEC_NATIVE_SE_COUNT,
      se_residual_range = PCODEC_NATIVE_SE_RESIDUAL_RANGE,
      se_missing = PCODEC_NATIVE_SE_COUNT,
      se_exception = PCODEC_NATIVE_SE_COUNT + 1L,
      se_center_block_rows = PCODEC_NATIVE_SE_CENTER_ROWS,
      block_centers_log2_residual = values$centres,
      exception_rows = nrow(values$exceptions), exception_precision = "float32",
      beta = "derived as z * standard_error",
      p_value = "derived as 2 * pnorm(-abs(z))"
    ),
    codec = list(
      name = "pcodec_native_standalone_z9_eaf8_se8_zstd_exceptions",
      library = "pcodec", pco_version = "1.0.3", abi = "standalone",
      page_rows = PCODEC_NATIVE_PAGE_ROWS,
      compression = paste0("Pcodec standalone streams; ", key_block_rows,
                           "-row key frames and ", block_rows, "-row value frames"),
      z_bits = 9L, eaf_bits = 8L, se_bits = PCODEC_NATIVE_SE_BITS,
      se_residual_range = PCODEC_NATIVE_SE_RESIDUAL_RANGE,
      p_storage = "omitted; derived from z",
      beta_storage = "omitted; derived from z and standard_error",
      exception_storage = "zstd level 19, 17-byte float32 records"
    ),
    files = files, logical_columns = c("z", "standard_error", "effect_allele_frequency"),
    derived_columns = list(beta = "z * standard_error", p_value = "2 * pnorm(-abs(z))"),
    variant_storage = "self_contained_identity_key",
    variant_identity = list(encoding = "global_position_plus_full_ref_alt_code",
                            position_storage = "within-block delta-coded uint32",
                            rsid = "not_stored", external_reference_required = FALSE),
    genome_build = "GRCh38",
    reference = metadata$reference %||% list(id = "none", build = "GRCh38",
                                             status = "identity key is self-contained"),
    harmonisation = metadata$harmonisation %||% list(method = "not recorded"),
    source_columns = metadata$source_columns %||% names(data),
    source = metadata$source %||% NULL,
    metadata = metadata
  )
  manifest$integrity <- list(
    algorithm = "sha256",
    files = stats::setNames(lapply(files, function(relative) {
      path <- file.path(output, relative)
      list(bytes = as.numeric(file.info(path)$size),
           sha256 = digest::digest(path, algo = "sha256", file = TRUE))
    }), unname(unlist(files)))
  )
  manifest$created_utc <- now_utc()
manifest$tolerances <- list(
    eaf_abs_max = 0.004, z_abs_max_central = 7 / (2 * (2^9 - 2)),
    se_relative_max = 2^(4 / PCODEC_NATIVE_SE_COUNT) - 1,
    beta_error_bound = "1.02 * (abs(SE) * z_abs_max_central + abs(Z) * abs(SE) * se_relative_max)",
    z_central_range = c(-3.5, 3.5),
    se_profile = "block-centred log2 residual quantisation",
    exception_precision = "float32", exact_values_in_exception_sidecar = FALSE
  )
  write_manifest(manifest, file.path(output, "manifest.json"))
  seal_pcodec_manifest(file.path(output, "manifest.json"))
  invisible(manifest)
}

pcodec_native_read_index <- function(store) {
  index_path <- file.path(store$path, store$manifest$files$index)
  index <- jsonlite::fromJSON(index_path, simplifyVector = FALSE)
  if (!identical(index$format, "CompreSSoR-native-index") ||
      !as.integer(index$version) %in% c(1L, 2L, 3L)) {
    stop("invalid native Pcodec index", call. = FALSE)
  }
  index
}

pcodec_native_read_blob <- function(path, offset, length) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  seek(connection, where = as.numeric(offset), origin = "start")
  blob <- readBin(connection, raw(), n = as.integer(length), endian = "little")
  if (length(blob) != as.integer(length)) stop("native Pcodec stream is truncated", call. = FALSE)
  blob
}

pcodec_native_read_stream_block <- function(store, index, stream, block) {
  spec <- index$streams[[stream]]
  location <- spec$blocks[[block]]
  blob <- pcodec_native_read_blob(
    file.path(store$path, spec$file), location$offset, location$length
  )
  values <- pcodec_native_decompress(blob, as.integer(location$values),
                           if (stream %in% c("position")) "u32" else
                             if (stream %in% c("z")) "u16" else "u8")
  if (identical(stream, "position") && identical(index$position_encoding, "delta_u32_within_block")) {
    key_blocks <- pcodec_native_index_blocks(index, "key")
    values <- cumsum(values) + as.numeric(key_blocks[[block]]$first_position)
  }
  values
}

pcodec_native_read_stream_all <- function(store, index, stream) {
  spec <- index$streams[[stream]]
  connection <- file(file.path(store$path, spec$file), open = "rb")
  on.exit(close(connection), add = TRUE)
  parts <- lapply(seq_along(spec$blocks), function(block) {
    location <- spec$blocks[[block]]
    seek(connection, where = as.numeric(location$offset), origin = "start")
    blob <- readBin(connection, raw(), n = as.integer(location$length), endian = "little")
    if (length(blob) != as.integer(location$length)) {
      stop("native Pcodec stream is truncated", call. = FALSE)
    }
    values <- pcodec_native_decompress(blob, as.integer(location$values),
      if (stream == "position") "u32" else if (stream == "z") "u16" else "u8")
    if (identical(stream, "position") && identical(index$position_encoding, "delta_u32_within_block")) {
      key_blocks <- pcodec_native_index_blocks(index, "key")
      values <- cumsum(values) + as.numeric(key_blocks[[block]]$first_position)
    }
    values
  })
  if (!length(parts)) {
    return(if (stream == "position") numeric() else integer())
  }
  unlist(parts, use.names = FALSE)
}

pcodec_native_block_matrix <- function(blocks, stream_blocks = blocks,
                                       first_position = FALSE) {
  n <- length(blocks)
  if (length(stream_blocks) != n) {
    stop("native Pcodec index and stream block counts differ", call. = FALSE)
  }
  columns <- if (first_position) 6L else 5L
  if (!n) return(matrix(numeric(), nrow = 0L, ncol = columns))
  values <- Map(function(block, stream_block) {
    out <- c(as.numeric(stream_block$offset), as.numeric(stream_block$length),
             as.numeric(stream_block$values), as.numeric(block$row_start),
             as.numeric(block$row_stop))
    if (first_position) out <- c(out, as.numeric(block$first_position))
    out
  }, blocks, stream_blocks)
  matrix(unlist(values, use.names = FALSE), nrow = n, ncol = columns,
         byrow = TRUE)
}

pcodec_native_read_native_codes <- function(store, index, streams, threads = 1L) {
  if (!is.loaded("compressor_read_pcodec_native_codes", PACKAGE = "CompreSSoR")) {
    stop("native Pcodec stream reader is not available in this build", call. = FALSE)
  }
  key_blocks <- pcodec_native_index_blocks(index, "key")
  value_blocks <- pcodec_native_index_blocks(index, "value")
  exception_blocks <- index$exceptions$blocks %||% list()
  files <- c(
    file.path(store$path, index$streams$position$file),
    file.path(store$path, index$streams$substitution$file),
    file.path(store$path, index$streams$z$file),
    file.path(store$path, index$streams$eaf$file),
    file.path(store$path, index$streams$se$file),
    file.path(store$path, index$exceptions$file)
  )
  .Call("compressor_read_pcodec_native_codes", files,
        pcodec_native_block_matrix(key_blocks, index$streams$position$blocks,
                                   first_position = TRUE),
        pcodec_native_block_matrix(key_blocks, index$streams$substitution$blocks,
                                   first_position = TRUE),
        pcodec_native_block_matrix(value_blocks, index$streams$z$blocks),
        pcodec_native_block_matrix(value_blocks, index$streams$eaf$blocks),
        pcodec_native_block_matrix(value_blocks, index$streams$se$blocks),
        do.call(rbind, lapply(exception_blocks, function(block) {
          c(as.numeric(block$offset), as.numeric(block$length),
            as.numeric(block$count), as.numeric(block$raw_length))
        })) %||% matrix(numeric(), nrow = 0L, ncol = 4L),
        as.numeric(store$manifest$n_rows %||% store$manifest$rows),
        as.character(streams), index$exceptions$codec %||% "raw",
        as.integer(threads),
        PACKAGE = "CompreSSoR")
}

pcodec_native_read_all_exceptions <- function(store, index) {
  locations <- index$exceptions$blocks
  if (!length(locations)) {
    return(data.frame(row = integer(), z = numeric(), log2se = numeric(),
                      eaf = numeric(), flags = integer()))
  }
  parts <- lapply(seq_along(locations), function(block) {
    pcodec_native_read_exception_block(store, index, block)
  })
  parts <- parts[vapply(parts, nrow, integer(1)) > 0L]
  if (!length(parts)) {
    return(data.frame(row = integer(), z = numeric(), log2se = numeric(),
                      eaf = numeric(), flags = integer()))
  }
  do.call(rbind, parts)
}

pcodec_native_full_read <- function(store, index, requested, need_identity,
                                     need_z, need_se, need_eaf, threads = 1L) {
  threads <- pcodec_validate_threads(threads)
  n <- as.integer(store$manifest$n_rows %||% store$manifest$rows)
  semantic <- store$manifest$semantic_codec %||% list()
  z_count <- as.integer(semantic$z_count %||% 510L)
  se_count <- as.integer(semantic$se_count %||% 62L)
  eaf_count <- as.integer(semantic$eaf_count %||% 255L)
  z_bits <- as.integer(semantic$z_bits %||% 9L)
  se_bits <- as.integer(semantic$se_bits %||% 6L)
  eaf_bits <- as.integer(semantic$eaf_bits %||% 8L)
  z_range <- as.numeric(unlist(semantic$z_range %||% c(-3.5, 3.5)))
  se_range <- as.numeric(unlist(semantic$se_residual_range %||% c(-1, 1)))
  needed <- unique(c(if (need_z) "z", if (need_se) "se", if (need_eaf) "eaf"))
  if (need_se) needed <- unique(c(needed, "eaf"))
  streams <- unique(c(needed, if (need_identity) c("position", "substitution")))
  native_stream_reader <- isTRUE(getOption(
    "CompreSSoR.pcodec.native_stream_reader", TRUE
  )) && is.loaded("compressor_read_pcodec_native_codes", PACKAGE = "CompreSSoR")
  if (native_stream_reader) {
    codes <- pcodec_native_read_native_codes(store, index, streams, threads = threads)
  } else {
    decoded_streams <- pcodec_parallel_lapply(
      streams, function(stream) pcodec_native_read_stream_all(store, index, stream),
      threads = threads
    )
    codes <- stats::setNames(decoded_streams, streams)
  }
  if ("z" %in% needed) z_code <- codes$z else z_code <- rep.int(z_count, n)
  if ("se" %in% needed) se_code <- codes$se else se_code <- rep.int(se_count, n)
  if ("eaf" %in% needed) eaf_code <- codes$eaf else eaf_code <- rep.int(0L, n)
  exceptions <- if (length(needed)) {
    if (native_stream_reader) {
      native_exceptions <- codes$exceptions
      data.frame(row = native_exceptions$row, z = native_exceptions$z,
                 log2se = native_exceptions$log2se,
                 eaf = native_exceptions$eaf,
                 flags = native_exceptions$flags)
    } else {
      pcodec_native_read_all_exceptions(store, index)
    }
  } else {
    data.frame(row = integer(), z = numeric(), log2se = numeric(),
               eaf = numeric(), flags = integer())
  }
  decoded <- if (length(needed)) {
    .Call("compressor_decode_native", as.integer(z_code), as.integer(se_code),
      as.integer(eaf_code), z_range[1], z_range[2], z_count, se_count, eaf_count,
      z_bits, se_bits, eaf_bits,
      as.integer(semantic$se_center_block_rows %||%
                   PCODEC_NATIVE_SE_CENTER_ROWS),
      as.numeric(unlist(semantic$block_centers_log2_residual)),
      as.integer(exceptions$row), as.numeric(exceptions$z),
      as.numeric(2^exceptions$log2se), as.numeric(exceptions$eaf),
      as.integer(exceptions$flags), isTRUE("beta" %in% requested),
      isTRUE("p_value" %in% requested), se_range[1], se_range[2],
      PACKAGE = "CompreSSoR")
  } else list(z = numeric(n), standard_error = numeric(n),
              effect_allele_frequency = numeric(n))
  output <- data.frame(row = seq_len(n) - 1L, stringsAsFactors = FALSE)
  if (need_identity) {
    position <- codes$position
    substitution <- codes$substitution
    if (any(c("global_position", "substitution") %in% requested)) {
      output$global_position <- position
      output$substitution <- as.integer(substitution)
    }
    if (is.null(requested) || any(c("chromosome", "base_pair_location",
                                    "effect_allele", "other_allele") %in% requested)) {
      output <- cbind(output, as.data.frame(pcodec_native_key_columns(position, substitution),
                                            stringsAsFactors = FALSE))
    }
  }
  if (need_z) output$z <- decoded$z
  if (need_se) output$standard_error <- decoded$standard_error
  if (need_eaf) output$effect_allele_frequency <- decoded$effect_allele_frequency
  if ("beta" %in% requested) output$beta <- decoded$beta
  if ("p_value" %in% requested) output$p_value <- decoded$p_value
  output[c("row", if (is.null(requested)) names(output)[-1L] else requested)]
}

pcodec_native_read_exception_block <- function(store, index, block) {
  location <- index$exceptions$blocks[[block]]
  if (!as.integer(location$count)) {
    return(data.frame(row = integer(), z = numeric(), log2se = numeric(),
                      eaf = numeric(), flags = integer()))
  }
  blob <- pcodec_native_read_blob(
    file.path(store$path, index$exceptions$file), location$offset, location$length
  )
  codec <- index$exceptions$codec %||% "raw"
  if (identical(codec, "zstd")) {
    raw_length <- as.integer(location$raw_length %||% (as.integer(location$count) * 17L))
    blob <- pcodec_native_zstd_decompress(blob, raw_length)
  } else if (!identical(codec, "raw")) {
    stop("unsupported native exception codec: ", codec, call. = FALSE)
  }
  pcodec_native_read_exception_bytes(blob, as.integer(location$count))
}

pcodec_native_block_ids_for_rows <- function(index, rows, kind = "value") {
  if (!length(rows)) return(integer())
  blocks <- pcodec_native_index_blocks(index, kind)
  stops <- vapply(blocks, function(block) as.numeric(block$row_stop), numeric(1))
  unique(findInterval(as.numeric(rows), stops) + 1L)
}

pcodec_native_key_columns <- function(position, substitution) {
  lengths <- compressor_grch38_chromosome_lengths
  offsets <- pcodec_native_offsets()
  chromosome_code <- findInterval(position, offsets)
  chromosome_code <- pmax(1L, pmin(length(lengths), chromosome_code))
  list(
    chromosome = names(lengths)[chromosome_code],
    base_pair_location = as.integer(position - offsets[chromosome_code] + 1),
    effect_allele = c("A", "C", "G", "T")[(bitwAnd(as.integer(substitution), 3L)) + 1L],
    other_allele = c("A", "C", "G", "T")[bitwShiftR(as.integer(substitution), 2L) + 1L]
  )
}

pcodec_native_target_keys <- function(variants) {
  parsed <- parse_canonical_variant_keys(variants)
  if (anyNA(parsed$chromosome) || anyNA(parsed$base_pair_location) ||
      anyNA(parsed$reference_allele) || anyNA(parsed$alternate_allele)) {
    stop("invalid canonical variant key", call. = FALSE)
  }
  canonical <- compressor_variant_key(
    parsed$chromosome, parsed$base_pair_location,
    parsed$reference_allele, parsed$alternate_allele
  )
  identity <- pcodec_native_identity(data.frame(
    chromosome = parsed$chromosome,
    base_pair_location = parsed$base_pair_location,
    other_allele = parsed$reference_allele,
    effect_allele = parsed$alternate_allele,
    stringsAsFactors = FALSE
  ))
  list(canonical = canonical, position = identity$global_position,
       substitution = identity$substitution)
}

pcodec_native_decode_values <- function(codes, exceptions, centre_id, centres,
                                         row_start, n, needed,
                                         semantic_codec = list()) {
  output <- list()
  z_count <- as.integer(semantic_codec$z_count %||% 510L)
  eaf_count <- as.integer(semantic_codec$eaf_count %||% 255L)
  se_count <- as.integer(semantic_codec$se_count %||% 62L)
  se_range <- as.numeric(unlist(semantic_codec$se_residual_range %||% c(-1, 1)))
  if (length(se_range) != 2L || !all(is.finite(se_range)) || se_range[2] <= se_range[1]) {
    stop("invalid native semantic SE residual range", call. = FALSE)
  }
  if ("z" %in% needed) {
    z <- rep(NA_real_, n)
    ok <- codes$z < z_count
    z[ok] <- -3.5 + (codes$z[ok] + 0.5) * (7 / z_count)
    output$z <- z
  }
  if ("eaf" %in% needed || "se" %in% needed) {
    output$eaf <- sin(pi * as.numeric(codes$eaf) / (2 * eaf_count))^2
  }
  if ("se" %in% needed) {
    se <- rep(NA_real_, n)
    ok <- codes$se < se_count
    safe <- pmin(1 - 1e-12, pmax(1e-12, output$eaf))
    residual <- se_range[1] + (codes$se + 0.5) * diff(se_range) / se_count
    se[ok] <- 2^(residual[ok] + centres[centre_id] -
      0.5 * log2(2 * safe[ok] * (1 - safe[ok])))
    output$se <- se
  }
  if (nrow(exceptions)) {
    local <- as.integer(exceptions$row) - as.integer(row_start) + 1L
    valid <- local >= 1L & local <= n
    local <- local[valid]
    exceptions <- exceptions[valid, , drop = FALSE]
    if ("z" %in% needed) {
      keep <- bitwAnd(exceptions$flags, 1L) != 0L
      output$z[local[keep]] <- exceptions$z[keep]
    }
    if ("se" %in% needed) {
      keep <- bitwAnd(exceptions$flags, 2L) != 0L
      output$se[local[keep]] <- 2^exceptions$log2se[keep]
    }
    if ("eaf" %in% needed) {
      keep <- bitwAnd(exceptions$flags, 4L) != 0L
      output$eaf[local[keep]] <- exceptions$eaf[keep]
    }
  }
  output
}

pcodec_native_empty_result <- function(columns) {
  result <- data.frame(row = integer(), stringsAsFactors = FALSE)
  all_columns <- c("global_position", "substitution", "chromosome", "base_pair_location", "effect_allele",
                   "other_allele", "z", "beta", "standard_error",
                   "effect_allele_frequency", "p_value")
  needed <- if (is.null(columns)) all_columns else unique(columns)
  for (column in setdiff(needed, names(result))) {
    result[[column]] <- switch(column,
      chromosome = character(), effect_allele = character(),
      other_allele = character(), base_pair_location = integer(),
      substitution = integer(), global_position = numeric(),
      numeric())
  }
  result[c("row", needed)]
}

pcodec_native_read_store <- function(store, region = NULL, variants = NULL,
                                      columns = NULL, threads = NULL) {
  if (!pcodec_native_available()) {
    stop("native Pcodec is not available in this build", call. = FALSE)
  }
  threads <- pcodec_native_default_threads(region = region, variants = variants,
                                            threads = threads)
  manifest <- store$manifest
  index <- pcodec_native_read_index(store)
  n <- as.integer(manifest$n_rows %||% manifest$rows)
  if (!is.null(columns) && !length(columns)) stop("columns must contain at least one column name", call. = FALSE)
  requested <- if (is.null(columns)) {
    c("chromosome", "base_pair_location", "effect_allele", "other_allele",
      "z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  } else unique(as.character(columns))
  allowed <- c("global_position", "substitution", "chromosome", "base_pair_location", "effect_allele", "other_allele",
               "z", "beta", "standard_error", "effect_allele_frequency", "p_value")
  unknown <- setdiff(requested, allowed)
  if (length(unknown)) stop("unknown output column(s): ", paste(unknown, collapse = ", "), call. = FALSE)

  row_targets <- NULL
  key_targets <- NULL
  if (!is.null(variants)) {
    if (is.character(variants)) {
      key_targets <- pcodec_native_target_keys(unique(trimws(variants)))
    } else {
      row_targets <- unique(as.integer(variants))
      if (anyNA(row_targets) || any(row_targets < 0L | row_targets >= n)) {
        stop("variants must be valid zero-based row IDs", call. = FALSE)
      }
    }
  }
  identity_needed <- is.null(columns) || any(c("global_position", "substitution",
                                                "chromosome", "base_pair_location",
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
    if (!chromosome %in% names(compressor_grch38_chromosome_lengths)) stop("unsupported region chromosome", call. = FALSE)
    offsets <- pcodec_native_offsets()
    lower <- offsets[match(chromosome, names(compressor_grch38_chromosome_lengths))] + bounds$start - 1
    upper <- offsets[match(chromosome, names(compressor_grch38_chromosome_lengths))] + bounds$end - 1
  }
  if (is.null(region) && is.null(variants)) {
    output <- pcodec_native_full_read(store, index, requested, identity_needed,
                                      need_z, need_se, need_eaf, threads = threads)
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
        selected_position[selected_index], selected_substitution[selected_index])
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

pcodec_native_validate_store <- function(store, full = FALSE) {
  result <- tryCatch({
    s <- if (inherits(store, "compressor_store")) store else open_compressor(store)
    m <- s$manifest
    index <- pcodec_native_read_index(s)
    errors <- character()
    files <- unname(unlist(m$files))
    missing <- files[!file.exists(file.path(s$path, files))]
    if (length(missing)) errors <- c(errors, paste("missing", missing))
    blocks <- index$blocks
    if (length(blocks)) {
      starts <- vapply(blocks, function(block) as.integer(block$row_start), integer(1))
      stops <- vapply(blocks, function(block) as.integer(block$row_stop), integer(1))
      if (starts[1] != 0L || utils::tail(stops, 1) != as.integer(m$n_rows) ||
          any(starts[-1] != stops[-length(stops)])) {
        errors <- c(errors, "native blocks do not cover rows contiguously")
      }
    }
    if (isTRUE(full) && !length(errors)) {
      pcodec_native_read_store(s, columns = c("z", "standard_error",
                                               "effect_allele_frequency"))
    }
    list(valid = !length(errors), errors = errors,
         rows = as.integer(m$n_rows), profile = m$profile, full = isTRUE(full))
  }, error = function(e) list(valid = FALSE, errors = conditionMessage(e),
                               rows = NA_integer_, profile = NA_character_, full = isTRUE(full)))
  result
}
