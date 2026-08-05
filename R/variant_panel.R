## Compact native Pcodec variant panel.
##
## A panel is deliberately not represented as a fake summary-statistics store.
## It contains the exact build-aware identity streams used by the native GWAS
## store and one lossless uint8 HM3 membership stream.

PCODEC_VARIANT_PANEL_FORMAT <- "0.1.0-pcodec-native-panel"
PCODEC_VARIANT_PANEL_INDEX_FORMAT <- "CompreSSoR-native-variant-panel-index"
PCODEC_VARIANT_PANEL_DEFAULT <- "core_hm3.cpr"

is_variant_panel_path <- function(path) {
  if (!dir.exists(path) || !file.exists(file.path(path, "manifest.json"))) return(FALSE)
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(FALSE)
  manifest <- tryCatch(
    jsonlite::read_json(file.path(path, "manifest.json"), simplifyVector = FALSE),
    error = function(e) NULL
  )
  is.list(manifest) && identical(manifest$format, "CompreSSoR-variant-panel")
}

bundled_variant_panel_path <- function(name = "core") {
  if (!name %in% c("core", "hm3")) return("")
  candidates <- c(
    system.file("extdata", "panels", PCODEC_VARIANT_PANEL_DEFAULT,
                package = "CompreSSoR"),
    file.path("inst", "extdata", "panels", PCODEC_VARIANT_PANEL_DEFAULT),
    file.path(".", "inst", "extdata", "panels", PCODEC_VARIANT_PANEL_DEFAULT)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[dir.exists(candidates)]
  if (length(hits)) normalizePath(hits[[1L]], mustWork = FALSE) else ""
}

variant_panel_manifest <- function(path) {
  manifest_path <- file.path(path, "manifest.json")
  if (!file.exists(manifest_path) || !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("variant panel is missing manifest.json: ", path, call. = FALSE)
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  if (!is.list(manifest) || !identical(manifest$format, "CompreSSoR-variant-panel") ||
      !identical(manifest$backend, "pcodec")) {
    stop("not a native CompreSSoR variant panel: ", path, call. = FALSE)
  }
  manifest
}

variant_panel_files <- function(manifest) {
  files <- manifest$files
  if (!is.list(files) || !length(files)) stop("variant panel manifest has no files", call. = FALSE)
  unname(unlist(files, use.names = FALSE))
}

variant_panel_index <- function(panel) {
  manifest <- panel$manifest
  index_path <- file.path(panel$path, manifest$files$index)
  if (!file.exists(index_path)) stop("variant panel index is missing", call. = FALSE)
  index <- jsonlite::read_json(index_path, simplifyVector = FALSE)
  if (!is.list(index) || !identical(index$format, PCODEC_VARIANT_PANEL_INDEX_FORMAT)) {
    stop("invalid native variant-panel index", call. = FALSE)
  }
  n <- as.numeric(manifest$n_rows %||% manifest$rows)
  if (length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n)) {
    stop("variant panel manifest row count is invalid", call. = FALSE)
  }
  pcodec_native_validate_index_partition(index$key_blocks, n, "variant-panel key blocks")
  pcodec_native_validate_index_partition(index$flag_blocks, n, "variant-panel flag blocks")
  same_rows <- function(left, right, label) {
    if (length(left) != length(right) || length(left) != 0L && any(vapply(
        seq_along(left), function(i) {
          as.numeric(left[[i]]$row_start) == as.numeric(right[[i]]$row_start) &&
            as.numeric(left[[i]]$row_stop) == as.numeric(right[[i]]$row_stop)
        }, logical(1)) == FALSE)) {
      stop(label, " are not row-aligned", call. = FALSE)
    }
  }
  same_rows(index$key_blocks, index$streams$position$blocks,
            "variant-panel position blocks")
  same_rows(index$key_blocks, index$streams$substitution$blocks,
            "variant-panel substitution blocks")
  same_rows(index$flag_blocks, index$streams$hm3$blocks,
            "variant-panel HM3 blocks")
  index
}

variant_panel_payload_integrity <- function(path, manifest) {
  files <- manifest$integrity$files
  if (!is.list(files) || !length(files)) stop("variant panel payload checksums are missing", call. = FALSE)
  expected_names <- sort(names(files))
  actual_names <- sort(variant_panel_files(manifest))
  if (!identical(expected_names, actual_names)) {
    stop("variant panel integrity file list does not match manifest files", call. = FALSE)
  }
  # JSON objects are represented as named lists whose names are the relative
  # paths in the manifest. Keep the loop explicit so a missing file fails with
  # its exact relative path.
  observed <- lapply(names(files), function(relative) {
    file <- file.path(path, relative)
    if (!file.exists(file) || dir.exists(file)) stop("variant panel payload is missing: ", relative, call. = FALSE)
    list(bytes = as.numeric(file.info(file)$size),
         sha256 = digest::digest(file, algo = "sha256", file = TRUE))
  })
  names(observed) <- names(files)
  for (relative in names(files)) {
    expected <- files[[relative]]
    observed_item <- observed[[relative]]
    if (!identical(as.numeric(expected$bytes), observed_item$bytes) ||
        !identical(tolower(as.character(expected$sha256)),
                   tolower(observed_item$sha256))) {
      stop("variant panel payload checksum mismatch: ", relative, call. = FALSE)
    }
  }
  payload_sha256 <- pcodec_payload_sha256(observed)
  if (!identical(tolower(as.character(manifest$integrity$payload_sha256)),
                 tolower(payload_sha256))) {
    stop("variant panel payload aggregate checksum mismatch", call. = FALSE)
  }
  invisible(observed)
}

validate_variant_panel <- function(panel, full = FALSE) {
  path <- if (inherits(panel, "compressor_variant_panel")) panel$path else panel
  if (length(path) != 1L || !is.character(path) || !dir.exists(path)) {
    stop("panel must be a native variant-panel directory", call. = FALSE)
  }
  manifest_path <- file.path(path, "manifest.json")
  manifest <- variant_panel_manifest(path)
  verify_pcodec_manifest(manifest_path)
  index <- variant_panel_index(list(path = path, manifest = manifest))
  files <- variant_panel_files(manifest)
  if (any(!file.exists(file.path(path, files)))) {
    stop("variant panel is missing one or more payload files", call. = FALSE)
  }
  if (isTRUE(full)) {
    variant_panel_payload_integrity(path, manifest)
    store <- list(path = path, manifest = manifest)
    position <- pcodec_native_read_stream_all(store, index, "position")
    substitution <- pcodec_native_read_stream_all(store, index, "substitution")
    hm3 <- pcodec_native_read_stream_all(store, index, "hm3")
    if (length(position) != manifest$n_rows || length(substitution) != manifest$n_rows ||
        length(hm3) != manifest$n_rows) {
      stop("variant panel stream lengths do not match the manifest", call. = FALSE)
    }
    if (length(position) > 1L && any(diff(position) < 0)) {
      stop("variant panel positions are not sorted", call. = FALSE)
    }
    if (length(substitution) && any(substitution < 0L | substitution > 15L |
        bitwShiftR(as.integer(substitution), 2L) == bitwAnd(as.integer(substitution), 3L))) {
      stop("variant panel substitution stream contains an invalid code", call. = FALSE)
    }
    if (length(hm3) && any(!hm3 %in% 0:1)) {
      stop("variant panel HM3 stream contains a non-binary flag", call. = FALSE)
    }
    if (!identical(as.integer(sum(hm3)), as.integer(manifest$hm3_rows))) {
      stop("variant panel HM3 count does not match the manifest", call. = FALSE)
    }
  }
  invisible(list(valid = TRUE, manifest = manifest, index = index,
                 payload_verified = isTRUE(full)))
}

open_variant_panel <- function(path, verify_payload = FALSE) {
  path <- normalizePath(path, mustWork = TRUE)
  manifest <- variant_panel_manifest(path)
  validate_variant_panel(path, full = verify_payload)
  structure(list(path = path, manifest = manifest),
            class = "compressor_variant_panel")
}

print.compressor_variant_panel <- function(x, ...) {
  m <- x$manifest
  cat("CompreSSoR native variant panel\n")
  cat("Path:       ", x$path, "\n", sep = "")
  cat("Universe:   ", m$panel_universe %||% "core", "\n", sep = "")
  cat("Build:      ", m$genome_build %||% m$build, "\n", sep = "")
  cat("Rows:       ", m$n_rows %||% m$rows, "\n", sep = "")
  cat("HM3 rows:   ", m$hm3_rows, "\n", sep = "")
  invisible(x)
}

write_variant_panel <- function(data, output, build = "GRCh38",
                                 source_provenance = NULL,
                                 block_rows = PCODEC_NATIVE_BLOCK_ROWS,
                                 key_block_rows = PCODEC_NATIVE_KEY_BLOCK_ROWS,
                                 overwrite = FALSE) {
  if (!pcodec_native_enabled()) {
    stop("native Pcodec is not enabled; install Rust/Cargo and reinstall CompreSSoR",
         call. = FALSE)
  }
  if (!is.data.frame(data) || !nrow(data)) stop("variant panel data must be non-empty", call. = FALSE)
  if (length(output) != 1L || !is.character(output) || !nzchar(output)) {
    stop("output must be one non-empty path", call. = FALSE)
  }
  build <- compressor_normalize_build(build)
  block_rows <- pcodec_native_validate_block_rows(block_rows, "variant-panel block_rows")
  key_block_rows <- pcodec_native_validate_block_rows(key_block_rows,
                                                       "variant-panel key_block_rows")
  if (file.exists(output) && !isTRUE(overwrite)) {
    stop("variant panel output already exists; set overwrite=TRUE: ", output, call. = FALSE)
  }
  if (!"hm3" %in% names(data)) {
    data$hm3 <- 0L
  }
  if (anyNA(data$hm3) || any(!as.numeric(data$hm3) %in% 0:1)) {
    stop("variant panel hm3 must be a lossless 0/1 flag", call. = FALSE)
  }
  # Keep only the identity and flag columns before normalisation. In particular
  # this prevents panel annotations or rsIDs from entering the artifact.
  identity_input <- data[, intersect(c("variant_id", "chromosome", "base_pair_location",
                                       "other_allele", "effect_allele"), names(data)),
                         drop = FALSE]
  identity_input$hm3 <- as.integer(data$hm3)
  normalized <- normalise_variant_set_columns(identity_input, build = build)
  if (!isTRUE(attr(normalized, "variant_set_canonical"))) {
    stop("variant panel identity must be canonical chromosome:position:REF:ALT", call. = FALSE)
  }
  required <- c("chromosome", "base_pair_location", "other_allele", "effect_allele")
  if (any(!vapply(normalized[required], function(x) !anyNA(x), logical(1)))) {
    stop("variant panel contains incomplete canonical identity", call. = FALSE)
  }
  identity <- pcodec_native_identity(normalized, build = build)
  order <- order(identity$global_position, identity$substitution, method = "radix")
  ordered_position <- identity$global_position[order]
  ordered_substitution <- identity$substitution[order]
  ordered_hm3 <- as.integer(normalized$hm3[order])
  if (anyDuplicated(paste(ordered_position, ordered_substitution, sep = ":"))) {
    stop("variant panel contains duplicate canonical identity keys", call. = FALSE)
  }
  parent <- dirname(normalizePath(output, mustWork = FALSE))
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile("compressor-variant-panel-", tmpdir = parent)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  streams <- list(
    position = pcodec_native_append_stream(
      pcodec_native_position_gaps(ordered_position, key_block_rows),
      file.path(staging, "position.pco"), "u32", key_block_rows
    ),
    substitution = pcodec_native_append_stream(
      ordered_substitution, file.path(staging, "substitution.pco"), "u8", key_block_rows
    ),
    hm3 = pcodec_native_append_stream(
      ordered_hm3, file.path(staging, "hm3.pco"), "u8", block_rows
    )
  )
  key_blocks <- pcodec_native_block_template(ordered_position, key_block_rows)
  flag_blocks <- pcodec_native_block_template(ordered_position, block_rows)
  index <- list(
    format = PCODEC_VARIANT_PANEL_INDEX_FORMAT, version = 1L,
    position_encoding = "delta_u32_within_block", rows = nrow(normalized),
    key_block_rows = key_block_rows, flag_block_rows = block_rows,
    key_blocks = key_blocks, flag_blocks = flag_blocks,
    streams = streams
  )
  jsonlite::write_json(index, file.path(staging, "native.index.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 17)
  files <- list(position = "position.pco", substitution = "substitution.pco",
                hm3 = "hm3.pco", index = "native.index.json")
  integrity_files <- stats::setNames(lapply(unname(unlist(files)), function(relative) {
    path <- file.path(staging, relative)
    list(bytes = as.numeric(file.info(path)$size),
         sha256 = digest::digest(path, algo = "sha256", file = TRUE))
  }), unname(unlist(files)))
  chromosomes <- compressor_chromosome_lengths(build)
  offsets <- pcodec_native_offsets(build)
  manifest <- list(
    format = "CompreSSoR-variant-panel", format_version = PCODEC_VARIANT_PANEL_FORMAT,
    backend = "pcodec", panel_universe = "core", genome_build = build, build = build,
    rows = nrow(normalized), n_rows = nrow(normalized), core_rows = nrow(normalized),
    hm3_rows = sum(ordered_hm3), hm3_flag = list(name = "hm3", dtype = "uint8",
                                                  values = c(0L, 1L), lossless = TRUE),
    identity = list(
      encoding = "native_global_position_plus_full_ref_alt_code",
      key = "chromosome:position:REF:ALT", sorted = TRUE, exact = TRUE,
      position_storage = "within-block delta-coded uint32; block first positions are in the index",
      substitution_storage = "directed REF-to-ALT code in uint8",
      external_reference_required = FALSE, chromosome_lengths = as.list(chromosomes),
      chromosome_offsets = as.list(offsets),
      effect_allele_is_alt = TRUE, other_allele_is_ref = TRUE
    ),
    codec = list(name = "pcodec_native_variant_panel_key_hm3_u8",
                 library = "pcodec", pco_version = "1.0.3", abi = "standalone",
                 page_rows = PCODEC_NATIVE_PAGE_ROWS,
                 key_block_rows = key_block_rows, flag_block_rows = block_rows),
    files = files, logical_columns = c("variant_id", "hm3"),
    reference = list(id = "none", status = "not_a_reference_genome_or_harmonisation"),
    preparation = list(method = "frozen_core_universe_plus_hm3_membership"),
    source_provenance = source_provenance %||% list(),
    tolerances = list(identity = "exact", hm3 = "lossless_binary_flag"),
    integrity = list(
      algorithm = "sha256", files = integrity_files,
      payload_sha256 = pcodec_payload_sha256(integrity_files),
      payload_definition = "sha256 over sorted relative payload paths, byte counts, and file sha256 values"
    ),
    created_utc = now_utc()
  )
  write_manifest(manifest, file.path(staging, "manifest.json"))
  seal_pcodec_manifest(file.path(staging, "manifest.json"))
  if (dir.exists(output)) unlink(output, recursive = TRUE, force = TRUE)
  if (!file.rename(staging, output)) stop("could not install variant panel: ", output, call. = FALSE)
  on.exit(NULL, add = TRUE)
  open_variant_panel(output, verify_payload = TRUE)
}

variant_panel_read_blocks <- function(store, index, stream, block_ids) {
  if (!length(block_ids)) return(integer())
  unlist(lapply(block_ids, function(block) {
    pcodec_native_read_stream_block(store, index, stream, block)
  }), use.names = FALSE)
}

variant_panel_stream_bytes <- function(index, stream, block_ids = NULL) {
  blocks <- index$streams[[stream]]$blocks
  if (is.null(block_ids)) block_ids <- seq_along(blocks)
  if (!length(block_ids)) return(0)
  sum(vapply(blocks[block_ids], function(block) as.numeric(block$length), numeric(1)))
}

read_variant_panel <- function(panel, chromosomes = NULL, rows = NULL, keys = NULL,
                               hm3_only = FALSE, verify_payload = FALSE,
                               identity_only = FALSE) {
  store <- if (inherits(panel, "compressor_variant_panel")) panel else {
    open_variant_panel(panel, verify_payload = verify_payload)
  }
  manifest <- store$manifest
  build <- compressor_normalize_build(manifest$genome_build %||% manifest$build)
  index <- variant_panel_index(store)
  n <- as.integer(manifest$n_rows %||% manifest$rows)
  key_blocks <- index$key_blocks
  flag_blocks <- index$flag_blocks
  key_block_ids <- seq_along(key_blocks)
  flag_block_ids <- seq_along(flag_blocks)
  if (!is.null(rows)) {
    rows <- as.numeric(rows)
    if (anyNA(rows) || any(rows < 1 | rows > n | rows != floor(rows))) {
      stop("rows must be valid one-based variant-panel row numbers", call. = FALSE)
    }
    row_zero <- as.integer(rows) - 1L
    key_stops <- vapply(key_blocks, function(block) as.numeric(block$row_stop), numeric(1))
    key_block_ids <- unique(findInterval(row_zero, key_stops) + 1L)
  } else if (!is.null(keys)) {
    targets <- pcodec_native_target_keys(unique(normalise_variant_key(keys)), build = build)
    selected <- vapply(key_blocks, function(block) {
      first <- as.numeric(block$first_position)
      last <- as.numeric(block$last_position)
      any(targets$position >= first & targets$position <= last)
    }, logical(1))
    key_block_ids <- which(selected)
  } else if (!is.null(chromosomes)) {
    chromosomes <- unique(normalise_chromosome(chromosomes))
    lengths <- compressor_chromosome_lengths(build)
    offsets <- pcodec_native_offsets(build)
    selected <- vapply(key_blocks, function(block) {
      first <- as.numeric(block$first_position)
      last <- as.numeric(block$last_position)
      any(vapply(chromosomes, function(chromosome) {
        if (!chromosome %in% names(lengths)) return(FALSE)
        lower <- offsets[[chromosome]]
        upper <- lower + lengths[[chromosome]] - 1
        first <= upper && last >= lower
      }, logical(1)))
    }, logical(1))
    key_block_ids <- which(selected)
  }
  if (!length(key_block_ids)) {
    position <- numeric()
    substitution <- integer()
    hm3 <- integer()
    source_rows <- integer()
    flag_block_ids <- integer()
  } else {
    position <- variant_panel_read_blocks(store, index, "position", key_block_ids)
    substitution <- as.integer(variant_panel_read_blocks(store, index, "substitution", key_block_ids))
    source_rows <- unlist(lapply(key_block_ids, function(block) {
      seq.int(as.integer(key_blocks[[block]]$row_start),
              as.integer(key_blocks[[block]]$row_stop) - 1L)
    }), use.names = FALSE)
    flag_stops <- vapply(flag_blocks, function(block) as.numeric(block$row_stop), numeric(1))
    flag_block_ids <- unique(findInterval(source_rows, flag_stops) + 1L)
    flag_values <- variant_panel_read_blocks(store, index, "hm3", flag_block_ids)
    flag_rows <- unlist(lapply(flag_block_ids, function(block) {
      seq.int(as.integer(flag_blocks[[block]]$row_start),
              as.integer(flag_blocks[[block]]$row_stop) - 1L)
    }), use.names = FALSE)
    hm3 <- as.integer(flag_values[match(source_rows, flag_rows)])
  }
  key_columns <- pcodec_native_key_columns(position, substitution, build = build)
  variant_id <- if (isTRUE(identity_only)) {
    rep(NA_character_, length(position))
  } else {
    compressor_variant_key(key_columns$chromosome,
                           key_columns$base_pair_location,
                           key_columns$reference_allele,
                           key_columns$alternate_allele, build = build)
  }
  keep <- rep(TRUE, length(position))
  if (!is.null(chromosomes)) {
    chromosomes <- unique(normalise_chromosome(chromosomes))
    keep <- keep & key_columns$chromosome %in% chromosomes
  }
  if (!is.null(rows)) {
    keep <- keep & source_rows %in% (as.integer(rows) - 1L)
  }
  if (!is.null(keys) && !isTRUE(identity_only)) {
    keys <- unique(normalise_variant_key(keys))
    keep <- keep & variant_id %in% keys
  }
  if (isTRUE(hm3_only)) keep <- keep & hm3 == 1L
  out <- if (isTRUE(identity_only)) {
    data.frame(global_position = as.numeric(position[keep]),
               substitution = as.integer(substitution[keep]), hm3 = hm3[keep],
               stringsAsFactors = FALSE)
  } else {
    data.frame(variant_id = variant_id[keep], hm3 = hm3[keep],
               stringsAsFactors = FALSE)
  }
  payload_bytes_read <- variant_panel_stream_bytes(index, "position", key_block_ids) +
    variant_panel_stream_bytes(index, "substitution", key_block_ids) +
    variant_panel_stream_bytes(index, "hm3", flag_block_ids)
  payload_bytes_total <- variant_panel_stream_bytes(index, "position") +
    variant_panel_stream_bytes(index, "substitution") +
    variant_panel_stream_bytes(index, "hm3")
  attr(out, "variant_panel_read") <- list(
    selective = !is.null(chromosomes) || !is.null(rows) || !is.null(keys),
    key_blocks = as.integer(key_block_ids), flag_blocks = as.integer(flag_block_ids),
    key_blocks_total = as.integer(length(key_blocks)),
    flag_blocks_total = as.integer(length(flag_blocks)),
    payload_bytes_read = as.numeric(payload_bytes_read),
    payload_bytes_total = as.numeric(payload_bytes_total)
  )
  attr(out, "variant_set_normalized_ids") <- TRUE
  attr(out, "variant_set_canonical") <- TRUE
  attr(out, "variant_set_identity") <- list(
    global_position = as.numeric(position[keep]),
    substitution = as.integer(substitution[keep]),
    code = compressor_identity_code(position[keep], substitution[keep])
  )
  attr(out, "variant_set_metadata") <- list(
    id = "bundled_core_hm3", name = if (isTRUE(hm3_only)) "hm3" else "core",
    source = "bundled", build = build, version = manifest$format_version,
    identity = "canonical_allele_aware", rows = nrow(out),
    core_rows = as.integer(manifest$core_rows), hm3_rows = as.integer(manifest$hm3_rows),
    hm3_only = isTRUE(hm3_only), manifest = manifest
  )
  out
}
