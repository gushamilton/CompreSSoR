## Optional Python/Pcodec backend.
##
## The durable format is written by inst/python/compressor_pcodec.py.  Keeping
## the subprocess boundary here means that the package remains installable on
## systems which use the Parquet backend, while Pcodec remains the primary
## format for users who configure its small Python environment.

.pcodec_state <- new.env(parent = emptyenv())

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

pcodec_python <- function() {
  candidates <- c(
    getOption("CompreSSoR.python", NULL),
    Sys.getenv("COMPRESSOR_PYTHON", unset = ""),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    resolved <- if (grepl("[/\\\\]", candidate)) candidate else Sys.which(candidate)
    if (nzchar(resolved) && file.exists(resolved)) return(resolved)
  }
  stop("Pcodec backend needs Python; set options(CompreSSoR.python = '/path/to/python') "
       , call. = FALSE)
}

pcodec_script <- function() {
  candidates <- c(
    system.file("python", "compressor_pcodec.py", package = "CompreSSoR"),
    file.path("inst", "python", "compressor_pcodec.py"),
    file.path(getwd(), "inst", "python", "compressor_pcodec.py")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("cannot find inst/python/compressor_pcodec.py", call. = FALSE)
  normalizePath(hit[1L], mustWork = TRUE)
}

pcodec_invocation <- function(refresh = FALSE) {
  python <- pcodec_python()
  cache_key <- paste(python, R.version$arch %||% "", sep = "|")
  if (!isTRUE(refresh) && identical(.pcodec_state$key %||% NULL, cache_key)) {
    return(.pcodec_state$invocation)
  }
  probes <- list(list(command = python, prefix = character()))
  if (identical(Sys.info()[["sysname"]], "Darwin") && file.exists("/usr/bin/arch")) {
    probes <- c(probes, list(
      list(command = "/usr/bin/arch", prefix = c("-arm64", python)),
      list(command = "/usr/bin/arch", prefix = c("-x86_64", python))
    ))
  }
  for (probe in probes) {
    status <- suppressWarnings(system2(
      probe$command,
      shQuote(c(
        probe$prefix, "-c",
        paste(
          "import sys; assert sys.version_info >= (3, 10);",
          "import numpy, pcodec, zstandard;",
          "from pcodec.wrapped import FileCompressor, FileDecompressor;",
          "assert hasattr(FileCompressor, 'write_header')"
        )
      )),
      stdout = FALSE, stderr = FALSE
    ))
    if (identical(as.integer(status), 0L)) {
      .pcodec_state$key <- cache_key
      .pcodec_state$invocation <- probe
      return(probe)
    }
  }
  stop("The selected Python cannot import numpy, pcodec, and zstandard; install them together in one environment or select another interpreter with COMPRESSOR_PYTHON",
       call. = FALSE)
}

pcodec_worker_enabled <- function() {
  isTRUE(getOption("CompreSSoR.persistent_worker", TRUE)) &&
    requireNamespace("processx", quietly = TRUE)
}

pcodec_stop_worker <- function() {
  worker <- .pcodec_state$worker %||% NULL
  if (!is.null(worker)) {
    if (isTRUE(tryCatch(worker$is_alive(), error = function(...) FALSE))) {
      try(worker$write_input('{"id":0,"command":"shutdown"}\n'), silent = TRUE)
      try(worker$poll_io(100L), silent = TRUE)
      if (isTRUE(tryCatch(worker$is_alive(), error = function(...) FALSE))) {
        try(worker$kill(), silent = TRUE)
      }
    }
    .pcodec_state$worker <- NULL
    .pcodec_state$worker_key <- NULL
    .pcodec_state$worker_checked_at <- NULL
  }
  invisible(NULL)
}

pcodec_worker_read_line <- function(worker, timeout = 30000L) {
  deadline <- unname(proc.time()[["elapsed"]]) + as.numeric(timeout) / 1000
  connection <- worker$get_output_connection()
  buffer <- ""
  # Read and frame the response on one low-level connection API. The explicit
  # accumulator handles replies split across arbitrary pipe chunks.
  spin_ms <- getOption("CompreSSoR.worker_spin_ms", 0)
  spin_ms <- max(0, min(10, as.numeric(spin_ms)))
  spin_deadline <- min(deadline, unname(proc.time()[["elapsed"]]) + spin_ms / 1000)
  repeat {
    chunk <- processx::conn_read_chars(connection, n = 65536L)
    if (length(chunk) && nzchar(chunk[[1L]])) {
      buffer <- paste0(buffer, paste0(chunk, collapse = ""))
      newline <- regexpr("\n", buffer, fixed = TRUE)[[1L]]
      if (newline > 0L) return(substr(buffer, 1L, newline - 1L))
    }
    if (!worker$is_alive()) {
      error <- paste(worker$read_error_lines(), collapse = "\n")
      if (!nzchar(error)) error <- "Pcodec worker stopped unexpectedly"
      stop(error, call. = FALSE)
    }
    now <- unname(proc.time()[["elapsed"]])
    remaining <- deadline - now
    if (remaining <= 0) stop("timed out waiting for the Pcodec worker", call. = FALSE)
    if (now < spin_deadline) next
    poll_ms <- getOption("CompreSSoR.worker_poll_ms", 1L)
    poll_ms <- max(1L, min(10L, as.integer(poll_ms)))
    processx::poll(
      list(connection),
      as.integer(min(poll_ms, ceiling(remaining * 1000)))
    )
  }
}

pcodec_worker <- function() {
  worker <- .pcodec_state$worker %||% NULL
  now <- unname(proc.time()[["elapsed"]])
  if (!is.null(worker) &&
      isTRUE(tryCatch(worker$is_alive(), error = function(...) FALSE)) &&
      now - (.pcodec_state$worker_checked_at %||% -Inf) < 1) {
    return(worker)
  }
  invocation <- pcodec_invocation()
  script <- pcodec_script()
  script_info <- file.info(script)
  worker_key <- paste(
    invocation$command, paste(invocation$prefix, collapse = "\037"), script,
    script_info$size, as.numeric(script_info$mtime), Sys.getpid(), sep = "|"
  )
  if (!is.null(worker) && identical(.pcodec_state$worker_key %||% NULL, worker_key) &&
      isTRUE(tryCatch(worker$is_alive(), error = function(...) FALSE))) {
    .pcodec_state$worker_checked_at <- now
    return(worker)
  }
  pcodec_stop_worker()
  worker <- processx::process$new(
    invocation$command,
    c(invocation$prefix, "-u", script, "serve"),
    stdin = "|", stdout = "|", stderr = "|",
    cleanup = TRUE, cleanup_tree = TRUE
  )
  ready <- tryCatch(
    jsonlite::fromJSON(pcodec_worker_read_line(worker), simplifyVector = FALSE),
    error = function(error) {
      try(worker$kill(), silent = TRUE)
      stop(conditionMessage(error), call. = FALSE)
    }
  )
  if (!identical(ready$format, "CompreSSoR-worker") ||
      !identical(as.integer(ready$version), 1L) || !isTRUE(ready$ready)) {
    try(worker$kill(), silent = TRUE)
    stop("Pcodec worker returned an invalid startup response", call. = FALSE)
  }
  .pcodec_state$worker <- worker
  .pcodec_state$worker_key <- worker_key
  .pcodec_state$worker_checked_at <- now
  .pcodec_state$worker_request <- 0L
  worker
}

pcodec_worker_request <- function(payload, timeout = 300000L) {
  worker <- pcodec_worker()
  request_id <- as.integer((.pcodec_state$worker_request %||% 0L) + 1L)
  .pcodec_state$worker_request <- request_id
  payload$id <- request_id
  request <- jsonlite::toJSON(
    payload, auto_unbox = TRUE, null = "null", digits = NA
  )
  tryCatch(
    worker$write_input(paste0(request, "\n")),
    error = function(error) {
      pcodec_stop_worker()
      stop("could not write to the Pcodec worker: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  response <- tryCatch(
    jsonlite::fromJSON(
      pcodec_worker_read_line(worker, timeout = timeout), simplifyVector = FALSE
    ),
    error = function(error) {
      pcodec_stop_worker()
      stop(conditionMessage(error), call. = FALSE)
    }
  )
  if (!identical(as.integer(response$id), request_id)) {
    pcodec_stop_worker()
    stop("Pcodec worker response ID did not match its request", call. = FALSE)
  }
  if (!isTRUE(response$ok)) {
    stop(response$error %||% "Pcodec worker request failed", call. = FALSE)
  }
  response$result
}

pcodec_open_store_cached <- function(store) {
  if (inherits(store, "compressor_store")) return(store)
  path <- normalizePath(store, mustWork = FALSE)
  manifest_path <- file.path(path, "manifest.json")
  checksum_path <- file.path(path, "manifest.sha256")
  info <- file.info(c(manifest_path, checksum_path))
  if (anyNA(info$size)) return(open_compressor(path))
  fingerprint <- paste(info$size, as.numeric(info$mtime), collapse = "|")
  cache <- .pcodec_state$store_cache %||% list()
  entry <- cache[[path]] %||% NULL
  if (!is.null(entry) && identical(entry$fingerprint, fingerprint)) {
    return(entry$store)
  }
  opened <- open_compressor(path)
  cache[[path]] <- list(fingerprint = fingerprint, store = opened)
  if (length(cache) > 32L) cache <- cache[seq.int(length(cache) - 31L, length(cache))]
  .pcodec_state$store_cache <- cache
  opened
}

reg.finalizer(.pcodec_state, function(environment) {
  try(pcodec_stop_worker(), silent = TRUE)
}, onexit = TRUE)

pcodec_tempdir <- function(preferred = NULL) {
  configured <- getOption("CompreSSoR.tempdir", NULL) %||%
    Sys.getenv("COMPRESSOR_TMPDIR", unset = "")
  candidates <- unique(c(configured, preferred, tempdir()))
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    candidate <- path.expand(candidate)
    if (dir.exists(candidate) && isTRUE(unname(file.access(candidate, 2L)) == 0L)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("no writable temporary directory is available for the Pcodec bridge", call. = FALSE)
}

pcodec_run <- function(args) {
  stdout <- tempfile("compressor-pcodec-stdout-")
  stderr <- tempfile("compressor-pcodec-stderr-")
  on.exit(unlink(c(stdout, stderr), force = TRUE), add = TRUE)
  invocation <- pcodec_invocation()
  status <- system2(invocation$command,
                    shQuote(c(invocation$prefix, pcodec_script(), args)),
                    stdout = stdout, stderr = stderr)
  if (!identical(as.integer(status), 0L)) {
    message <- paste(c(readLines(stderr, warn = FALSE), readLines(stdout, warn = FALSE)), collapse = "\n")
    if (!nzchar(message)) message <- "Python Pcodec command failed"
    stop(message, call. = FALSE)
  }
  invisible(readLines(stdout, warn = FALSE))
}

pcodec_write_store <- function(data, output, metadata = list()) {
  required <- c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                "beta", "standard_error", "effect_allele_frequency", "z")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Pcodec input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(data)) stop("cannot write an empty Pcodec store", call. = FALSE)

  bridge_tmp <- pcodec_tempdir(dirname(output))
  input_path <- tempfile("compressor-pcodec-input-", tmpdir = bridge_tmp, fileext = ".tsv")
  metadata_path <- tempfile("compressor-pcodec-metadata-", tmpdir = bridge_tmp, fileext = ".json")
  on.exit(unlink(c(input_path, metadata_path), force = TRUE), add = TRUE)
  payload <- data.frame(
    chromosome = as.character(data$chromosome),
    base_pair_location = as.integer(data$base_pair_location),
    effect_allele = toupper(as.character(data$effect_allele)),
    other_allele = toupper(as.character(data$other_allele)),
    beta = as.numeric(data$beta),
    standard_error = as.numeric(data$standard_error),
    effect_allele_frequency = as.numeric(data$effect_allele_frequency),
    z = as.numeric(data$z),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  data.table::fwrite(payload, input_path, sep = "\t", na = ".", quote = FALSE,
                     col.names = TRUE, showProgress = FALSE)
  jsonlite::write_json(metadata, metadata_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  pcodec_run(c("write", "--input", input_path, "--output", normalizePath(output, mustWork = FALSE),
               "--metadata", metadata_path))

  manifest_path <- file.path(output, "manifest.json")
  manifest <- read_manifest(manifest_path)
  manifest$format <- "CompreSSoR"
  manifest$backend <- "pcodec"
  manifest$n_rows <- as.integer(manifest$rows %||% nrow(data))
  manifest$block_rows <- as.integer(manifest$value_block_rows %||% 65536L)
  manifest$genome_build <- metadata$genome_build %||% "GRCh38"
  manifest$reference <- metadata$reference %||% list(
    id = "none", build = "GRCh38", status = "identity key is self-contained"
  )
  manifest$harmonisation <- metadata$harmonisation %||% list(method = "not recorded")
  manifest$logical_columns <- c("z", "standard_error", "effect_allele_frequency")
  manifest$derived_columns <- list(beta = "z * standard_error", p_value = "2 * pnorm(-abs(z))")
  manifest$source_columns <- metadata$source_columns %||% names(data)
  manifest$codec <- list(
    name = if (identical(manifest$format_version, "0.3.0-pcodec")) {
      "pcodec_wrapped_z9_eaf8_se6_paged"
    } else {
      "pcodec_z9_eaf8_se6_block_framed"
    },
    library = "pcodec",
    key = "recursive_escaped_position_gap_plus_full_ref_alt_code",
    compression = "Pcodec per independent numerical stream; block-framed Zstandard exceptions",
    z_bits = 9L, eaf_bits = 8L, se_bits = 6L,
    p_storage = "omitted; derived from z",
    beta_storage = "omitted; derived from z and standard_error"
  )
  manifest$variant_storage <- "self_contained_identity_key"
  manifest$variant_identity <- list(
    encoding = manifest$identity$encoding,
    rsid = "not_stored",
    external_reference_required = FALSE
  )
  manifest$benchmark <- metadata$benchmark %||% NULL
  manifest$benchmark_comparisons <- metadata$benchmark_comparisons %||% NULL
  manifest$source <- metadata$source %||% NULL
  manifest$created_utc <- now_utc()
  manifest$tolerances <- list(
    eaf_abs_max = 0.004,
    z_abs_max_central = 7 / (2 * (2^9 - 2)),
    z_central_range = c(-3.5, 3.5),
    se_profile = "block-centred log2 residual quantisation; see codec metadata",
    exception_precision = "float32",
    exception_layout = "one frame per value block",
    se_log2_half_step = 1 / 62,
    se_transform_multiplicative_half_step = 2^(1 / 62) - 1,
    exact_values_in_exception_sidecar = FALSE
  )
  write_manifest(manifest, manifest_path)
  seal_pcodec_manifest(manifest_path)
  invisible(manifest)
}

read_pcodec_binary_bridge <- function(path) {
  metadata <- jsonlite::fromJSON(file.path(path, "bridge.json"), simplifyVector = FALSE)
  if (!identical(metadata$format, "CompreSSoR-binary-bridge") ||
      !identical(as.integer(metadata$version), 1L)) {
    stop("unsupported Pcodec binary bridge", call. = FALSE)
  }
  n <- as.integer(metadata$rows)
  requested <- unlist(metadata$requested_columns, use.names = FALSE)
  identity_encoding <- metadata$identity$encoding %||% NULL
  compact_identity <- identical(
    identity_encoding, "global_position_substitution"
  )
  native_bridge <- (
    compact_identity || (
      identical(metadata$codec$encoding %||% NULL, "semantic_codes") &&
        isTRUE(getOption("CompreSSoR.native_bridge", TRUE))
    )
  ) &&
    is.loaded("compressor_read_pcodec_bridge", PACKAGE = "CompreSSoR")
  if (native_bridge) {
    codec <- metadata$codec %||% list(
      block_centers_log2_residual = numeric(), z_min = -3.5, z_max = 3.5,
      z_count = 510L, se_count = 62L, eaf_count = 255L,
      block_rows = 65536L
    )
    chromosome_lengths <- if (compact_identity) {
      as.numeric(unlist(
        metadata$identity$chromosome_lengths, use.names = FALSE
      ))
    } else {
      numeric()
    }
    return(.Call(
      "compressor_read_pcodec_bridge",
      normalizePath(path, mustWork = TRUE), as.character(requested), as.numeric(n),
      as.numeric(unlist(codec$block_centers_log2_residual, use.names = FALSE)),
      as.numeric(codec$z_min), as.numeric(codec$z_max),
      as.integer(codec$z_count), as.integer(codec$se_count), as.integer(codec$eaf_count),
      as.integer(codec$block_rows), chromosome_lengths, PACKAGE = "CompreSSoR"
    ))
  }
  read_column <- function(spec) {
    connection <- file(file.path(path, spec$file), open = "rb")
    on.exit(close(connection), add = TRUE)
    count <- as.integer(spec$length %||% n)
    switch(spec$dtype,
      int32 = readBin(connection, integer(), n = count, size = 4L, endian = "little"),
      uint8 = readBin(connection, integer(), n = count, size = 1L, signed = FALSE,
                      endian = "little"),
      uint16 = readBin(connection, integer(), n = count, size = 2L, signed = FALSE,
                       endian = "little"),
      float32 = readBin(connection, numeric(), n = count, size = 4L, endian = "little"),
      float64 = readBin(connection, numeric(), n = count, size = 8L, endian = "little"),
      stop("unsupported Pcodec bridge dtype: ", spec$dtype, call. = FALSE)
    )
  }
  values <- lapply(metadata$files, read_column)
  expected_lengths <- vapply(metadata$files, function(spec) as.integer(spec$length %||% n), integer(1))
  if (any(vapply(values, length, integer(1)) != expected_lengths)) {
    stop("truncated Pcodec binary bridge", call. = FALSE)
  }
  if (identical(metadata$codec$encoding %||% NULL, "semantic_codes")) {
    codec <- metadata$codec
    placeholder <- function(name, value) {
      if (is.null(values[[name]])) rep.int(as.integer(value), n) else as.integer(values[[name]])
    }
    decoded <- .Call(
      "compressor_decode_native",
      placeholder("z_code", codec$z_count),
      placeholder("se_code", codec$se_count),
      placeholder("eaf_code", 0L),
      as.numeric(codec$z_min), as.numeric(codec$z_max),
      as.integer(codec$z_count), as.integer(codec$se_count), as.integer(codec$eaf_count),
      as.integer(codec$z_bits), as.integer(codec$se_bits), as.integer(codec$eaf_bits),
      as.integer(codec$block_rows),
      as.numeric(unlist(codec$block_centers_log2_residual, use.names = FALSE)),
      as.integer(values$exception_row), as.numeric(values$exception_z),
      as.numeric(values$exception_se), as.numeric(values$exception_eaf),
      as.integer(values$exception_flags), FALSE, FALSE,
      PACKAGE = "CompreSSoR"
    )
    values <- c(values[intersect(c("chromosome", "base_pair_location",
                                  "effect_allele", "other_allele"), names(values))],
                decoded)
  }
  if ("chromosome" %in% names(values)) {
    chromosome_levels <- c(as.character(1:22), "X", "Y")
    values$chromosome <- chromosome_levels[values$chromosome]
  }
  allele_levels <- c("A", "C", "G", "T")
  if ("effect_allele" %in% names(values)) {
    values$effect_allele <- allele_levels[values$effect_allele + 1L]
  }
  if ("other_allele" %in% names(values)) {
    values$other_allele <- allele_levels[values$other_allele + 1L]
  }
  if ("beta" %in% unlist(metadata$requested_columns, use.names = FALSE) &&
      is.null(values$beta)) {
    values$beta <- values$z * values$standard_error
  }
  if ("p_value" %in% unlist(metadata$requested_columns, use.names = FALSE) &&
      is.null(values$p_value)) {
    values$p_value <- 2 * stats::pnorm(-abs(values$z))
  }
  out <- as.data.frame(values[requested], stringsAsFactors = FALSE,
                       check.names = FALSE)
  row.names(out) <- NULL
  out
}

pcodec_read_store <- function(store, region = NULL, variants = NULL, columns = NULL) {
  store <- if (inherits(store, "compressor_store")) store else open_compressor(store)
  m <- store$manifest
  if (!is.null(columns) && !length(columns)) {
    stop("columns must contain at least one column name", call. = FALSE)
  }
  key_variants <- !is.null(variants) && is.character(variants)
  if (!is.null(variants) && !is.numeric(variants) && !key_variants) {
    stop("Pcodec variants must be zero-based row IDs or canonical chromosome:position:REF:ALT keys",
         call. = FALSE)
  }
  n <- as.integer(m$n_rows %||% m$rows)
  if (!is.null(variants) && !key_variants) {
    numeric_rows <- as.numeric(variants)
    if (any(!is.finite(numeric_rows) | numeric_rows != floor(numeric_rows))) {
      stop("variants must contain finite whole-number zero-based row IDs", call. = FALSE)
    }
    variants <- unique(as.integer(numeric_rows))
  }
  if (!is.null(variants) && !key_variants &&
      any(is.na(variants) | variants < 0L | variants >= n)) {
    stop("variants must be valid zero-based row IDs", call. = FALSE)
  }
  if (key_variants) {
    variants <- unique(trimws(variants))
    if (any(is.na(variants) | !nzchar(variants))) {
      stop("canonical variant keys must be non-empty strings", call. = FALSE)
    }
  }
  if (!is.null(variants) && !length(variants)) {
    out <- data.frame(row = integer(), chromosome = character(), base_pair_location = integer(),
                      effect_allele = character(), other_allele = character(), z = numeric(),
                      beta = numeric(), standard_error = numeric(),
                      effect_allele_frequency = numeric(), p_value = numeric(),
                      stringsAsFactors = FALSE)
  } else {
    bounds <- read_region_bounds(region)
    output_path <- tempfile(
      "compressor-pcodec-read-", tmpdir = pcodec_tempdir(),
      fileext = ".bridge"
    )
    on.exit(unlink(output_path, recursive = TRUE, force = TRUE), add = TRUE)
    if (pcodec_worker_enabled()) {
      pcodec_worker_request(list(
        command = "read",
        store = normalizePath(store$path, mustWork = TRUE),
        output = output_path,
        rows = if (!key_variants) unname(variants) else NULL,
        keys = if (key_variants) unname(variants) else NULL,
        chromosome = if (!is.null(bounds)) as.character(bounds$chromosome) else NULL,
        start = if (!is.null(bounds)) as.integer(bounds$start) else NULL,
        end = if (!is.null(bounds)) as.integer(bounds$end) else NULL,
        columns = if (!is.null(columns)) unname(unique(columns)) else NULL,
        compact_identity = isTRUE(getOption("CompreSSoR.native_bridge", TRUE))
      ))
    } else {
      args <- c("read", "--store", normalizePath(store$path, mustWork = TRUE),
                "--output", output_path, "--output-format", "binary")
      if (key_variants) {
        keys_path <- tempfile(
          "compressor-pcodec-keys-", tmpdir = pcodec_tempdir(dirname(store$path))
        )
        on.exit(unlink(keys_path, force = TRUE), add = TRUE)
        writeLines(variants, keys_path, useBytes = TRUE)
        args <- c(args, "--keys-file", keys_path)
      } else if (!is.null(variants)) {
        args <- c(args, "--rows", paste(variants, collapse = ","))
      }
      if (!is.null(bounds)) {
        args <- c(args, "--chromosome", as.character(bounds$chromosome),
                  "--start", as.character(as.integer(bounds$start)),
                  "--end", as.character(as.integer(bounds$end)))
      }
      if (!is.null(columns)) {
        args <- c(args, "--columns", paste(unique(columns), collapse = ","))
      }
      if (!isTRUE(getOption("CompreSSoR.native_bridge", TRUE))) {
        args <- c(args, "--expanded-identity-bridge")
      }
      pcodec_run(args)
    }
    out <- read_pcodec_binary_bridge(output_path)
  }
  if ("row" %in% names(out)) out <- out[order(out$row), , drop = FALSE]
  if ("row" %in% names(out)) out$row <- NULL
  if (is.null(columns)) return(out)
  missing <- setdiff(columns, names(out))
  if (length(missing)) stop("requested columns are not present: ", paste(missing, collapse = ", "), call. = FALSE)
  out[columns]
}

pcodec_read_stores <- function(stores, variants, columns, threads = 1L) {
  if (!length(stores) || length(stores) != length(variants)) {
    stop("stores and variants must have the same non-zero length", call. = FALSE)
  }
  if (!length(columns) || anyNA(columns) || any(!nzchar(columns))) {
    stop("columns must contain at least one non-empty column name", call. = FALSE)
  }
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) ||
      !is.finite(threads) || threads < 1 || threads != floor(threads)) {
    stop("threads must be one positive integer", call. = FALSE)
  }
  stores <- lapply(stores, function(store) {
    opened <- pcodec_open_store_cached(store)
    if (!identical(opened$manifest$backend, "pcodec")) {
      stop("batched reads require Pcodec CompreSSoR stores", call. = FALSE)
    }
    opened
  })
  variants <- lapply(variants, function(keys) {
    if (!is.character(keys) || anyNA(keys) || any(!nzchar(trimws(keys)))) {
      stop("each variants element must contain canonical variant keys", call. = FALSE)
    }
    unique(trimws(keys))
  })
  bridge_tmp <- pcodec_tempdir()
  output_path <- tempfile("compressor-pcodec-batch-", tmpdir = bridge_tmp,
                          fileext = ".bridge")
  on.exit(unlink(output_path, recursive = TRUE, force = TRUE), add = TRUE)
  reads <- lapply(seq_along(stores), function(index) {
    list(
      store = jsonlite::unbox(normalizePath(stores[[index]]$path, mustWork = TRUE)),
      keys = unname(variants[[index]]),
      columns = unname(unique(columns))
    )
  })
  coalesce <- isTRUE(getOption("CompreSSoR.coalesce_batch_reads", TRUE))
  # A batch gets one Python process and one binary bridge for all reads.  The
  # persistent processx worker is ideal for single requests, but large JSON
  # writes to its stdin can stall on macOS before Python receives the command.
  request_path <- tempfile("compressor-pcodec-batch-", tmpdir = bridge_tmp,
                           fileext = ".json")
  on.exit(unlink(request_path, force = TRUE), add = TRUE)
  jsonlite::write_json(
    list(reads = reads, coalesce = coalesce), request_path,
    auto_unbox = FALSE
  )
  pcodec_run(c("batch-read", "--request", request_path, "--output", output_path,
               "--threads", as.character(as.integer(threads))))
  response <- jsonlite::fromJSON(
    file.path(output_path, "batch.json"), simplifyVector = FALSE
  )
  if (!identical(response$format, "CompreSSoR-batch-bridge") ||
      length(response$reads) != length(stores)) {
    stop("Pcodec batch bridge response is malformed", call. = FALSE)
  }
  bridge_ids <- vapply(seq_along(stores), function(index) {
    bridge <- as.integer(response$reads[[index]]$bridge)
    if (length(bridge) != 1L || is.na(bridge) || bridge < 0L ||
        bridge >= length(stores)) {
      stop("Pcodec batch bridge index is malformed", call. = FALSE)
    }
    bridge
  }, integer(1))
  unique_bridges <- unique(bridge_ids)
  decoded <- lapply(unique_bridges, function(bridge) {
    read_pcodec_binary_bridge(file.path(output_path, bridge))
  })
  names(decoded) <- as.character(unique_bridges)
  lapply(bridge_ids, function(bridge) {
    out <- decoded[[as.character(bridge)]]
    if ("row" %in% names(out)) {
      out <- out[order(out$row), , drop = FALSE]
      out$row <- NULL
    }
    missing <- setdiff(columns, names(out))
    if (length(missing)) {
      stop("requested columns are not present: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
    out[columns]
  })
}

pcodec_validate_store <- function(store, full = FALSE) {
  tryCatch({
    s <- if (inherits(store, "compressor_store")) store else open_compressor(store)
    args <- c("validate", "--store", normalizePath(s$path, mustWork = TRUE))
    if (isTRUE(full)) args <- c(args, "--full")
    output <- pcodec_run(args)
    jsonlite::fromJSON(paste(output, collapse = "\n"), simplifyVector = FALSE)
  }, error = function(e) {
    list(valid = FALSE, errors = conditionMessage(e), rows = NA_integer_,
         profile = NA_character_, full = isTRUE(full))
  })
}
