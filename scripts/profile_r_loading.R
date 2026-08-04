#!/usr/bin/env Rscript

# Profile the production R/Python bridge under an OS-warm cache. This isolates
# implementation costs; it is not the cold-cache release benchmark.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

store_path <- Sys.getenv(
  "COMPRESSOR_PROFILE_STORE",
  "/Volumes/crucial_x9/CompreSSoR-benchmarks/finngen-full-pcodec-v03-release.cpr"
)
output_path <- Sys.getenv("COMPRESSOR_PROFILE_OUTPUT", "")
runs <- as.integer(Sys.getenv("COMPRESSOR_PROFILE_RUNS", "5"))
stopifnot(dir.exists(store_path), runs >= 5L)

options(
  CompreSSoR.tempdir = Sys.getenv("COMPRESSOR_PROFILE_TMPDIR", tempdir()),
  CompreSSoR.native_bridge = TRUE,
  CompreSSoR.persistent_worker = TRUE
)
python <- Sys.getenv("COMPRESSOR_PYTHON", unset = "")
stopifnot(nzchar(python), file.exists(python))

store <- open_compressor(store_path)
namespace <- asNamespace("CompreSSoR")
worker_request <- get("pcodec_worker_request", namespace)
bridge_read <- get("read_pcodec_binary_bridge", namespace)
stop_worker <- get("pcodec_stop_worker", namespace)

column_sets <- list(
  analysis = c("z", "beta", "standard_error", "effect_allele_frequency"),
  public_full = c(
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "z", "beta", "standard_error", "effect_allele_frequency"
  )
)

touch <- function(data) {
  numeric_columns <- vapply(data, is.numeric, logical(1L))
  character_columns <- vapply(data, is.character, logical(1L))
  numeric_checksum <- sum(vapply(
    data[numeric_columns], function(column) sum(column, na.rm = TRUE), numeric(1L)
  ))
  character_checksum <- sum(vapply(
    data[character_columns], function(column) sum(nchar(column)), numeric(1L)
  ))
  numeric_checksum + character_checksum
}

records <- list()
record_index <- 0L
for (projection in names(column_sets)) {
  columns <- column_sets[[projection]]
  stop_worker()
  for (run in seq_len(runs)) {
    bridge <- tempfile(
      paste0("compressor-r-profile-", projection, "-"),
      tmpdir = getOption("CompreSSoR.tempdir"), fileext = ".bridge"
    )
    request_seconds <- system.time({
      worker_request(list(
        command = "read",
        store = normalizePath(store$path, mustWork = TRUE),
        output = bridge,
        rows = NULL,
        keys = NULL,
        chromosome = NULL,
        start = NULL,
        end = NULL,
        columns = columns
      ))
    })[["elapsed"]]
    bridge_bytes <- sum(file.info(list.files(
      bridge, recursive = TRUE, full.names = TRUE, all.files = TRUE,
      no.. = TRUE
    ))$size)
    bridge_metadata <- fromJSON(
      file.path(bridge, "bridge.json"), simplifyVector = FALSE
    )
    decoded <- NULL
    materialise_seconds <- system.time({
      decoded <- bridge_read(bridge)
    })[["elapsed"]]
    checksum <- NA_real_
    touch_seconds <- system.time({
      checksum <- touch(decoded)
    })[["elapsed"]]
    record_index <- record_index + 1L
    records[[record_index]] <- data.table(
      path = "split_bridge",
      projection = projection,
      run = run,
      python_decode_bridge_seconds = as.numeric(request_seconds),
      r_materialise_seconds = as.numeric(materialise_seconds),
      touch_seconds = as.numeric(touch_seconds),
      total_load_seconds = as.numeric(request_seconds + materialise_seconds),
      load_and_touch_seconds = as.numeric(
        request_seconds + materialise_seconds + touch_seconds
      ),
      rows = nrow(decoded),
      columns = ncol(decoded),
      bridge_bytes = bridge_bytes,
      source_bytes_read = as.numeric(
        bridge_metadata$source_bytes_read %||% NA_real_
      ),
      object_bytes = as.numeric(object.size(decoded)),
      checksum = checksum
    )
    unlink(bridge, recursive = TRUE, force = TRUE)
    rm(decoded)
    gc()
  }

  stop_worker()
  for (run in seq_len(runs)) {
    decoded <- NULL
    load_seconds <- system.time({
      decoded <- read_sumstats(store, columns = columns)
    })[["elapsed"]]
    checksum <- NA_real_
    touch_seconds <- system.time({
      checksum <- touch(decoded)
    })[["elapsed"]]
    record_index <- record_index + 1L
    records[[record_index]] <- data.table(
      path = "public_api",
      projection = projection,
      run = run,
      python_decode_bridge_seconds = NA_real_,
      r_materialise_seconds = NA_real_,
      touch_seconds = as.numeric(touch_seconds),
      total_load_seconds = as.numeric(load_seconds),
      load_and_touch_seconds = as.numeric(load_seconds + touch_seconds),
      rows = nrow(decoded),
      columns = ncol(decoded),
      bridge_bytes = NA_real_,
      source_bytes_read = NA_real_,
      object_bytes = as.numeric(object.size(decoded)),
      checksum = checksum
    )
    rm(decoded)
    gc()
  }
}
stop_worker()

records <- rbindlist(records, fill = TRUE)
checks <- records[, .(checksums = uniqueN(checksum)), by = .(projection)]
stopifnot(all(checks$checksums == 1L))
summary <- records[, .(
  runs = .N,
  median_load_seconds = median(total_load_seconds),
  min_load_seconds = min(total_load_seconds),
  max_load_seconds = max(total_load_seconds),
  median_load_and_touch_seconds = median(load_and_touch_seconds),
  median_python_decode_bridge_seconds = median(
    python_decode_bridge_seconds, na.rm = TRUE
  ),
  median_r_materialise_seconds = median(r_materialise_seconds, na.rm = TRUE),
  object_bytes = max(object_bytes),
  bridge_bytes = if (all(is.na(bridge_bytes))) NA_real_ else max(
    bridge_bytes, na.rm = TRUE
  ),
  source_bytes_read = if (all(is.na(source_bytes_read))) NA_real_ else median(
    source_bytes_read, na.rm = TRUE
  )
), by = .(path, projection)]

result <- list(
  schema_version = "1.0.0",
  benchmark_scope = "OS-warm component profile; not a cold-cache result",
  measured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_version = as.character(packageVersion("CompreSSoR")),
  store = normalizePath(store_path),
  rows = store$manifest$n_rows,
  repetitions = runs,
  summary = summary,
  runs = records
)
if (nzchar(output_path)) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  write_json(result, output_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
print(summary)
