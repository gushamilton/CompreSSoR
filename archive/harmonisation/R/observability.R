# Bounded, opt-in observability for the archived reference-backed workflow.
#
# This layer deliberately records phase metadata rather than data values.  The
# elapsed clock is proc.time()[["elapsed"]], which is monotonic for the R
# process and does not require a wall-clock timestamp for duration arithmetic.

compressor_harmonisation_observability_schema <-
  "CompreSSoR.harmonisation.observability.v1"

compressor_observability_clock <- function() {
  now <- proc.time()
  cpu <- suppressWarnings(sum(as.numeric(now[c("user.self", "sys.self")]),
                              na.rm = TRUE))
  wall <- suppressWarnings(as.numeric(now[["elapsed"]]))
  list(
    wall = if (is.finite(wall)) wall else NA_real_,
    cpu = if (is.finite(cpu)) cpu else NA_real_
  )
}

compressor_observability_scalar <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value)) NA_real_ else value
}

compressor_observability_read_rss_mb <- function() {
  # /proc is a read-only, bounded observation and is available on the Linux
  # compute nodes where the slow path is normally run.  Return NA elsewhere
  # rather than invoking a platform-specific process or GC side effect.
  status <- "/proc/self/status"
  if (!file.exists(status)) return(NA_real_)
  lines <- tryCatch(readLines(status, warn = FALSE), error = function(e) character())
  hit <- grep("^VmRSS:[[:space:]]*[0-9]+[[:space:]]*kB$", lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(sub("^VmRSS:[[:space:]]*([0-9]+).*$", "\\1", hit[1L])))
  if (is.finite(kb)) kb / 1024 else NA_real_
}

compressor_observability_normalise <- function(value = NULL) {
  if (is.null(value)) value <- list(level = "summary")
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    value <- list(level = if (isTRUE(value)) "summary" else "off")
  }
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    value <- list(level = value)
  }
  if (is.function(value)) value <- list(level = "events", callback = value)
  if (!is.list(value) || length(value) == 0L) {
    stop("observability must be NULL, logical, a level, a function, or a named list",
         call. = FALSE)
  }
  level <- tolower(as.character(value$level %||% "summary"))
  level[level %in% c("quiet", "none", "disabled")] <- "off"
  level[level %in% c("progress", "event", "logging")] <- "events"
  level <- match.arg(level, c("off", "summary", "events", "debug"))
  callback <- value$callback %||% value$progress_callback %||% NULL
  if (!is.null(callback) && !is.function(callback)) {
    stop("observability$callback must be a function", call. = FALSE)
  }
  jsonl <- value$jsonl %||% value$log %||% value$progress_log %||% NULL
  if (!is.null(jsonl) &&
      (!is.character(jsonl) || length(jsonl) != 1L || is.na(jsonl) || !nzchar(jsonl))) {
    stop("observability$jsonl must be one non-empty file path", call. = FALSE)
  }
  memory <- value$memory %||% value$sample_memory %||% FALSE
  if (length(memory) != 1L || !is.logical(memory) || is.na(memory)) {
    stop("observability$memory must be TRUE or FALSE", call. = FALSE)
  }
  if (identical(level, "off")) {
    callback <- NULL
    jsonl <- NULL
    memory <- FALSE
  }
  if (identical(level, "summary") && (!is.null(callback) || !is.null(jsonl))) {
    level <- "events"
  }
  list(level = level, callback = callback, jsonl = jsonl, memory = isTRUE(memory))
}

compressor_observability_create <- function(value = NULL, workers = 1L) {
  config <- compressor_observability_normalise(value)
  state <- new.env(parent = emptyenv())
  state$config <- config
  state$enabled <- !identical(config$level, "off")
  state$events_enabled <- state$enabled && config$level %in% c("events", "debug")
  state$memory_enabled <- state$enabled && isTRUE(config$memory)
  state$started <- compressor_observability_clock()
  state$ended <- NULL
  state$status <- "running"
  state$finalized <- FALSE
  state$next_phase_id <- 0L
  state$next_event_id <- 0L
  state$active <- list()
  state$phases <- list()
  state$event_count <- 0L
  state$logging_errors <- 0L
  state$peak_rss_mb <- NA_real_
  state$rows_in <- NA_real_
  state$rows_out <- NA_real_
  state$workers_requested <- max(1L, as.integer(workers)[1L])
  state$workers_effective <- state$workers_requested
  state$log_connection <- NULL
  if (state$events_enabled && !is.null(config$jsonl)) {
    dir.create(dirname(normalizePath(config$jsonl, mustWork = FALSE)),
               recursive = TRUE, showWarnings = FALSE)
    state$log_connection <- tryCatch(
      file(config$jsonl, open = "at", encoding = "UTF-8"),
      error = function(e) stop("could not open observability JSONL log: ",
                               conditionMessage(e), call. = FALSE)
    )
  }
  state
}

compressor_observability_child <- function(parent, workers = 1L) {
  if (is.null(parent)) return(NULL)
  child <- compressor_observability_create(parent$config, workers = workers)
  # Forked workers must never share a callback or file connection with the
  # parent. Their bounded phase records are returned and merged by the parent.
  child$events_enabled <- FALSE
  child$log_connection <- NULL
  child
}

compressor_observability_sample_memory <- function(state) {
  if (is.null(state) || !isTRUE(state$memory_enabled)) return(NA_real_)
  current <- compressor_observability_read_rss_mb()
  if (is.finite(current)) {
    state$peak_rss_mb <- max(state$peak_rss_mb %||% NA_real_, current, na.rm = TRUE)
    if (!is.finite(state$peak_rss_mb)) state$peak_rss_mb <- current
  }
  current
}

compressor_observability_emit <- function(state, event, phase = NULL,
                                          rows_in = NA_real_, rows_out = NA_real_,
                                          rows_dropped = NA_real_, workers = NULL,
                                          phase_elapsed = 0, phase_cpu = 0,
                                          extra = list()) {
  if (is.null(state) || !isTRUE(state$events_enabled)) return(invisible(NULL))
  clock <- compressor_observability_clock()
  elapsed <- clock$wall - state$started$wall
  cpu <- clock$cpu - state$started$cpu
  state$next_event_id <- state$next_event_id + 1L
  state$event_count <- state$event_count + 1L
  memory <- compressor_observability_sample_memory(state)
  event_data <- c(
    list(
      schema = compressor_harmonisation_observability_schema,
      event_id = state$next_event_id,
      event = as.character(event),
      phase = phase %||% NA_character_,
      elapsed_seconds = max(0, compressor_observability_scalar(elapsed)),
      cpu_seconds = max(0, compressor_observability_scalar(cpu)),
      phase_elapsed_seconds = max(0, compressor_observability_scalar(phase_elapsed)),
      phase_cpu_seconds = max(0, compressor_observability_scalar(phase_cpu)),
      rows_in = compressor_observability_scalar(rows_in),
      rows_out = compressor_observability_scalar(rows_out),
      rows_dropped = compressor_observability_scalar(rows_dropped),
      workers = as.integer(workers %||% state$workers_effective),
      rss_mb = memory
    ),
    extra
  )
  emit_one <- function(sink) {
    tryCatch(sink(event_data), error = function(e) {
      state$logging_errors <- state$logging_errors + 1L
      invisible(NULL)
    })
  }
  if (is.function(state$config$callback)) emit_one(state$config$callback)
  if (!is.null(state$log_connection)) {
    tryCatch({
      writeLines(jsonlite::toJSON(event_data, auto_unbox = TRUE, null = "null",
                                  digits = 16, pretty = FALSE),
                 state$log_connection, useBytes = TRUE)
      flush(state$log_connection)
    }, error = function(e) state$logging_errors <- state$logging_errors + 1L)
  }
  invisible(event_data)
}

compressor_observability_start_phase <- function(state, name, rows_in = NA_real_,
                                                 workers = NULL) {
  if (is.null(state) || !isTRUE(state$enabled)) return(NULL)
  name <- as.character(name)[1L]
  if (!nzchar(name)) stop("observability phase name must be non-empty", call. = FALSE)
  state$next_phase_id <- state$next_phase_id + 1L
  id <- as.character(state$next_phase_id)
  token <- list(id = id, name = name, started = compressor_observability_clock(),
                rows_in = rows_in, workers = workers %||% state$workers_effective)
  state$active[[id]] <- token
  compressor_observability_emit(state, "phase_start", name, rows_in = rows_in,
                                workers = token$workers)
  id
}

compressor_observability_finish_phase <- function(state, id, rows_out = NA_real_,
                                                  rows_dropped = NULL,
                                                  status = "completed", extra = list()) {
  if (is.null(state) || is.null(id) || !isTRUE(state$enabled)) return(invisible(NULL))
  token <- state$active[[as.character(id)]]
  if (is.null(token)) return(invisible(NULL))
  state$active[[as.character(id)]] <- NULL
  ended <- compressor_observability_clock()
  elapsed <- max(0, ended$wall - token$started$wall)
  cpu <- max(0, ended$cpu - token$started$cpu)
  input <- compressor_observability_scalar(token$rows_in)
  output <- compressor_observability_scalar(rows_out)
  dropped <- rows_dropped
  if (is.null(dropped) && is.finite(input) && is.finite(output)) {
    dropped <- max(0, input - output)
  }
  dropped <- compressor_observability_scalar(dropped)
  memory <- compressor_observability_sample_memory(state)
  previous <- state$phases[[token$name]]
  if (is.null(previous)) previous <- list(
    calls = 0L, elapsed_seconds = 0, cpu_seconds = 0,
    rows_in = 0, rows_out = 0, rows_dropped = 0,
    workers = token$workers, status = "completed"
  )
  previous$calls <- as.integer(previous$calls %||% 0L) + 1L
  previous$elapsed_seconds <- as.numeric(previous$elapsed_seconds %||% 0) + elapsed
  previous$cpu_seconds <- as.numeric(previous$cpu_seconds %||% 0) + cpu
  if (is.finite(input)) previous$rows_in <- as.numeric(previous$rows_in %||% 0) + input
  if (is.finite(output)) previous$rows_out <- as.numeric(previous$rows_out %||% 0) + output
  if (is.finite(dropped)) previous$rows_dropped <- as.numeric(previous$rows_dropped %||% 0) + dropped
  previous$workers <- max(as.integer(previous$workers %||% 1L), as.integer(token$workers))
  previous$status <- as.character(status)
  if (length(extra)) previous$metadata <- utils::modifyList(previous$metadata %||% list(), extra)
  if (is.finite(memory)) previous$peak_rss_mb <- max(previous$peak_rss_mb %||% memory,
                                                       memory, na.rm = TRUE)
  state$phases[[token$name]] <- previous
  compressor_observability_emit(
    state, "phase_end", token$name, rows_in = input, rows_out = output,
    rows_dropped = dropped, workers = token$workers,
    phase_elapsed = elapsed, phase_cpu = cpu,
    extra = c(list(status = status), extra)
  )
  invisible(previous)
}

compressor_observability_set_workers <- function(state, requested, effective) {
  if (is.null(state)) return(invisible(NULL))
  state$workers_requested <- max(1L, as.integer(requested)[1L])
  state$workers_effective <- max(1L, as.integer(effective)[1L])
  invisible(NULL)
}

compressor_observability_set_rows <- function(state, input = NULL, output = NULL) {
  if (is.null(state)) return(invisible(NULL))
  if (!is.null(input)) state$rows_in <- compressor_observability_scalar(input)
  if (!is.null(output)) state$rows_out <- compressor_observability_scalar(output)
  invisible(NULL)
}

compressor_observability_merge_phases <- function(state, phases) {
  if (is.null(state) || !isTRUE(state$enabled) || !is.list(phases)) return(invisible(NULL))
  for (name in names(phases)) {
    part <- phases[[name]]
    if (!is.list(part)) next
    current <- state$phases[[name]]
    if (is.null(current)) {
      state$phases[[name]] <- part
      next
    }
    for (field in c("calls", "elapsed_seconds", "cpu_seconds", "rows_in",
                    "rows_out", "rows_dropped")) {
      left <- compressor_observability_scalar(current[[field]] %||% 0)
      right <- compressor_observability_scalar(part[[field]] %||% 0)
      current[[field]] <- if (is.finite(left) && is.finite(right)) left + right else left
    }
    current$workers <- max(as.integer(current$workers %||% 1L),
                           as.integer(part$workers %||% 1L))
    if (is.finite(part$peak_rss_mb %||% NA_real_)) {
      current$peak_rss_mb <- max(current$peak_rss_mb %||% part$peak_rss_mb,
                                 part$peak_rss_mb, na.rm = TRUE)
    }
    current$status <- part$status %||% current$status
    current$metadata <- utils::modifyList(current$metadata %||% list(),
                                           part$metadata %||% list())
    state$phases[[name]] <- current
  }
  invisible(NULL)
}

compressor_observability_snapshot <- function(state) {
  if (is.null(state) || !isTRUE(state$enabled)) return(NULL)
  ended <- state$ended %||% compressor_observability_clock()
  elapsed <- max(0, ended$wall - state$started$wall)
  cpu <- max(0, ended$cpu - state$started$cpu)
  phases <- state$phases
  result <- list(
    schema = compressor_harmonisation_observability_schema,
    clock = "proc.time.elapsed",
    unit = "seconds",
    status = state$status,
    elapsed_seconds = elapsed,
    cpu_seconds = cpu,
    rows = list(
      input = if (is.finite(state$rows_in)) as.integer(state$rows_in) else NA_integer_,
      output = if (is.finite(state$rows_out)) as.integer(state$rows_out) else NA_integer_,
      dropped = if (is.finite(state$rows_in) && is.finite(state$rows_out))
        max(0L, as.integer(state$rows_in - state$rows_out)) else NA_integer_
    ),
    workers = list(requested = as.integer(state$workers_requested),
                   effective = as.integer(state$workers_effective)),
    phases = phases
  )
  if (isTRUE(state$memory_enabled)) {
    result$memory <- list(peak_rss_mb = if (is.finite(state$peak_rss_mb))
      state$peak_rss_mb else NA_real_, sampling = "phase_boundaries")
  }
  result$events <- list(emitted = as.integer(state$event_count),
                        logging_errors = as.integer(state$logging_errors))
  result
}

compressor_observability_finalize <- function(state, status = "completed", output_rows = NULL) {
  if (is.null(state) || isTRUE(state$finalized)) return(compressor_observability_snapshot(state))
  if (!is.null(output_rows)) compressor_observability_set_rows(state, output = output_rows)
  active_ids <- names(state$active)
  if (length(active_ids)) {
    for (id in active_ids) {
      compressor_observability_finish_phase(state, id,
                                            rows_out = state$rows_out,
                                            status = status)
    }
  }
  state$status <- as.character(status)[1L]
  state$ended <- compressor_observability_clock()
  compressor_observability_emit(state, "finalized", phase = NA_character_,
                                rows_in = state$rows_in, rows_out = state$rows_out,
                                rows_dropped = if (is.finite(state$rows_in) &&
                                                   is.finite(state$rows_out)) {
                                  max(0, state$rows_in - state$rows_out)
                                } else NA_real_,
                                extra = list(status = state$status))
  state$finalized <- TRUE
  compressor_observability_snapshot(state)
}

compressor_observability_close <- function(state) {
  if (is.null(state) || is.null(state$log_connection)) return(invisible(NULL))
  try(close(state$log_connection), silent = TRUE)
  state$log_connection <- NULL
  invisible(NULL)
}
