variant_set_aliases <- c(
  common = "common", common_variants = "common",
  core = "core", core_variants = "core",
  tag = "tag", tags = "tag", tag_variants = "tag",
  hm3 = "hm3", hm3_variants = "hm3", hapmap3 = "hm3"
)

variant_set_environment_names <- function(name) {
  upper <- toupper(name)
  c(
    paste0("COMPRESSOR_", upper, "_VARIANTS"),
    paste0("COMPRESSOR_VARIANT_SET_", upper),
    paste0("COMPRESSOR_PANEL_", upper)
  )
}

variant_set_hash_environment_names <- function(name) {
  upper <- toupper(name)
  c(
    paste0("COMPRESSOR_", upper, "_VARIANTS_SHA256"),
    paste0("COMPRESSOR_VARIANT_SET_", upper, "_SHA256"),
    paste0("COMPRESSOR_PANEL_", upper, "_SHA256")
  )
}

variant_set_cache_dir <- function(cache_dir = NULL) {
  cache_dir <- cache_dir %||% Sys.getenv("COMPRESSOR_VARIANT_SET_CACHE", unset = "")
  if (!nzchar(cache_dir)) cache_dir <- Sys.getenv("COMPRESSOR_PANEL_CACHE", unset = "")
  if (!nzchar(cache_dir)) cache_dir <- tools::R_user_dir("CompreSSoR", which = "cache")
  normalizePath(path.expand(cache_dir), mustWork = FALSE)
}

variant_set_sha256 <- function(value) {
  value <- tolower(trimws(as.character(value %||% "")))
  if (length(value) != 1L || !grepl("^[0-9a-f]{64}$", value)) return(NA_character_)
  value
}

parse_named_variant_set_request <- function(variant_set) {
  if (!is.character(variant_set) || length(variant_set) != 1L || is.na(variant_set)) {
    return(NULL)
  }
  request <- trimws(variant_set)
  pieces <- regexec("^([^@#]+)(?:[@#](?:sha256[:=]?)?([0-9A-Fa-f]{64}))?$",
                    request, perl = TRUE)
  match <- regmatches(request, pieces)[[1L]]
  if (!length(match)) return(NULL)
  name <- tolower(trimws(match[[2L]]))
  canonical_name <- unname(variant_set_aliases[name])
  if (is.na(canonical_name)) return(NULL)
  pin <- if (length(match) >= 3L && nzchar(match[[3L]])) {
    variant_set_sha256(match[[3L]])
  } else {
    NA_character_
  }
  list(name = canonical_name, request = request, sha256 = pin)
}

variant_set_manifest <- function(path) {
  candidates <- if (dir.exists(path)) {
    c(file.path(path, "panel_manifest.json"),
      file.path(dirname(path), "panel_manifest.json"),
      file.path(path, "manifest.json"))
  } else {
    c(file.path(dirname(path), "panel_manifest.json"),
      paste0(path, ".manifest.json"),
      sub("[.]tsv(?:[.]gz)?$", ".manifest.json", path, ignore.case = TRUE))
  }
  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates) || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  tryCatch(jsonlite::read_json(candidates[[1L]], simplifyVector = FALSE),
           error = function(e) NULL)
}

variant_set_manifest_panel <- function(manifest, name = NULL) {
  if (is.null(manifest) || !is.list(manifest)) return(list())
  panels <- manifest$panels
  if (!is.list(panels)) return(list())
  candidates <- unique(c(name, names(panels)))
  hit <- candidates[!is.na(candidates) & nzchar(candidates) & candidates %in% names(panels)]
  if (length(hit)) panels[[hit[[1L]]]] else list()
}

variant_set_digest_target <- function(path, name = NULL, manifest = NULL) {
  if (!dir.exists(path)) return(path)
  panel <- variant_set_manifest_panel(manifest %||% variant_set_manifest(path), name)
  prepared <- panel$prepared %||% NULL
  if (is.character(prepared) && length(prepared) == 1L && file.exists(prepared)) {
    return(prepared)
  }
  if (is.character(prepared) && length(prepared) == 1L) {
    candidate <- file.path(path, basename(prepared))
    if (file.exists(candidate)) return(candidate)
    candidate <- file.path(dirname(path), basename(prepared))
    if (file.exists(candidate)) return(candidate)
  }
  manifest_path <- if (file.exists(file.path(path, "panel_manifest.json"))) {
    file.path(path, "panel_manifest.json")
  } else if (file.exists(file.path(dirname(path), "panel_manifest.json"))) {
    file.path(dirname(path), "panel_manifest.json")
  } else {
    file.path(path, "manifest.json")
  }
  if (file.exists(manifest_path)) manifest_path else path
}

variant_set_declared_sha256 <- function(path, name = NULL, manifest = NULL) {
  manifest <- manifest %||% variant_set_manifest(path)
  panel <- variant_set_manifest_panel(manifest, name)
  declared <- panel$prepared_sha256 %||% manifest$sha256 %||% NA_character_
  variant_set_sha256(declared)
}

verify_variant_set_hash <- function(path, expected, name = NULL, manifest = NULL) {
  expected <- variant_set_sha256(expected)
  if (is.na(expected)) return(list(sha256 = NA_character_, verified = FALSE, target = NULL))
  manifest <- manifest %||% variant_set_manifest(path)
  target <- variant_set_digest_target(path, name = name, manifest = manifest)
  if (!file.exists(target) || dir.exists(target)) {
    stop("hash-pinned variant set is missing its digest target: ", target, call. = FALSE)
  }
  observed <- tolower(reference_file_digest(target, "sha256"))
  if (!identical(observed, expected)) {
    stop("variant-set SHA-256 mismatch for ", name %||% "panel",
         ": expected ", expected, ", observed ", observed, call. = FALSE)
  }
  list(sha256 = observed, verified = TRUE,
       target = normalizePath(target, mustWork = FALSE))
}

download_cached_variant_set <- function(url, name, expected, cache_dir) {
  expected <- variant_set_sha256(expected)
  if (is.na(expected)) {
    stop("remote named variant sets require a SHA-256 pin", call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  extension <- sub("^[^?]*", "", basename(sub("[?].*$", "", url)))
  if (!nzchar(extension) || !grepl("[.]", extension)) extension <- ".tsv.gz"
  target <- file.path(cache_dir, paste0(name, "-", expected, extension))
  if (file.exists(target)) {
    verified <- verify_variant_set_hash(target, expected, name = name)
    return(list(path = target, cached = TRUE, verified = verified$verified,
                source_url = url, cache_dir = cache_dir))
  }
  temporary <- tempfile(paste0("compressor-panel-", name, "-"), tmpdir = cache_dir)
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  status <- tryCatch(utils::download.file(url, temporary, mode = "wb", quiet = TRUE,
                                          method = "libcurl"),
                     error = function(e) stop("could not download variant set: ",
                                              conditionMessage(e), call. = FALSE))
  if (!identical(status, 0L) || !file.exists(temporary)) {
    stop("could not download variant set: ", url, call. = FALSE)
  }
  verified <- verify_variant_set_hash(temporary, expected, name = name)
  if (!file.rename(temporary, target)) stop("could not install cached variant set: ", target,
                                             call. = FALSE)
  list(path = target, cached = FALSE, verified = verified$verified,
       source_url = url, cache_dir = cache_dir)
}

resolve_named_variant_set <- function(variant_set, cache_dir = NULL) {
  request <- parse_named_variant_set_request(variant_set)
  if (is.null(request)) return(NULL)
  name <- request$name
  env_names <- variant_set_environment_names(name)
  env_var <- env_names[[1L]]
  path <- ""
  for (candidate in env_names) {
    value <- Sys.getenv(candidate, unset = "")
    if (nzchar(value)) {
      path <- value
      env_var <- candidate
      break
    }
  }
  expected <- request$sha256
  if (is.na(expected)) {
    for (candidate in variant_set_hash_environment_names(name)) {
      value <- variant_set_sha256(Sys.getenv(candidate, unset = ""))
      if (!is.na(value)) {
        expected <- value
        break
      }
    }
  }
  if (!nzchar(path)) {
    root <- Sys.getenv("COMPRESSOR_VARIANT_SET_DIR", unset = "")
    if (nzchar(root)) {
      candidates <- switch(name,
        common = c("common.parquet", "common.cpr", "common.tsv.gz", "common.bim.gz"),
        core = c("core.parquet", "core.cpr", "core.tsv.gz", "core.bim.gz"),
        tag = c("tag.parquet", "tag.cpr", "tags.parquet", "tag.tsv.gz", "tag.bim.gz"),
        hm3 = c("hm3.parquet", "hm3.cpr", "hm3.tsv.gz", "hm3.bim.gz")
      )
      hits <- file.path(path.expand(root), candidates)
      hits <- hits[file.exists(hits)]
      if (length(hits)) path <- hits[[1L]]
    }
  }
  if (!nzchar(path)) {
    stop("variant_set='", name, "' needs ", env_var,
         " (or COMPRESSOR_VARIANT_SET_DIR with a named panel file)", call. = FALSE)
  }
  source_url <- NULL
  cached <- FALSE
  if (grepl("^https?://", path, ignore.case = TRUE)) {
    cached_panel <- download_cached_variant_set(path, name, expected, variant_set_cache_dir(cache_dir))
    source_url <- cached_panel$source_url
    cached <- cached_panel$cached
    path <- cached_panel$path
  } else {
    path <- path.expand(path)
    if (!file.exists(path) && !dir.exists(path)) {
      stop(env_var, " points to a missing variant-set file: ", path, call. = FALSE)
    }
  }
  manifest <- variant_set_manifest(path)
  declared <- variant_set_declared_sha256(path, name = name, manifest = manifest)
  if (is.na(expected)) expected <- declared
  if (name %in% c("core", "hm3") && is.na(expected)) {
    stop("named ", name, " panel is not hash-pinned; set ",
         variant_set_hash_environment_names(name)[[1L]],
         " or provide name@sha256:<digest>", call. = FALSE)
  }
  verification <- if (!is.na(expected)) {
    verify_variant_set_hash(path, expected, name = name, manifest = manifest)
  } else {
    list(sha256 = NA_character_, verified = FALSE, target = NULL)
  }
  list(name = name, requested = request$request, path = normalizePath(path, mustWork = FALSE),
       env_var = env_var, expected_sha256 = expected,
       sha256 = verification$sha256 %||% declared,
       hash_pinned = !is.na(expected), hash_verified = isTRUE(verification$verified),
       digest_target = verification$target, source_url = source_url,
       cache_dir = if (cached) variant_set_cache_dir(cache_dir) else NULL,
       cache_hit = if (cached) cached else NULL)
}

variant_set_provenance <- function(path = NULL, name = NULL, named = NULL,
                                   rows = NA_integer_, columns_read = NULL,
                                   identity = "unknown", manifest = NULL) {
  manifest <- manifest %||% if (!is.null(path)) variant_set_manifest(path) else NULL
  panel <- variant_set_manifest_panel(manifest, name %||% named$name %||% NULL)
  digest_target <- if (!is.null(path)) {
    variant_set_digest_target(path, name = name %||% named$name %||% NULL,
                              manifest = manifest)
  } else {
    NULL
  }
  observed <- if (!is.null(digest_target) && file.exists(digest_target) &&
                  !dir.exists(digest_target)) {
    tolower(reference_file_digest(digest_target, "sha256"))
  } else {
    NA_character_
  }
  expected <- named$expected_sha256 %||% NA_character_
  if (length(expected) != 1L || is.na(expected) || !nzchar(expected)) {
    expected <- panel$prepared_sha256 %||% NA_character_
  }
  expected <- variant_set_sha256(expected)
  if (!is.na(expected) && !is.na(observed) && !identical(expected, observed)) {
    stop("variant-set SHA-256 mismatch for ", name %||% named$name %||% "panel",
         ": expected ", expected, ", observed ", observed, call. = FALSE)
  }
  list(
    id = named$name %||% name %||% "local",
    name = named$name %||% name %||% NULL,
    requested = named$requested %||% NULL,
    source = if (!is.null(named$source_url)) "url" else if (!is.null(path)) "file" else "in_memory",
    local_path = if (!is.null(path)) normalizePath(path, mustWork = FALSE) else NULL,
    source_url = named$source_url %||% panel$source_url %||% NULL,
    cache_dir = named$cache_dir %||% NULL,
    cache_hit = named$cache_hit %||% NULL,
    build = manifest$build %||% panel$build %||% "GRCh38",
    version = manifest$format_version %||% panel$version %||% NULL,
    source_sha256 = panel$source_sha256 %||% NULL,
    sha256 = observed %||% NA_character_,
    expected_sha256 = expected,
    hash_pinned = !is.na(expected),
    hash_verified = isTRUE(named$hash_verified) ||
      (!is.na(expected) && !is.na(observed) && identical(expected, observed)),
    digest_target = if (!is.null(digest_target)) normalizePath(digest_target, mustWork = FALSE) else NULL,
    identity = identity,
    rows = as.integer(rows),
    columns_read = columns_read %||% NULL,
    manifest = if (!is.null(manifest)) manifest else NULL
  )
}
normalise_variant_set_columns <- function(data) {
  data <- clean_input_names(as.data.frame(data, stringsAsFactors = FALSE,
                                           check.names = FALSE))
  alias_map <- list(
    variant_id = c("variant_id", "variant_id_grch38", "vid", "SNPID", "SNP", "snp", "ID", "id"),
    rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID", "ID"),
    chromosome = c("chromosome", "chr", "CHR", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "BP", "bp"),
    other_allele = c("other_allele", "NEA", "REF", "ref", "A2"),
    effect_allele = c("effect_allele", "EA", "ALT", "alt", "A1")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  has_coordinate_identity <- all(c("chromosome", "base_pair_location",
                                   "other_allele", "effect_allele") %in% names(data))
  if (!"variant_id" %in% names(data) && !"rsid" %in% names(data) &&
      !has_coordinate_identity) {
    stop("variant set needs variant_id, rsid, or chromosome/position/allele columns", call. = FALSE)
  }
  if (!"variant_id" %in% names(data)) data$variant_id <- NA_character_
  if (!"rsid" %in% names(data)) data$rsid <- NA_character_
  if (!"chromosome" %in% names(data)) data$chromosome <- NA_character_
  if (!"base_pair_location" %in% names(data)) data$base_pair_location <- NA_integer_
  if (!"other_allele" %in% names(data)) data$other_allele <- NA_character_
  if (!"effect_allele" %in% names(data)) data$effect_allele <- NA_character_
  data$variant_id <- toupper(sub("^chr", "", trimws(as.character(data$variant_id)),
                                ignore.case = TRUE))
  data$rsid <- as.character(data$rsid)
  data$chromosome <- normalise_chromosome(data$chromosome)
  data$base_pair_location <- suppressWarnings(as.integer(data$base_pair_location))
  data$other_allele <- toupper(trimws(as.character(data$other_allele)))
  data$effect_allele <- toupper(trimws(as.character(data$effect_allele)))
  missing_text <- function(x) is.na(x) | !nzchar(x) | x %in% c(".", "NA")
  data$variant_id[missing_text(data$variant_id)] <- NA_character_
  data$rsid[missing_text(data$rsid)] <- NA_character_
  data$chromosome[missing_text(data$chromosome)] <- NA_character_
  data$other_allele[missing_text(data$other_allele)] <- NA_character_
  data$effect_allele[missing_text(data$effect_allele)] <- NA_character_

  canonical_pattern <- "^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$"
  supplied_canonical <- grepl(canonical_pattern, data$variant_id)
  parsed_chromosome <- sub(":.*$", "", data$variant_id)
  parsed_position <- suppressWarnings(as.integer(sub("^[^:]+:([^:]+):.*$", "\\1", data$variant_id)))
  parsed_other <- sub("^[^:]+:[^:]+:([^:]+):.*$", "\\1", data$variant_id)
  parsed_effect <- sub("^[^:]+:[^:]+:[^:]+:([^:]+)$", "\\1", data$variant_id)
  fill_chromosome <- supplied_canonical & is.na(data$chromosome)
  fill_position <- supplied_canonical & is.na(data$base_pair_location)
  fill_other <- supplied_canonical & is.na(data$other_allele)
  fill_effect <- supplied_canonical & is.na(data$effect_allele)
  data$chromosome[fill_chromosome] <- parsed_chromosome[fill_chromosome]
  data$base_pair_location[fill_position] <- parsed_position[fill_position]
  data$other_allele[fill_other] <- parsed_other[fill_other]
  data$effect_allele[fill_effect] <- parsed_effect[fill_effect]

  valid_columns <- data$chromosome %in% c(as.character(1:22), "X", "Y") &
    is.finite(data$base_pair_location) & data$base_pair_location >= 1L &
    !is.na(data$other_allele) & data$other_allele %in% c("A", "C", "G", "T") &
    !is.na(data$effect_allele) & data$effect_allele %in% c("A", "C", "G", "T") &
    data$other_allele != data$effect_allele
  canonical_from_columns <- rep(NA_character_, nrow(data))
  if (any(valid_columns)) {
    canonical_from_columns[valid_columns] <- compressor_variant_key(
      data$chromosome[valid_columns], data$base_pair_location[valid_columns],
      data$other_allele[valid_columns], data$effect_allele[valid_columns]
    )
  }
  # Coordinates plus REF/ALT are authoritative when present. This turns a
  # panel carrying rsIDs or a project-specific ID into the same directed,
  # allele-aware identity used by the store.
  use_canonical <- !is.na(canonical_from_columns)
  contradictory <- supplied_canonical & use_canonical &
    data$variant_id != canonical_from_columns
  if (any(contradictory)) {
    stop("variant set has conflicting canonical ID and chromosome/allele identity at row(s): ",
         paste(utils::head(which(contradictory), 5L), collapse = ", "), call. = FALSE)
  }
  data$variant_id[use_canonical] <- canonical_from_columns[use_canonical]
  canonical_values <- !is.na(data$variant_id) & grepl(canonical_pattern, data$variant_id)
  attr(data, "variant_set_canonical") <- all(canonical_values | is.na(data$variant_id)) &&
    any(canonical_values)
  attr(data, "variant_set_normalized_ids") <- isTRUE(attr(data, "variant_set_canonical"))
  attr(data, "variant_set_canonical_keys") <- unique(data$variant_id[canonical_values])
  data
}

read_variant_set <- function(variant_set, chromosomes = NULL) {
  named <- resolve_named_variant_set(variant_set)
  if (!is.data.frame(variant_set) && !is.null(named)) variant_set <- named$path
  if (is.data.frame(variant_set)) {
    out <- normalise_variant_set_columns(variant_set)
    metadata <- list(
      id = "in_memory", name = NULL, source = "in_memory",
      sha256 = digest::digest(out, algo = "sha256", serialize = TRUE),
      expected_sha256 = NA_character_, hash_pinned = FALSE,
      hash_verified = FALSE, identity = if (isTRUE(attr(out, "variant_set_canonical"))) {
        "canonical_allele_aware"
      } else {
        "legacy"
      }, rows = nrow(out)
    )
  } else {
    if (length(variant_set) != 1L || !is.character(variant_set) ||
        (!file.exists(variant_set) && !dir.exists(variant_set))) {
      stop("variant_set must be a data.frame, panel file, or chromosome-shard directory",
           call. = FALSE)
    }
    if (dir.exists(variant_set) && file.exists(file.path(variant_set, "manifest.json"))) {
      store <- open_compressor(variant_set)
      out <- read_sumstats(store, columns = c("chromosome", "base_pair_location",
                                               "effect_allele", "other_allele"))
      out$variant_id <- compressor_variant_key(
        out$chromosome, out$base_pair_location,
        out$other_allele, out$effect_allele
      )
      out$rsid <- NA_character_
      metadata <- variant_set_provenance(
        path = variant_set, name = named$name %||% NULL, named = named,
        rows = nrow(out),
        columns_read = c("chromosome", "base_pair_location", "effect_allele", "other_allele"),
        identity = "canonical_allele_aware"
      )
      attr(out, "variant_set_normalized_ids") <- TRUE
      attr(out, "variant_set_canonical") <- TRUE
      attr(out, "variant_set_metadata") <- metadata
      return(out)
    } else if (dir.exists(variant_set)) {
      # A staged panel directory contains one canonical TSV(.gz) shard per
      # chromosome. Read only the shards represented in the already
      # harmonised input; this is especially valuable for chr-local jobs.
      shard_dir <- if (dir.exists(file.path(variant_set, "by_chrom"))) {
        file.path(variant_set, "by_chrom")
      } else {
        variant_set
      }
      shard_paths <- list.files(
        shard_dir,
        pattern = "^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)[.]tsv(?:[.]gz)?$",
        full.names = TRUE, ignore.case = TRUE
      )
      if (!length(shard_paths)) {
        stop("variant-set directory has no chr*.tsv or chr*.tsv.gz shards: ",
             variant_set, call. = FALSE)
      }
      shard_chromosomes <- sub(
        "^chr", "", sub("[.]tsv(?:[.]gz)?$", "", basename(shard_paths),
                           ignore.case = TRUE), ignore.case = TRUE
      )
      chromosome_order <- c(as.character(1:22), "X", "Y")
      shard_paths <- shard_paths[order(match(shard_chromosomes, chromosome_order),
                                       shard_chromosomes)]
      shard_chromosomes <- shard_chromosomes[order(match(shard_chromosomes, chromosome_order),
                                                   shard_chromosomes)]
      requested_chromosomes <- if (is.null(chromosomes)) {
        NULL
      } else {
        unique(normalise_chromosome(chromosomes))
      }
      selected <- if (is.null(requested_chromosomes)) {
        rep(TRUE, length(shard_paths))
      } else {
        shard_chromosomes %in% requested_chromosomes
      }
      shard_paths <- shard_paths[selected]
      shard_chromosomes <- shard_chromosomes[selected]
      if (!length(shard_paths)) {
        stop("variant-set directory has no shards for the requested chromosomes",
             call. = FALSE)
      }
      pieces <- lapply(shard_paths, read_variant_set)
      canonical <- all(vapply(pieces, function(x) {
        isTRUE(attr(x, "variant_set_canonical"))
      }, logical(1L)))
      out <- if (canonical) {
        data.frame(
          variant_id = unlist(lapply(pieces, function(x) x$variant_id),
                               use.names = FALSE),
          stringsAsFactors = FALSE
        )
      } else {
        do.call(rbind, pieces)
      }
      if (canonical) {
        attr(out, "variant_set_normalized_ids") <- TRUE
        attr(out, "variant_set_canonical") <- TRUE
      }
      manifest_path <- file.path(dirname(shard_dir), "panel_manifest.json")
      manifest <- if (file.exists(manifest_path) && requireNamespace("jsonlite", quietly = TRUE)) {
        tryCatch(jsonlite::read_json(manifest_path, simplifyVector = FALSE),
                 error = function(e) NULL)
      } else {
        NULL
      }
      panel_name <- sub("_by_chrom$", "", basename(shard_dir))
      metadata <- variant_set_provenance(
        path = variant_set, name = named$name %||% NULL, named = named,
        rows = nrow(out), columns_read = "variant_id",
        identity = if (canonical) "canonical_allele_aware" else "legacy",
        manifest = manifest
      )
      metadata$panel_name <- panel_name
      metadata$chromosomes <- shard_chromosomes
      metadata$shards <- basename(shard_paths)
      attr(out, "variant_set_metadata") <- metadata
      return(out)
    }
    columns_read <- NULL
    if (grepl("[.]bim(?:[.]gz)?$", variant_set, ignore.case = TRUE)) {
      bim_args <- list(data.table = FALSE, header = FALSE,
                       col.names = c("chromosome", "variant_id", "genetic_distance",
                                     "base_pair_location", "other_allele", "effect_allele"),
                       showProgress = FALSE)
      if (grepl("[.]gz$", variant_set, ignore.case = TRUE)) {
        bim_args$cmd <- compressed_read_command(variant_set)
        out <- do.call(data.table::fread, bim_args)
      } else {
        bim_args$input <- variant_set
        out <- do.call(data.table::fread, bim_args)
      }
    } else if (grepl("[.]parquet$", variant_set, ignore.case = TRUE)) {
      out <- arrow::read_parquet(variant_set)
    } else {
      out <- read_delimited_variant_set(variant_set)
    }
    columns_read <- attr(out, "variant_set_columns_read")
    canonical_ids <- attr(out, "variant_set_canonical_ids")
    if (!is.null(canonical_ids) && nrow(out) > 100000L) {
      # Large canonical dictionaries are already allele-aware identity keys.
      # Do not expand six million keys into derived chromosome, position, and
      # allele columns: membership only needs the key, and the regex/substr
      # work otherwise dominates panel preparation. Small panels retain the
      # richer historical representation for callers that inspect them.
      out[[1L]] <- as.character(canonical_ids)
      names(out)[1L] <- "variant_id"
      attr(out, "variant_set_canonical_ids") <- NULL
      attr(out, "variant_set_normalized_ids") <- TRUE
      attr(out, "variant_set_canonical") <- TRUE
    } else {
      out <- normalise_variant_set_columns(out)
      if (!is.null(canonical_ids)) {
        out$variant_id <- normalise_variant_key(out$variant_id)
        attr(out, "variant_set_normalized_ids") <- TRUE
        attr(out, "variant_set_canonical") <- TRUE
      }
    }
    metadata <- variant_set_provenance(
      path = variant_set, name = named$name %||% NULL, named = named,
      rows = nrow(out), columns_read = columns_read %||% names(out),
      identity = if (isTRUE(attr(out, "variant_set_canonical"))) {
        "canonical_allele_aware"
      } else {
        "legacy"
      }
    )
  }
  attr(out, "variant_set_metadata") <- metadata
  out
}

normalise_variant_key <- function(x) {
  x <- as.character(x)
  toupper(sub("^chr", "", x, ignore.case = TRUE))
}

variant_set_membership <- function(data, panel) {
  panel_values <- panel$variant_id[!is.na(panel$variant_id) & nzchar(panel$variant_id)]
  panel_id <- if (isTRUE(attr(panel, "variant_set_normalized_ids"))) {
    unique(as.character(panel_values))
  } else {
    unique(normalise_variant_key(panel_values))
  }
  data_id <- if ("variant_id" %in% names(data)) {
    normalise_variant_key(data$variant_id)
  } else {
    rep(NA_character_, nrow(data))
  }
  data_key <- rep(NA_character_, nrow(data))
  data_chromosome <- data$chromosome %||% rep(NA_character_, nrow(data))
  data_position <- data$base_pair_location %||% rep(NA_real_, nrow(data))
  data_other <- data$other_allele %||% rep(NA_character_, nrow(data))
  data_effect <- data$effect_allele %||% rep(NA_character_, nrow(data))
  valid_identity <- data_chromosome %in% c(as.character(1:22), "X", "Y") &
    is.finite(data_position) & data_position >= 1L &
    !is.na(data_other) & data_other %in% c("A", "C", "G", "T") &
    !is.na(data_effect) & data_effect %in% c("A", "C", "G", "T") &
    data_other != data_effect
  if (any(valid_identity)) {
    data_key[valid_identity] <- compressor_variant_key(
      data_chromosome[valid_identity], data_position[valid_identity],
      data_other[valid_identity], data_effect[valid_identity]
    )
  }
  # Canonical panels are allele-aware identity dictionaries. Do not fall back
  # to coordinate or rsID matching, which could admit the wrong allele and
  # forces needless construction of large coordinate vectors.
  data_id_match <- data_id %in% panel_id & (!valid_identity | data_key == data_id)
  if (isTRUE(attr(panel, "variant_set_canonical"))) {
    return(data_id_match | data_key %in% panel_id)
  }
  canonical_id <- grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$",
                         panel_id)
  if (any(canonical_id)) return(data_id_match & data_id %in% panel_id[canonical_id] |
                                data_key %in% panel_id[canonical_id])
  panel_rsid <- unique(normalise_variant_key(panel$rsid[!is.na(panel$rsid) & nzchar(panel$rsid)]))
  panel_rsid <- unique(c(panel_rsid, panel_id[grepl("^rs", panel_id, ignore.case = TRUE)]))
  keep <- data_id %in% panel_id
  if ("rsid" %in% names(data)) keep <- keep | normalise_variant_key(data$rsid) %in% panel_rsid
  panel_coord <- paste(panel$chromosome, panel$base_pair_location, sep = ":")
  data_coord <- paste(data_chromosome, data_position, sep = ":")
  keep | (!is.na(data_position) & data_coord %in% panel_coord)
}

validate_pvalue_region_arguments <- function(pvalue_threshold, region_padding) {
  if (length(pvalue_threshold) != 1L || !is.numeric(pvalue_threshold) ||
      !is.finite(pvalue_threshold) || pvalue_threshold <= 0 ||
      pvalue_threshold > 1) {
    stop("pvalue_threshold must be one finite number in (0, 1]", call. = FALSE)
  }
  if (length(region_padding) != 1L || !is.numeric(region_padding) ||
      !is.finite(region_padding) || region_padding < 0 ||
      region_padding != floor(region_padding) ||
      region_padding > .Machine$integer.max) {
    stop("region_padding must be one non-negative whole number", call. = FALSE)
  }
  invisible(TRUE)
}

merge_pvalue_regions <- function(seeds) {
  if (!nrow(seeds)) {
    return(data.frame(chromosome = character(), start = integer(), end = integer(),
                      seed_snps = integer(), stringsAsFactors = FALSE))
  }
  chromosome_order <- c(as.character(1:22), "X", "Y", "MT")
  split_seeds <- split(seeds, seeds$chromosome, drop = TRUE)
  split_seeds <- split_seeds[order(match(names(split_seeds), chromosome_order),
                                   names(split_seeds), na.last = TRUE)]
  merged <- lapply(split_seeds, function(group) {
    group <- group[order(group$start, group$end), , drop = FALSE]
    starts <- integer()
    ends <- integer()
    counts <- integer()
    current_start <- group$start[[1L]]
    current_end <- group$end[[1L]]
    current_count <- 1L
    if (nrow(group) > 1L) {
      for (row in seq.int(2L, nrow(group))) {
        if (group$start[[row]] <= current_end + 1) {
          current_end <- max(current_end, group$end[[row]])
          current_count <- current_count + 1L
        } else {
          starts <- c(starts, current_start)
          ends <- c(ends, current_end)
          counts <- c(counts, current_count)
          current_start <- group$start[[row]]
          current_end <- group$end[[row]]
          current_count <- 1L
        }
      }
    }
    data.frame(chromosome = group$chromosome[[1L]],
               start = c(starts, current_start),
               end = c(ends, current_end),
               seed_snps = c(counts, current_count),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, merged)
  row.names(out) <- NULL
  out
}

pvalue_region_selection <- function(data, pvalue_threshold = 1e-5,
                                     region_padding = 50000L) {
  validate_pvalue_region_arguments(pvalue_threshold, region_padding)
  if (!"z" %in% names(data)) {
    stop("p-value region selection needs the pre-encoding harmonised z statistic",
         call. = FALSE)
  }
  n <- nrow(data)
  chromosomes <- as.character(data$chromosome)
  positions <- as.numeric(data$base_pair_location)
  # This is deliberately evaluated on the harmonised in-memory statistic. It
  # runs before any bounded semantic encoding and is never reconstructed from
  # a lossy stored Z value.
  p_value <- 2 * stats::pnorm(-abs(as.numeric(data$z)))
  valid_coordinate <- !is.na(chromosomes) & nzchar(chromosomes) &
    is.finite(positions) & positions >= 1
  seed <- valid_coordinate & is.finite(p_value) & p_value < pvalue_threshold
  seeds <- data.frame(
    chromosome = chromosomes[seed],
    start = as.integer(pmax(1, positions[seed] -
                             as.numeric(region_padding))),
    end = as.integer(pmin(.Machine$integer.max, positions[seed] +
                            as.numeric(region_padding))),
    stringsAsFactors = FALSE
  )
  regions <- merge_pvalue_regions(seeds)
  keep <- rep(FALSE, n)
  if (nrow(regions)) {
    for (region_chromosome in unique(regions$chromosome)) {
      rows <- which(valid_coordinate & chromosomes == region_chromosome)
      if (!length(rows)) next
      region_rows <- regions[regions$chromosome == region_chromosome, , drop = FALSE]
      index <- findInterval(positions[rows], region_rows$start)
      inside <- index > 0L & positions[rows] <=
        region_rows$end[pmax(1L, index)]
      keep[rows[inside]] <- TRUE
    }
  }
  metadata <- list(
    tag = "core",
    method = "pvalue_regions",
    p_value_source = "derived_from_z",
    threshold_source = "pre_encoding_harmonised",
    threshold_statistic = "p_value_from_harmonised_z",
    threshold_operator = "<",
    threshold_semantics = list(
      source = "pre_encoding_harmonised",
      statistic = "p_value_from_harmonised_z",
      derivation = "2 * pnorm(-abs(z))",
      operator = "<",
      value = as.numeric(pvalue_threshold),
      encoding = "not_encoded"
    ),
    pvalue_threshold = as.numeric(pvalue_threshold),
    padding_bp = as.integer(region_padding),
    seed_snps = as.integer(sum(seed)),
    regions = as.integer(nrow(regions)),
    input_rows = as.integer(n),
    kept_rows = as.integer(sum(keep)),
    dropped_rows = as.integer(sum(!keep))
  )
  list(keep = keep, regions = regions, metadata = metadata)
}

empty_selection_regions <- function() {
  data.frame(chromosome = character(), start = integer(), end = integer(),
             seed_snps = integer(), stringsAsFactors = FALSE)
}

compact_panel_provenance <- function(panel) {
  metadata <- attr(panel, "variant_set_metadata") %||% list()
  fields <- c("id", "name", "requested", "source", "local_path", "source_url",
              "cache_dir", "cache_hit", "build", "version", "source_sha256",
              "sha256", "expected_sha256", "hash_pinned", "hash_verified",
              "digest_target", "identity", "rows", "columns_read", "panel_name",
              "chromosomes", "shards")
  metadata[intersect(fields, names(metadata))]
}

panel_selection_result <- function(data, panel, selection_name,
                                   panel_keep = NULL, method = "panel_membership") {
  panel_keep <- panel_keep %||% variant_set_membership(data, panel)
  if (!any(panel_keep)) {
    stop(selection_name, " selection retained no variants", call. = FALSE)
  }
  panel_metadata <- compact_panel_provenance(panel)
  panel_name <- panel_metadata$name %||% "custom"
  metadata <- list(
    selection = selection_name,
    tag = selection_name,
    method = method,
    input_rows = as.integer(nrow(data)),
    kept_rows = as.integer(sum(panel_keep)),
    dropped_rows = as.integer(sum(!panel_keep)),
    panel_name = panel_name,
    panel_hash = panel_metadata$sha256 %||% NA_character_,
    panel_hash_pinned = isTRUE(panel_metadata$hash_pinned),
    panel_provenance = panel_metadata
  )
  list(
    data = data[panel_keep, , drop = FALSE],
    keep = panel_keep,
    panel = panel,
    regions = empty_selection_regions(),
    metadata = metadata
  )
}

select_full_variants <- function(data) {
  keep <- rep(TRUE, nrow(data))
  list(
    data = data,
    keep = keep,
    panel = NULL,
    regions = empty_selection_regions(),
    metadata = list(
      selection = "full", tag = "full", method = "all_variants",
      input_rows = as.integer(nrow(data)), kept_rows = as.integer(nrow(data)),
      dropped_rows = 0L
    )
  )
}

select_core_variants <- function(data, variant_set = "core") {
  panel <- read_variant_set(variant_set, chromosomes = unique(data$chromosome))
  panel_selection_result(data, panel, "core")
}

select_hm3_variants <- function(data, variant_set = "hm3") {
  panel <- read_variant_set(variant_set, chromosomes = unique(data$chromosome))
  panel_selection_result(data, panel, "hm3")
}

select_core_plus_variants <- function(data, variant_set = "core",
                                      pvalue_threshold = 1e-5,
                                      region_padding = 50000L) {
  panel <- read_variant_set(variant_set, chromosomes = unique(data$chromosome))
  panel_keep <- variant_set_membership(data, panel)
  regions <- pvalue_region_selection(data, pvalue_threshold = pvalue_threshold,
                                     region_padding = region_padding)
  keep <- panel_keep | regions$keep
  if (!any(keep)) {
    stop("core_plus selection retained no variants; check the core panel or adjust pvalue_threshold/region_padding",
         call. = FALSE)
  }
  panel_metadata <- compact_panel_provenance(panel)
  panel_name <- panel_metadata$name %||% "custom"
  # Update the region metadata in place. Appending a second `tag`/`method`
  # entry leaves duplicate list names, and `$tag` then returns the earlier
  # p-value-only `core` value when the selection is written to a manifest.
  metadata <- regions$metadata
  metadata$selection <- "core_plus"
  metadata$tag <- "core_plus"
  metadata$method <- "panel_union_pvalue_regions"
  metadata$core_variant_rows <- as.integer(sum(panel_keep))
  metadata$panel_name <- panel_name
  metadata$panel_hash <- panel_metadata$sha256 %||% NA_character_
  metadata$panel_hash_pinned <- isTRUE(panel_metadata$hash_pinned)
  metadata$panel_provenance <- panel_metadata
  metadata$kept_rows <- as.integer(sum(keep))
  metadata$dropped_rows <- as.integer(sum(!keep))
  list(
    data = data[keep, , drop = FALSE],
    keep = keep,
    panel = panel,
    regions = regions$regions,
    metadata = metadata
  )
}

select_pvalue_regions <- function(data, pvalue_threshold = 1e-5,
                                  region_padding = 50000L) {
  regions <- pvalue_region_selection(data, pvalue_threshold = pvalue_threshold,
                                     region_padding = region_padding)
  if (!regions$metadata$kept_rows) {
    stop("pvalue_regions selection retained no variants; adjust pvalue_threshold or region_padding",
         call. = FALSE)
  }
  regions$metadata$selection <- "pvalue_regions"
  regions$metadata$regions_table <- regions$regions
  list(
    data = data[regions$keep, , drop = FALSE],
    keep = regions$keep,
    panel = NULL,
    regions = regions$regions,
    metadata = regions$metadata
  )
}

select_variant_rows <- function(data, selection = c("full", "core", "hm3", "core_plus"),
                                variant_set = NULL, pvalue_threshold = 1e-5,
                                region_padding = 50000L) {
  selection <- match.arg(selection)
  switch(selection,
    full = select_full_variants(data),
    core = select_core_variants(data, variant_set = variant_set %||% "core"),
    hm3 = select_hm3_variants(data, variant_set = variant_set %||% "hm3"),
    core_plus = select_core_plus_variants(
      data, variant_set = variant_set %||% "core",
      pvalue_threshold = pvalue_threshold, region_padding = region_padding
    )
  )
}

# Descriptive aliases used by callers that treat selection as a public
# preparation step. They all return the same stable list(data, keep, panel,
# regions, metadata) contract.
full_variant_selection <- select_full_variants
core_variant_selection <- select_core_variants
hm3_variant_selection <- select_hm3_variants
core_plus_variant_selection <- select_core_plus_variants
select_full <- select_full_variants
select_core <- select_core_variants
select_hm3 <- select_hm3_variants
select_core_plus <- select_core_plus_variants
select_variant_set <- select_variant_rows

read_delimited_variant_set <- function(path) {
  command <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    compressed_read_command(path)
  } else {
    NULL
  }
  header_args <- list(data.table = FALSE, nrows = 0L, showProgress = FALSE)
  header_args$check.names <- FALSE
  if (is.null(command)) header_args$input <- path else header_args$cmd <- command
  header <- do.call(data.table::fread, header_args)
  alias_groups <- list(
    variant_id = c("variant_id", "variant_id_grch38", "vid", "SNPID", "SNP", "snp", "ID", "id"),
    rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID"),
    chromosome = c("chromosome", "chr", "CHR", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "BP", "bp"),
    other_allele = c("other_allele", "NEA", "REF", "ref", "A2"),
    effect_allele = c("effect_allele", "EA", "ALT", "alt", "A1")
  )
  selected <- unname(vapply(alias_groups, function(options) {
    hit <- intersect(options, names(header))
    if (length(hit)) hit[[1L]] else NA_character_
  }, character(1L)))
  selected <- selected[!is.na(selected)]
  if (!length(selected)) {
    stop("variant set has no recognizable identity columns: ", path, call. = FALSE)
  }
  # Canonical dictionaries only need their allele-aware key for membership.
  # Probe that one column first; this avoids materialising chromosome, position,
  # allele and rsID columns for the multi-million-row core/HM3 files. Legacy
  # rsID/coordinate panels fall through to the small identity-column set.
  canonical_pattern <- "^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$"
  variant_id_column <- selected[selected %in% alias_groups$variant_id][1L]
  if (!is.na(variant_id_column)) {
    id_args <- list(data.table = FALSE, showProgress = FALSE, check.names = FALSE,
                    select = variant_id_column)
    if (is.null(command)) id_args$input <- path else id_args$cmd <- command
    id_data <- do.call(data.table::fread, id_args)
    id_values <- toupper(sub("^chr", "", trimws(as.character(id_data[[variant_id_column]])),
                          ignore.case = TRUE))
    present <- !is.na(id_values) & nzchar(id_values) & id_values != "."
    if (any(present) && all(grepl(canonical_pattern, id_values[present]))) {
      # Retain the normalized values in the one column already read. This
      # avoids keeping both the raw and normalized six-million-row strings
      # alive while the caller decides whether the large-panel fast path is
      # applicable.
      id_data[[variant_id_column]] <- id_values
      names(id_data)[1L] <- "variant_id"
      attr(id_data, "variant_set_canonical_ids") <- id_values
      attr(id_data, "variant_set_columns_read") <- variant_id_column
      return(id_data)
    }
  }
  read_args <- list(data.table = FALSE, showProgress = FALSE, check.names = FALSE,
                    select = selected)
  if (is.null(command)) read_args$input <- path else read_args$cmd <- command
  out <- do.call(data.table::fread, read_args)
  attr(out, "variant_set_columns_read") <- selected
  out
}

prepare_sumstats_data <- function(raw, reference, mode = c("qc", "convert", "all", "core", "hm3", "pvalue_regions", "core_plus"),
                                  variant_set = NULL, strict = FALSE, chrom_threads = 1L,
                                  drop_unresolved = TRUE, input_build = "GRCh38", chain = NULL,
                                  pvalue_threshold = 1e-5, region_padding = 50000L,
                                  pre_harmonised = FALSE, inherited_alignment = NULL,
                                  target_build = "GRCh38") {
  observability <- .compressor_harmonise_context$observability %||% NULL
  requested_mode <- match.arg(mode)
  effective_mode <- if (requested_mode == "all") "qc" else requested_mode
  target_build <- compressor_normalize_build(target_build)
  input_build <- compressor_normalize_build(input_build)
  if (effective_mode %in% c("core", "hm3", "core_plus") && is.null(variant_set)) {
    # The explicit named defaults are resolved lazily. This keeps ordinary
    # full/QC conversion independent of external panel assets while making
    # core/HM3 selection deterministic when the corresponding named panel is
    # configured.
    variant_set <- if (effective_mode == "hm3") "hm3" else "core"
  }
  alignment_reference <- if (isTRUE(pre_harmonised)) {
    NULL
  } else if (effective_mode == "convert" ||
                            (effective_mode == "pvalue_regions" && is.null(reference))) {
    NULL
  } else {
    reference
  }
  if (!isTRUE(pre_harmonised) && effective_mode != "convert" && effective_mode != "pvalue_regions" &&
      is.null(alignment_reference)) {
    stop("reference is required for harmonised CompreSSoR conversion; use mode='convert' explicitly to bypass it",
         call. = FALSE)
  }
  if (isTRUE(pre_harmonised)) {
    # A harmonise_sumstats() result already has canonical coordinates, alleles,
    # orientation and Z. Keep that object on the fast path: panel/p-value
    # selection is post-harmonisation and no reference or coordinate work is
    # repeated here.
    alignment <- inherited_alignment %||% list()
    alignment$data <- raw
    alignment$reference_hash <- alignment$reference_hash %||% NA_character_
    alignment$reference_rows <- alignment$reference_rows %||% NA_integer_
    alignment$reference_metadata <- alignment$reference_metadata %||% NULL
    alignment$alignment_stats <- alignment$alignment_stats %||% list(
      input_rows = nrow(raw), output_rows = nrow(raw), aligned_rows = nrow(raw),
      method = "pre_harmonised"
    )
  } else {
    had_coordinate <- !is.na(raw$chromosome) & nzchar(raw$chromosome) &
      !is.na(raw$base_pair_location) & raw$base_pair_location >= 1L
    liftover_phase <- compressor_observability_start_phase(
      observability, "liftover", rows_in = nrow(raw)
    )
    if (identical(input_build, target_build)) {
      raw$.compressor_liftover_status <- rep("not_needed", nrow(raw))
      attr(raw, "genome_build") <- target_build
    } else if (identical(target_build, "GRCh38")) {
      raw <- lift_table_to_grch38(raw, input_build = input_build, chain = chain)
    } else {
      stop("target_build='GRCh37' requires input already in GRCh37; reverse liftover is not supported", call. = FALSE)
    }
    compressor_observability_finish_phase(
      observability, liftover_phase, rows_out = nrow(raw)
    )
    if (!is.null(alignment_reference)) {
      build <- toupper(gsub("[. -]", "", as.character(input_build %||% "GRCh38")))
      eligible_alias <- if (build %in% c("GRCH38", "HG38", "38")) {
        rep(TRUE, nrow(raw))
      } else {
        !had_coordinate
      }
      raw <- resolve_sumstats_aliases(raw, alignment_reference, eligible = eligible_alias)
      resolved_alias <- raw$.compressor_reference_alias_status == "resolved"
      raw$.compressor_liftover_status[resolved_alias] <- "resolved_by_grch38_alias"
    }
    partitioned_reference <- !is.null(alignment_reference) &&
      is_partitioned_reference(alignment_reference)
    alignment <- if (!is.null(alignment_reference) &&
                    (as.integer(chrom_threads) > 1L || partitioned_reference)) {
      harmonise_by_chromosome(raw, alignment_reference, strict = strict,
                              preserve = !isTRUE(drop_unresolved),
                              chrom_threads = chrom_threads)
    } else {
      harmonise_to_reference(raw, alignment_reference, strict = strict,
                             preserve = !isTRUE(drop_unresolved))
    }
  }
  data <- alignment$data
  if (!isTRUE(pre_harmonised)) {
    alias_status <- raw$.compressor_reference_alias_status %||% rep("not_attempted", nrow(raw))
    liftover_status <- raw$.compressor_liftover_status %||% rep("not_needed", nrow(raw))
    alignment$alignment_stats$alias_resolution <- as.list(table(alias_status, useNA = "ifany"))
    alignment$alignment_stats$liftover <- as.list(table(liftover_status, useNA = "ifany"))
  }
  input_rows <- nrow(data)
  selected <- if (effective_mode == "core") {
    select_core_variants(data, variant_set = variant_set)
  } else if (effective_mode == "hm3") {
    select_hm3_variants(data, variant_set = variant_set)
  } else if (effective_mode == "core_plus") {
    select_core_plus_variants(data, variant_set = variant_set,
                              pvalue_threshold = pvalue_threshold,
                              region_padding = region_padding)
  } else if (effective_mode == "pvalue_regions") {
    select_pvalue_regions(data, pvalue_threshold = pvalue_threshold,
                          region_padding = region_padding)
  } else if (!is.null(variant_set)) {
    panel <- read_variant_set(variant_set, chromosomes = unique(data$chromosome))
    panel_selection_result(data, panel, "full", method = "panel_membership")
  } else {
    select_full_variants(data)
  }
  data <- selected$data
  panel_metadata <- if (!is.null(selected$panel)) {
    attr(selected$panel, "variant_set_metadata") %||% list()
  } else {
    NULL
  }
  filter_stats <- if (is.null(panel_metadata)) {
    list(name = "all", input_rows = input_rows, kept_rows = nrow(data),
         dropped_rows = input_rows - nrow(data))
  } else {
    c(
      list(name = panel_metadata$name %||% effective_mode),
      panel_metadata,
      list(panel_name = selected$metadata$panel_name %||% "custom",
           input_rows = input_rows, kept_rows = sum(selected$keep),
           dropped_rows = sum(!selected$keep))
    )
  }
  alignment$alignment_stats$variant_set <- filter_stats
  selection <- NULL
  if (effective_mode %in% c("pvalue_regions", "core_plus")) {
    selection <- selected$metadata
    selection$regions_table <- selected$regions
    selection_stats <- selected$metadata
    selection_stats$regions_table <- NULL
    alignment$alignment_stats$pvalue_regions <- selection_stats
  }
  list(data = data, alignment = alignment, requested_mode = requested_mode,
       effective_mode = effective_mode,
       selection = selection,
       selection_result = selected,
    genome_build = if (isTRUE(pre_harmonised)) {
         inherited_alignment$genome_build %||% target_build
       } else if (effective_mode == "convert" || is.null(reference)) {
         # A reference-free conversion normalises the schema but cannot prove
         # the assembly. Keep the build unknown so a caller cannot silently
         # publish GRCh38 coordinates as if they had been verified.
         "unknown"
       } else {
         alignment$reference_metadata$build %||% "GRCh38"
       })
}
