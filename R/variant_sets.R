resolve_named_variant_set <- function(variant_set) {
  if (!is.character(variant_set) || length(variant_set) != 1L ||
      is.na(variant_set)) return(NULL)
  name <- tolower(trimws(variant_set))
  aliases <- c(common = "common", common_variants = "common",
               core = "core", core_variants = "core",
               tag = "tag", tags = "tag", tag_variants = "tag",
               hm3 = "hm3", hm3_variants = "hm3")
  canonical_name <- unname(aliases[name])
  if (is.na(canonical_name)) return(NULL)
  env_name <- switch(canonical_name,
    common = "COMPRESSOR_COMMON_VARIANTS",
    core = "COMPRESSOR_CORE_VARIANTS",
    tag = "COMPRESSOR_TAG_VARIANTS",
    hm3 = "COMPRESSOR_HM3_VARIANTS"
  )
  path <- Sys.getenv(env_name, unset = "")
  if (!nzchar(path)) {
    alternate_env <- paste0("COMPRESSOR_VARIANT_SET_", toupper(canonical_name))
    path <- Sys.getenv(alternate_env, unset = "")
    if (nzchar(path)) env_name <- alternate_env
  }
  if (!nzchar(path)) {
    root <- Sys.getenv("COMPRESSOR_VARIANT_SET_DIR", unset = "")
    if (nzchar(root)) {
      candidates <- switch(canonical_name,
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
    stop("variant_set='", canonical_name, "' needs ", env_name,
         " (or COMPRESSOR_VARIANT_SET_DIR with a named panel file)", call. = FALSE)
  }
  path <- path.expand(path)
  if (!file.exists(path)) {
    stop(env_name, " points to a missing variant-set file: ", path, call. = FALSE)
  }
  list(name = canonical_name, path = path, env_var = env_name)
}
normalise_variant_set_columns <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  alias_map <- list(
    variant_id = c("variant_id", "variant_id_grch38", "vid", "SNPID", "SNP", "snp", "ID", "id"),
    rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID", "ID"),
    chromosome = c("chromosome", "chr", "CHR", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "BP", "bp"),
    other_allele = c("other_allele", "NEA", "REF", "ref", "A2"),
    effect_allele = c("effect_allele", "EA", "ALT", "alt", "A1")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  if (!"variant_id" %in% names(data) && all(c("chromosome", "base_pair_location", "other_allele", "effect_allele") %in% names(data))) {
    data$variant_id <- paste(
      sub("^chr", "", as.character(data$chromosome), ignore.case = TRUE),
      data$base_pair_location, toupper(data$other_allele), toupper(data$effect_allele), sep = ":"
    )
  }
  if (!"variant_id" %in% names(data) && !"rsid" %in% names(data)) {
    stop("variant set needs variant_id, rsid, or chromosome/position/allele columns", call. = FALSE)
  }
  if (!"variant_id" %in% names(data)) data$variant_id <- NA_character_
  if (!"rsid" %in% names(data)) data$rsid <- NA_character_
  if (!"chromosome" %in% names(data)) data$chromosome <- NA_character_
  if (!"base_pair_location" %in% names(data)) data$base_pair_location <- NA_integer_
  canonical_id <- toupper(sub("^chr", "", as.character(data$variant_id), ignore.case = TRUE))
  canonical <- grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$",
                      canonical_id)
  if (any(canonical)) {
    parsed_chromosome <- sub(":.*$", "", canonical_id)
    parsed_position <- suppressWarnings(as.integer(sub("^[^:]+:([^:]+):.*$", "\\1", canonical_id)))
    fill_chromosome <- canonical & (is.na(data$chromosome) | !nzchar(data$chromosome))
    fill_position <- canonical & is.na(data$base_pair_location)
    data$chromosome[fill_chromosome] <- parsed_chromosome[fill_chromosome]
    data$base_pair_location[fill_position] <- parsed_position[fill_position]
  }
  data$variant_id <- as.character(data$variant_id)
  data$rsid <- as.character(data$rsid)
  data$chromosome <- sub("^chr", "", as.character(data$chromosome), ignore.case = TRUE)
  data$base_pair_location <- suppressWarnings(as.integer(data$base_pair_location))
  data
}

read_variant_set <- function(variant_set, chromosomes = NULL) {
  named <- resolve_named_variant_set(variant_set)
  if (!is.data.frame(variant_set) && !is.null(named)) variant_set <- named$path
  if (is.data.frame(variant_set)) {
    out <- normalise_variant_set_columns(variant_set)
    metadata <- list(id = "in_memory", rows = nrow(out))
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
      metadata <- list(
        id = named$name %||% "local",
        name = named$name %||% NULL,
        local_path = normalizePath(variant_set, mustWork = FALSE),
        sha256 = reference_file_digest(file.path(variant_set, "manifest.json"), "sha256"),
        rows = nrow(out),
        columns_read = c("chromosome", "base_pair_location", "effect_allele", "other_allele")
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
      panel_metadata <- manifest$panels[[panel_name]] %||% list()
      metadata <- list(
        id = named$name %||% "local",
        name = named$name %||% NULL,
        local_path = normalizePath(variant_set, mustWork = FALSE),
        sha256 = panel_metadata$prepared_sha256 %||% NA_character_,
        rows = nrow(out),
        columns_read = "variant_id",
        chromosomes = shard_chromosomes,
        shards = basename(shard_paths)
      )
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
    metadata <- list(
      id = named$name %||% "local",
      name = named$name %||% NULL,
      local_path = normalizePath(variant_set, mustWork = FALSE),
      sha256 = reference_file_digest(variant_set, "sha256"),
      rows = nrow(out),
      columns_read = columns_read %||% names(out)
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
  data_id <- normalise_variant_key(data$variant_id)
  data_key <- rep(NA_character_, nrow(data))
  valid_identity <- !is.na(data$chromosome) & nzchar(data$chromosome) &
    is.finite(data$base_pair_location) & data$base_pair_location >= 1L &
    !is.na(data$other_allele) & data$other_allele %in% c("A", "C", "G", "T") &
    !is.na(data$effect_allele) & data$effect_allele %in% c("A", "C", "G", "T") &
    data$other_allele != data$effect_allele
  if (any(valid_identity)) {
    data_key[valid_identity] <- compressor_variant_key(
      data$chromosome[valid_identity], data$base_pair_location[valid_identity],
      data$other_allele[valid_identity], data$effect_allele[valid_identity]
    )
  }
  # Canonical panels are allele-aware identity dictionaries. Do not fall back
  # to coordinate or rsID matching, which could admit the wrong allele and
  # forces needless construction of large coordinate vectors.
  if (isTRUE(attr(panel, "variant_set_canonical"))) {
    return(data_id %in% panel_id | data_key %in% panel_id)
  }
  canonical_id <- grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$",
                         panel_id)
  if (any(canonical_id)) return(data_id %in% panel_id[canonical_id] |
                                data_key %in% panel_id[canonical_id])
  panel_rsid <- unique(normalise_variant_key(panel$rsid[!is.na(panel$rsid) & nzchar(panel$rsid)]))
  panel_rsid <- unique(c(panel_rsid, panel_id[grepl("^rs", panel_id, ignore.case = TRUE)]))
  keep <- data_id %in% panel_id
  if ("rsid" %in% names(data)) keep <- keep | normalise_variant_key(data$rsid) %in% panel_rsid
  panel_coord <- paste(panel$chromosome, panel$base_pair_location, sep = ":")
  data_coord <- paste(data$chromosome, data$base_pair_location, sep = ":")
  keep | (!is.na(data$base_pair_location) & data_coord %in% panel_coord)
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
  n <- nrow(data)
  chromosomes <- as.character(data$chromosome)
  positions <- as.numeric(data$base_pair_location)
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

read_delimited_variant_set <- function(path) {
  command <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    compressed_read_command(path)
  } else {
    NULL
  }
  header_args <- list(data.table = FALSE, nrows = 0L, showProgress = FALSE)
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
    id_args <- list(data.table = FALSE, showProgress = FALSE,
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
  read_args <- list(data.table = FALSE, showProgress = FALSE, select = selected)
  if (is.null(command)) read_args$input <- path else read_args$cmd <- command
  out <- do.call(data.table::fread, read_args)
  attr(out, "variant_set_columns_read") <- selected
  out
}

prepare_sumstats_data <- function(raw, reference, mode = c("qc", "convert", "all", "core", "hm3", "pvalue_regions", "core_plus"),
                                  variant_set = NULL, strict = FALSE, chrom_threads = 1L,
                                  drop_unresolved = TRUE, input_build = "GRCh38", chain = NULL,
                                  pvalue_threshold = 1e-5, region_padding = 50000L,
                                  pre_harmonised = FALSE, inherited_alignment = NULL) {
  requested_mode <- match.arg(mode)
  effective_mode <- if (requested_mode == "all") "qc" else requested_mode
  if (effective_mode %in% c("core", "hm3", "core_plus") && is.null(variant_set)) {
    stop("mode='", requested_mode, "' requires variant_set = a core/HM3 panel file or data.frame", call. = FALSE)
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
    raw <- lift_table_to_grch38(raw, input_build = input_build, chain = chain)
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
  filter_stats <- list(name = "all", input_rows = nrow(data), kept_rows = nrow(data), dropped_rows = 0L)
  panel_keep <- rep(TRUE, nrow(data))
  if (!is.null(variant_set)) {
    panel <- read_variant_set(variant_set, chromosomes = unique(data$chromosome))
    panel_keep <- variant_set_membership(data, panel)
    filter_stats <- c(
      list(name = attr(panel, "variant_set_metadata")$name %||% effective_mode),
      attr(panel, "variant_set_metadata") %||% list(),
      list(input_rows = nrow(data), kept_rows = sum(panel_keep), dropped_rows = sum(!panel_keep))
    )
    if (effective_mode %in% c("core", "hm3") ||
        !effective_mode %in% c("core_plus", "pvalue_regions")) {
      data <- data[panel_keep, , drop = FALSE]
      if (!nrow(data)) stop("variant_set filter retained no variants", call. = FALSE)
    }
  }
  alignment$alignment_stats$variant_set <- filter_stats
  selection <- NULL
  if (effective_mode %in% c("pvalue_regions", "core_plus")) {
    selected <- pvalue_region_selection(data, pvalue_threshold = pvalue_threshold,
                                        region_padding = region_padding)
    if (effective_mode == "pvalue_regions" && !selected$metadata$kept_rows) {
      stop("pvalue_regions selection retained no variants; adjust pvalue_threshold or region_padding",
           call. = FALSE)
    }
    if (effective_mode == "core_plus") {
      keep <- panel_keep | selected$keep
      if (!any(keep)) {
        stop("core_plus selection retained no variants; check the core panel or adjust pvalue_threshold/region_padding",
             call. = FALSE)
      }
      selected$metadata$tag <- "core_plus"
      selected$metadata$method <- "core_plus_pvalue_regions"
      selected$metadata$core_variant_rows <- as.integer(sum(panel_keep))
      selected$metadata$kept_rows <- as.integer(sum(keep))
      selected$metadata$dropped_rows <- as.integer(sum(!keep))
      data <- data[keep, , drop = FALSE]
    } else {
      data <- data[selected$keep, , drop = FALSE]
    }
    selection <- c(selected$metadata, list(regions_table = selected$regions))
    alignment$alignment_stats$pvalue_regions <- selected$metadata
  }
  list(data = data, alignment = alignment, requested_mode = requested_mode,
       effective_mode = effective_mode,
       selection = selection,
       genome_build = if (isTRUE(pre_harmonised)) {
         inherited_alignment$genome_build %||% "GRCh38"
       } else if (effective_mode == "convert" || is.null(reference)) {
         "unknown"
       } else {
         alignment$reference_metadata$build %||% "GRCh38"
       })
}
