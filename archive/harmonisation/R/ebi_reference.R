# EBI GWAS Catalog harmonisation reference support.
#
utils::globalVariables(c("rsid", "rsid_code", "rsid_other"))

# Chromosomes present in the local EBI/Ensembl release-95 VCF set.
ebi_reference_chromosomes <- function() c(as.character(1:22), "X", "Y", "MT")

# The EBI release-95 files are GRCh38 Ensembl variation resources. They are
# used here as the canonical identity universe; population panels remain
# optional analysis inputs and are not used to define variant_index.

ebi_local_source_paths <- function(source, chromosomes) {
  chromosomes <- as.character(chromosomes)
  if (dir.exists(source)) {
    files <- list.files(source, pattern = "homo_sapiens-chr[^/]+[.]vcf[.]gz$",
                        full.names = TRUE, ignore.case = TRUE)
  } else {
    files <- as.character(source)
  }
  out <- character(length(chromosomes))
  names(out) <- chromosomes
  for (chr in chromosomes) {
    expected <- paste0("homo_sapiens-chr", chr, ".vcf.gz")
    hit <- files[basename(files) == expected]
    if (!length(hit)) {
      hit <- files[grepl(paste0("(^|[/])chr", chr, "[.]vcf([.]gz|[.]bgz)?$"),
                         files, ignore.case = TRUE)]
    }
    if (length(hit)) out[chr] <- normalizePath(hit[1L], mustWork = TRUE)
  }
  out
}

reference_alias_rows <- function(data) {
  if (!nrow(data)) {
    return(data.frame(rsid = character(), variant_id = character(),
                       variant_index = integer(), chromosome = character(),
                       base_pair_location = integer(), stringsAsFactors = FALSE))
  }
  ids <- strsplit(ifelse(is.na(data$rsid), "", as.character(data$rsid)),
                  ";", fixed = TRUE)
  counts <- lengths(ids)
  if (!any(counts)) {
    return(data.frame(rsid = character(), variant_id = character(),
                      variant_index = integer(), chromosome = character(),
                      base_pair_location = integer(), stringsAsFactors = FALSE))
  }
  row <- rep(seq_len(nrow(data)), counts)
  aliases <- data.frame(
    rsid = trimws(unlist(ids, use.names = FALSE)),
    variant_id = data$variant_id[row],
    variant_index = data$variant_index[row],
    chromosome = data$chromosome[row],
    base_pair_location = data$base_pair_location[row],
    stringsAsFactors = FALSE
  )
  aliases <- aliases[!is.na(aliases$rsid) & nzchar(aliases$rsid) &
                       aliases$rsid != ".", , drop = FALSE]
  aliases <- aliases[!duplicated(paste(aliases$rsid, aliases$variant_id, sep = "\r"),
                                 fromLast = FALSE), , drop = FALSE]
  row.names(aliases) <- NULL
  aliases
}

ebi_compact_reference_data <- function(data, previous_position = NA_integer_) {
  positions <- suppressWarnings(as.integer(data$base_pair_location))
  if (!length(positions)) {
    return(data.frame(
      chromosome_code = integer(), position_delta = integer(),
      reference_allele = character(), alternate_allele = character(),
      rsid_code = integer(), rsid_other = character(),
      stringsAsFactors = FALSE
    ))
  }
  position_delta <- positions
  if (!is.na(previous_position)) position_delta[1L] <- positions[1L] - previous_position
  if (length(positions) > 1L) position_delta[-1L] <- diff(positions)
  if (anyNA(position_delta) || any(position_delta < 0L)) {
    stop("EBI VCF positions must be non-decreasing within each chromosome", call. = FALSE)
  }

  rsids <- as.character(data$rsid)
  numeric_candidate <- !is.na(rsids) & grepl("^rs[0-9]+$", rsids, ignore.case = TRUE)
  suffix <- suppressWarnings(as.numeric(sub("^rs", "", rsids,
                                             ignore.case = TRUE)))
  numeric_rsid <- numeric_candidate & !is.na(suffix) &
    suffix <= .Machine$integer.max
  rsid_code <- rep(NA_integer_, length(rsids))
  rsid_code[numeric_rsid] <- as.integer(suffix[numeric_rsid])
  rsid_other <- rsids
  rsid_other[numeric_rsid] <- NA_character_
  data.frame(
    chromosome_code = as.integer(chromosome_sort_key(data$chromosome)),
    position_delta = as.integer(position_delta),
    reference_allele = as.character(data$reference_allele),
    alternate_allele = as.character(data$alternate_allele),
    rsid_code = rsid_code,
    rsid_other = rsid_other,
    stringsAsFactors = FALSE
  )
}

ebi_vcf_stream <- function(path, chunk_lines = 100000L, callback) {
  if (length(chunk_lines) != 1L || !is.numeric(chunk_lines) ||
      is.na(chunk_lines) || chunk_lines < 1L) {
    stop("chunk_lines must be a positive integer", call. = FALSE)
  }
  chunk_lines <- as.integer(chunk_lines)
  if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzip <- Sys.which("gzip")
    if (!nzchar(gzip)) {
      stop("validating a compressed EBI reference requires the 'gzip' executable",
           call. = FALSE)
    }
    status <- suppressWarnings(system2(gzip, c("-t", shQuote(normalizePath(path, mustWork = TRUE))),
                                       stdout = FALSE, stderr = FALSE))
    if (!identical(as.integer(status), 0L)) {
      stop("compressed EBI reference failed gzip integrity validation: ", path,
           call. = FALSE)
    }
  }
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  pending_lines <- character()
  repeat {
    incoming <- readLines(con, n = chunk_lines, warn = FALSE)
    at_eof <- !length(incoming)
    lines <- c(pending_lines, incoming)
    pending_lines <- character()
    lines <- sub("\r$", "", lines)
    lines <- lines[!startsWith(lines, "#")]
    # Never split one genomic locus across callbacks. This makes duplicate
    # canonical identities and all of their rsID aliases deterministic even
    # when a locus falls exactly on a chunk boundary.
    if (!at_eof && length(lines)) {
      locus_fields <- strsplit(lines, "\t", fixed = TRUE)
      locus <- vapply(locus_fields, function(x) {
        if (length(x) >= 2L) paste(x[[1L]], x[[2L]], sep = "\r") else ""
      }, character(1L))
      hold <- locus == utils::tail(locus, 1L)
      pending_lines <- lines[hold]
      lines <- lines[!hold]
    }
    if (!length(lines)) {
      if (at_eof) break
      next
    }
    fields <- strsplit(lines, "\t", fixed = TRUE)
    short <- lengths(fields) < 5L
    if (any(short)) {
      fields[short] <- strsplit(trimws(lines[short]), "[[:space:]]+")
    }
    usable <- lengths(fields) >= 5L
    if (!any(usable)) next
    fields <- fields[usable]
    field <- function(index) {
      vapply(fields, function(x) {
        if (length(x) >= index) x[[index]] else NA_character_
      }, character(1))
    }
    chromosome <- sub("^chr", "", field(1L), ignore.case = TRUE)
    position <- suppressWarnings(as.integer(field(2L)))
    rsid <- field(3L)
    rsid[is.na(rsid) | rsid == "."] <- NA_character_
    reference <- toupper(field(4L))
    alternates <- strsplit(toupper(field(5L)), ",", fixed = TRUE)
    alt_counts <- lengths(alternates)
    rows <- rep(seq_along(alternates), alt_counts)
    if (!length(rows)) next
    data <- data.frame(
      chromosome = chromosome[rows],
      base_pair_location = position[rows],
      reference_allele = reference[rows],
      alternate_allele = unlist(alternates, use.names = FALSE),
      rsid = rsid[rows],
      stringsAsFactors = FALSE
    )
    valid <- !is.na(data$chromosome) & nzchar(data$chromosome) &
      !is.na(data$base_pair_location) & data$base_pair_location >= 1L &
      grepl("^[ACGTN]+$", data$reference_allele) &
      grepl("^[ACGTN]+$", data$alternate_allele) &
      data$reference_allele != data$alternate_allele
    if (!any(valid)) next
    data <- data[valid, , drop = FALSE]
    data$variant_key <- paste(data$chromosome, data$base_pair_location,
                              data$reference_allele, data$alternate_allele,
                              sep = ":")
    if (anyDuplicated(data$variant_key)) {
      groups <- split(seq_len(nrow(data)), factor(
        data$variant_key, levels = unique(data$variant_key)
      ))
      data <- do.call(rbind, lapply(groups, function(index) {
        row <- data[index[1L], , drop = FALSE]
        aliases <- unique(data$rsid[index])
        aliases <- aliases[!is.na(aliases) & nzchar(aliases) & aliases != "."]
        row$rsid <- if (length(aliases)) paste(aliases, collapse = ";") else NA_character_
        row
      }))
    }
    data$variant_id <- data$variant_key
    data$variant_type <- reference_variant_type(
      data$reference_allele, data$alternate_allele
    )
    row.names(data) <- NULL
    callback(data)
    if (at_eof) break
  }
  invisible(NULL)
}

#' Build the complete EBI/Ensembl GRCh38 canonical variant dictionary
#'
#' Processes already-downloaded EBI release-95 chromosome VCFs in bounded
#' chunks and writes one Parquet master file with one global, stable
#' variant_index. The output directory is the single canonical reference
#' asset used by CompreSSoR. Downloading the EBI files is deliberately outside
#' the package.
#'
#' @param output Output directory containing the master dictionary.
#' @param source Local directory/vector of EBI chromosome VCFs.
#' @param chromosomes Chromosomes to include, defaulting to 1-22, X, Y and MT.
#' @param chunk_rows Number of VCF records read per bounded chunk.
#' @param row_group_rows Target number of rows per physical Parquet row group.
#' @param overwrite Whether existing partition files may be replaced.
#' @return A CompreSSoR reference descriptor, invisibly.
#' @export
build_ebi_reference <- function(output, source,
                                chromosomes = ebi_reference_chromosomes(),
                                chunk_rows = 100000L,
                                row_group_rows = 500000L,
                                overwrite = FALSE) {
  require_parquet_backend("build_ebi_reference()")
  if (length(output) != 1L || !is.character(output)) {
    stop("output must be one directory path", call. = FALSE)
  }
  if (file.exists(output) && !dir.exists(output)) {
    stop("output exists and is not a directory: ", output, call. = FALSE)
  }
  if (length(chunk_rows) != 1L || !is.numeric(chunk_rows) ||
      is.na(chunk_rows) || chunk_rows < 1L) {
    stop("chunk_rows must be a positive integer", call. = FALSE)
  }
  chunk_rows <- as.integer(chunk_rows)
  if (length(row_group_rows) != 1L || !is.numeric(row_group_rows) ||
      is.na(row_group_rows) || row_group_rows < 1L) {
    stop("row_group_rows must be a positive integer", call. = FALSE)
  }
  row_group_rows <- as.integer(row_group_rows)
  if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                              no.. = TRUE)) && !isTRUE(overwrite)) {
    stop("output directory is not empty; use overwrite = TRUE", call. = FALSE)
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  master_path <- file.path(output, "variants.parquet")
  chromosomes <- as.character(chromosomes)

  if (missing(source) || is.null(source)) {
    stop("source must be a local directory or vector of downloaded EBI chromosome VCFs",
         call. = FALSE)
  }
  paths <- ebi_local_source_paths(source, chromosomes)
  if (any(!nzchar(paths))) {
    missing <- names(paths)[!nzchar(paths)]
    stop("missing EBI chromosome VCF files: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (isTRUE(overwrite)) {
    unlink(master_path, force = TRUE)
    unlink(list.files(output, pattern = "^chr[^/]+[.]parquet$",
                      full.names = TRUE, ignore.case = TRUE))
    unlink(file.path(output, "aliases"), recursive = TRUE, force = TRUE)
    unlink(file.path(output, "manifest.json"), force = TRUE)
  }

  writer <- NULL
  sink <- NULL
  close_writer <- function() {
    if (!is.null(writer)) try(writer$Close(), silent = TRUE)
    if (!is.null(sink)) try(sink$close(), silent = TRUE)
    writer <<- NULL
    sink <<- NULL
  }
  on.exit(close_writer(), add = TRUE)

  partition_rows <- vector("list", length(paths))
  offset <- 0L
  block_index <- 0L
  for (i in seq_along(paths)) {
    chr <- names(paths)[i]
    message("Building canonical reference partition chr", chr)
    start_offset <- offset
    previous_position <- NA_integer_
    boundary_keys <- character()
    pending <- list()
    pending_rows <- 0L
    pending_first_position <- NA_integer_
    flush_pending <- function() {
      if (!pending_rows) return(invisible(NULL))
      compact <- do.call(rbind, pending)
      compact$position_delta[1L] <- pending_first_position
      compact$position_block <- as.integer(block_index)
      compact <- compact[c("chromosome_code", "position_block", "position_delta",
                           "reference_allele", "alternate_allele", "rsid_code",
                           "rsid_other")]
      table <- arrow::arrow_table(compact)
      if (is.null(writer)) {
        props <- arrow::ParquetWriterProperties$create(
          column_names = names(compact), compression = "zstd",
          compression_level = 9L, use_dictionary = TRUE,
          write_statistics = TRUE
        )
        sink <<- arrow::FileOutputStream$create(master_path)
        writer <<- arrow::ParquetFileWriter$create(table$schema, sink, props)
      }
      writer$WriteTable(table, chunk_size = nrow(compact))
      pending <<- list()
      pending_rows <<- 0L
      pending_first_position <<- NA_integer_
      block_index <<- block_index + 1L
      invisible(NULL)
    }
    process_chunk <- function(data) {
      if (!nrow(data)) return(invisible(NULL))
      if (length(boundary_keys) && !is.na(previous_position)) {
        duplicate_boundary <- data$base_pair_location == previous_position &
          data$variant_key %in% boundary_keys
        if (any(duplicate_boundary)) data <- data[!duplicate_boundary, , drop = FALSE]
      }
      if (!nrow(data)) return(invisible(NULL))
      compact <- ebi_compact_reference_data(data, previous_position = previous_position)
      previous_position <<- utils::tail(data$base_pair_location, 1L)
      boundary_keys <<- data$variant_key[
        data$base_pair_location == previous_position
      ]
      if (!pending_rows) pending_first_position <<- data$base_pair_location[1L]
      pending[[length(pending) + 1L]] <<- compact
      pending_rows <<- pending_rows + nrow(compact)
      offset <<- offset + nrow(data)
      if (pending_rows >= row_group_rows) flush_pending()
      invisible(NULL)
    }
    ebi_vcf_stream(paths[i], chunk_lines = chunk_rows, callback = process_chunk)
    flush_pending()
    partition_rows[[i]] <- list(
      chromosome = chr,
      source = normalizePath(paths[i], mustWork = FALSE),
      source_sha256 = reference_file_digest(paths[i], "sha256"),
      rows = as.integer(offset - start_offset),
      first_index = if (offset > start_offset) start_offset else NA_integer_,
      last_index = if (offset > start_offset) offset - 1L else NA_integer_
    )
  }

  if (is.null(writer)) stop("EBI source produced no valid variants", call. = FALSE)
  close_writer()
  dataset_sha256 <- reference_file_digest(master_path, "sha256")
  metadata <- list(
    id = "ebi_ensembl95_grch38_all_v2",
    build = "GRCh38",
    source = "local EBI/Ensembl GRCh38 VCF files",
    source_release = "Ensembl 95",
    local_path = normalizePath(output, mustWork = FALSE),
    variants_file = "variants.parquet",
    sha256 = dataset_sha256,
    rows = as.integer(offset),
    key = "chrom:pos:ref:alt",
    rsid = "primary_annotation",
    schema_version = "compact_delta_v2",
    storage = list(
      chromosome = "chromosome_code (1-22, X=23, Y=24, MT=25)",
      position = "position_delta, anchored at each physical Parquet row group",
      alleles = "reference_allele and alternate_allele strings",
      rsid = "rsid_code for numeric rs IDs, rsid_other otherwise",
      derived = c("variant_index", "variant_type", "variant_id")
    ),
    snps_only = FALSE,
    variant_types = c("SNV", "MNV", "INDEL"),
    partitioned = TRUE,
    chunk_rows = chunk_rows,
    row_group_rows = row_group_rows,
    partitions = partition_rows,
    created_utc = now_utc()
  )
  write_manifest(metadata, file.path(output, "manifest.json"))
  descriptor <- list(id = metadata$id, build = metadata$build,
                     source = metadata$source, variants = output,
                     cache_dir = dirname(output), metadata = metadata)
  structure(descriptor, class = "compressor_reference_descriptor")
}


reference_alias_table <- function(reference, ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids) & ids != "."]
  empty <- data.frame(rsid = character(), variant_id = character(),
                      variant_index = integer(), chromosome = character(),
                      base_pair_location = integer(), stringsAsFactors = FALSE)
  if (!length(ids)) return(empty)
  resolved <- resolve_reference(reference)
  path <- resolved$variants
  if (is.data.frame(path) && !nrow(path)) return(empty)
  if (is.data.frame(path)) {
    ref <- reference_table(resolved)
    if (is.null(ref) || !nrow(ref)) return(empty)
    aliases <- reference_alias_rows(ref)
    return(aliases[aliases$rsid %in% ids, , drop = FALSE])
  }
  if (is.character(path) && length(path) == 1L && dir.exists(path)) {
    master <- file.path(path, "variants.parquet")
    if (file.exists(master)) {
      out <- tryCatch({
        dataset <- arrow::open_dataset(master, format = "parquet")
        numeric_candidate <- grepl("^rs[0-9]+$", ids, ignore.case = TRUE)
        numeric_suffix <- suppressWarnings(as.numeric(sub(
          "^rs", "", ids, ignore.case = TRUE
        )))
        numeric_ids <- ids[numeric_candidate & !is.na(numeric_suffix) &
                             numeric_suffix <= .Machine$integer.max]
        other_ids <- ids[!(numeric_candidate & !is.na(numeric_suffix) &
                            numeric_suffix <= .Machine$integer.max)]
        chromosome_codes <- integer()
        if (length(numeric_ids)) {
          codes <- as.integer(sub("^rs", "", numeric_ids, ignore.case = TRUE))
          matched <- as.data.frame(dplyr::collect(dplyr::select(
            dplyr::filter(dataset, rsid_code %in% codes), chromosome_code
          )))
          chromosome_codes <- c(chromosome_codes, matched$chromosome_code)
        }
        if (length(other_ids)) {
          exact <- as.data.frame(dplyr::collect(dplyr::select(
            dplyr::filter(dataset, rsid_other %in% other_ids), chromosome_code
          )))
          chromosome_codes <- c(chromosome_codes, exact$chromosome_code)

          # Duplicate canonical records retain all synonyms in one
          # semicolon-delimited annotation. Scan only the compact rsid_other
          # and chromosome-code columns to locate those uncommon aliases, then
          # decode only the chromosomes that contain a requested synonym.
          candidates <- as.data.frame(dplyr::collect(dplyr::select(
            dplyr::filter(dataset, !is.na(rsid_other)),
            chromosome_code, rsid_other
          )))
          if (nrow(candidates)) {
            synonym_hit <- vapply(strsplit(as.character(candidates$rsid_other),
                                           ";", fixed = TRUE),
                                  function(value) any(value %in% other_ids),
                                  logical(1L))
            chromosome_codes <- c(chromosome_codes,
                                  candidates$chromosome_code[synonym_hit])
          }
        }
        chromosome_codes <- unique(chromosome_codes)
        if (!length(chromosome_codes)) return(empty)

        # Delta-coded positions need the preceding rows in each chromosome to
        # reconstruct coordinates. Alias lookup is deliberately the slower,
        # exceptional path; chromosome harmonisation reads and decodes one
        # complete chromosome in one sequential scan.
        decoded <- lapply(chromosome_codes, function(code) {
          reference_table(resolved, chromosome = compact_chromosome_labels(code))
        })
        ref <- do.call(rbind, decoded)
        aliases <- reference_alias_rows(ref)
        aliases[aliases$rsid %in% ids, , drop = FALSE]
      }, error = function(e) empty)
      return(out)
    }
    alias_dir <- resolved$metadata$aliases_dir %||% "aliases"
    if (!grepl("^[/~]", alias_dir)) alias_dir <- file.path(path, alias_dir)
    alias_files <- list.files(alias_dir, pattern = "[.]parquet$",
                              full.names = TRUE, ignore.case = TRUE)
    if (!length(alias_files)) return(empty)
    out <- tryCatch({
      dataset <- arrow::open_dataset(alias_dir, format = "parquet")
      as.data.frame(dplyr::collect(dplyr::filter(dataset, rsid %in% ids)))
    }, error = function(e) {
      parts <- lapply(alias_files, function(file) {
        part <- as.data.frame(arrow::read_parquet(file))
        part[part$rsid %in% ids, , drop = FALSE]
      })
      if (!length(parts)) empty else do.call(rbind, parts)
    })
    if (!nrow(out)) return(empty)
    if (!"chromosome" %in% names(out)) {
      parsed <- parse_canonical_variant_keys(out$variant_id)
      out$chromosome <- parsed$chromosome
      out$base_pair_location <- parsed$base_pair_location
    }
    return(out[!duplicated(paste(out$rsid, out$variant_id, sep = "\r")), ,
               drop = FALSE])
  }
  ref <- reference_table(resolved)
  if (is.null(ref) || !nrow(ref)) return(empty)
  aliases <- reference_alias_rows(ref)
  aliases[aliases$rsid %in% ids, , drop = FALSE]
}

resolve_sumstats_aliases <- function(data, reference, eligible = NULL) {
  if (is.null(reference) || !nrow(data)) return(data)
  if (is.null(eligible)) eligible <- rep(TRUE, nrow(data))
  if (length(eligible) != nrow(data)) stop("eligible alias rows do not match input", call. = FALSE)
  ids <- as.character(data$rsid)
  missing_rsid <- is.na(ids) | !nzchar(ids) | ids == "."
  ids[missing_rsid] <- as.character(data$variant_id[missing_rsid])
  lookup <- reference_alias_table(reference, ids)
  has_id <- !is.na(ids) & nzchar(ids) & ids != "."
  needs_coordinate <- is.na(data$chromosome) | !nzchar(data$chromosome) |
    is.na(data$base_pair_location)
  status <- rep("not_attempted", nrow(data))
  status[has_id & !needs_coordinate] <- "coordinate"
  status[has_id & needs_coordinate & eligible] <- "unresolved"
  status[has_id & needs_coordinate & !eligible] <- "liftover_failed"
  if (!nrow(lookup)) {
    data$.compressor_reference_alias_status <- status
    return(data)
  }
  for (i in which(status == "unresolved")) {
    candidates <- lookup[lookup$rsid == ids[i], , drop = FALSE]
    variants <- unique(candidates$variant_id)
    if (length(variants) != 1L) {
      status[i] <- if (length(variants)) "ambiguous" else "unresolved"
      next
    }
    hit <- candidates[match(variants[1L], candidates$variant_id), , drop = FALSE]
    if (is.na(data$chromosome[i]) || !nzchar(data$chromosome[i])) {
      data$chromosome[i] <- hit$chromosome
    }
    if (is.na(data$base_pair_location[i])) {
      data$base_pair_location[i] <- hit$base_pair_location
    }
    status[i] <- "resolved"
  }
  data$.compressor_reference_alias_status <- status
  data
}
