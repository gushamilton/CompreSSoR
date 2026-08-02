normalise_reference_columns <- function(data) {
  data <- clean_input_names(data)
  alias_map <- list(
    chromosome = c("chromosome", "chr", "CHR", "#chrom", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "bp"),
    effect_allele = c("effect_allele", "ea", "EA", "a1", "A1", "alt", "ALT"),
    other_allele = c("other_allele", "oa", "NEA", "nea", "a2", "A2", "ref", "REF"),
    variant_id = c("variant_id", "SNPID", "SNP_ID", "SNP", "snp", "variant", "rsid", "rsID", "ID"),
    effect_allele_frequency = c("effect_allele_frequency", "eaf", "EAF", "af", "AF", "effect_af")
  )
  for (target in names(alias_map)) data <- rename_first_alias(data, target, alias_map[[target]])
  required <- c("chromosome", "base_pair_location", "effect_allele", "other_allele")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Missing required reference columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data$chromosome <- sub("^chr", "", as.character(data$chromosome), ignore.case = TRUE)
  data$base_pair_location <- as.integer(data$base_pair_location)
  data$effect_allele <- toupper(as.character(data$effect_allele))
  data$other_allele <- toupper(as.character(data$other_allele))
  if (!"variant_id" %in% names(data)) {
    data$variant_id <- paste(data$chromosome, data$base_pair_location,
                             data$other_allele, data$effect_allele, sep = "_")
  }
  data$variant_id <- as.character(data$variant_id)
  missing_id <- is.na(data$variant_id) | !nzchar(data$variant_id) | data$variant_id == "."
  if (any(missing_id)) {
    data$variant_id[missing_id] <- paste(data$chromosome[missing_id],
                                         data$base_pair_location[missing_id],
                                         data$other_allele[missing_id],
                                         data$effect_allele[missing_id], sep = "_")
  }
  if (!"effect_allele_frequency" %in% names(data)) data$effect_allele_frequency <- NA_real_
  data$effect_allele_frequency <- as.numeric(data$effect_allele_frequency)
  bad_eaf <- !is.na(data$effect_allele_frequency) &
    (data$effect_allele_frequency < 0 | data$effect_allele_frequency > 1)
  if (any(bad_eaf)) stop("reference effect_allele_frequency must be between 0 and 1", call. = FALSE)
  data
}

complement_allele <- function(x) {
  x <- toupper(as.character(x))
  valid <- !is.na(x) & nchar(x) == 1L & x %in% c("A", "T", "C", "G")
  out <- rep(NA_character_, length(x))
  out[valid] <- chartr("ATCG", "TAGC", x[valid])
  out
}

reference_normalized_cache_path <- function(resolved) {
  if (is.data.frame(resolved$variants)) return(NULL)
  path <- resolved$variants
  if (!is.character(path) || length(path) != 1L || grepl("[.]parquet$", path, ignore.case = TRUE)) {
    return(NULL)
  }
  sha256 <- resolved$metadata$sha256 %||% NULL
  if (is.null(sha256) || length(sha256) != 1L || !nzchar(sha256)) return(NULL)
  file.path(reference_cache_dir(resolved$cache_dir %||% NULL),
            paste0("normalized-", tolower(sha256), ".parquet"))
}

write_normalized_reference_cache <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("compressor-normalized-reference-", tmpdir = dirname(path), fileext = ".parquet")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  tryCatch({
    arrow::write_parquet(data, temporary, compression = "zstd", compression_level = 7,
                         write_statistics = TRUE, use_dictionary = TRUE,
                         chunk_size = 65536L)
    if (file.exists(path)) unlink(path, force = TRUE)
    file.rename(temporary, path)
  }, error = function(e) NULL)
  invisible(path)
}

reference_table <- function(reference) {
  if (is.null(reference)) return(NULL)
  resolved <- resolve_reference(reference)
  if (is.data.frame(resolved$variants)) {
    out <- normalise_reference_columns(resolved$variants)
  } else {
    path <- resolved$variants
    if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
      stop("resolved reference does not contain a readable variants asset", call. = FALSE)
    }
    if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
      out <- normalise_reference_columns(arrow::read_parquet(path))
    } else if (grepl("[.]bim(?:[.]gz)?$", path, ignore.case = TRUE)) {
      out <- data.table::fread(
        path, data.table = FALSE, header = FALSE,
        col.names = c("chromosome", "variant_id", "genetic_distance",
                      "base_pair_location", "other_allele", "effect_allele"),
        showProgress = FALSE
      )
      out <- normalise_reference_columns(out)
    } else {
      normalized_cache <- reference_normalized_cache_path(resolved)
      if (!is.null(normalized_cache) && file.exists(normalized_cache)) {
        out <- normalise_reference_columns(arrow::read_parquet(normalized_cache))
      } else {
        if (grepl("[.]gz$", path, ignore.case = TRUE)) {
          command <- paste("gzip -dc", shQuote(normalizePath(path)))
          raw_reference <- data.table::fread(cmd = command, data.table = FALSE, showProgress = FALSE)
        } else {
          raw_reference <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
        }
        out <- normalise_reference_columns(raw_reference)
        if (!is.null(normalized_cache)) write_normalized_reference_cache(out, normalized_cache)
      }
    }
  }
  attr(out, "reference_metadata") <- resolved$metadata %||% list(
    id = resolved$id %||% "local",
    build = resolved$build %||% "GRCh38"
  )
  metadata <- attr(out, "reference_metadata")
  normalized_cache <- reference_normalized_cache_path(resolved)
  if (!is.null(normalized_cache) && file.exists(normalized_cache)) {
    metadata$normalized_cache_path <- normalizePath(normalized_cache, mustWork = FALSE)
    attr(out, "reference_metadata") <- metadata
  }
  if (!is.null(metadata$sha256)) attr(out, "reference_sha256") <- metadata$sha256
  out
}

reference_hash <- function(reference_data) {
  if (is.null(reference_data)) return(NA_character_)
  known_hash <- attr(reference_data, "reference_sha256")
  if (!is.null(known_hash) && length(known_hash) == 1L && nzchar(known_hash)) return(known_hash)
  digest::digest(reference_data, algo = "sha256", serialize = TRUE)
}

harmonise_to_reference <- function(data, reference, strict = FALSE, preserve = TRUE,
                                   reference_data = NULL, reference_key = NULL) {
  input_rows <- nrow(data)
  ref <- if (!is.null(reference_data)) reference_data else reference_table(reference)
  if (is.null(ref)) {
    data$harmonisation_status <- "unreferenced"
    data$harmonisation_flip <- FALSE
    data$input_duplicate <- duplicated(data$variant_id)
    return(list(data = data, reference_hash = NA_character_,
                reference_rows = NA_integer_, reference_metadata = NULL,
                alignment_stats = list(
                  input_rows = input_rows, output_rows = nrow(data),
                  aligned_rows = 0L, duplicate_input_rows = sum(data$input_duplicate),
                  unmatched_rows = 0L, incompatible_rows = 0L, ambiguous_rows = 0L,
                  dropped_duplicates = 0L, dropped_unmatched = 0L,
                  dropped_incompatible = 0L, dropped_ambiguous = 0L
                )))
  }
  if (anyNA(ref$chromosome) || anyNA(ref$base_pair_location) ||
      any(ref$base_pair_location < 1L, na.rm = TRUE)) {
    stop("reference requires chromosome and positive base-pair locations", call. = FALSE)
  }
  if (anyNA(ref$effect_allele) || anyNA(ref$other_allele) ||
      any(ref$effect_allele == "" | ref$other_allele == "", na.rm = TRUE)) {
    stop("reference alleles must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(ref$variant_id)) stop("reference contains duplicated variant_id values", call. = FALSE)

  input_duplicate <- duplicated(data$variant_id)
  if (isTRUE(strict) && any(input_duplicate)) {
    stop("input contains duplicated variant_id values", call. = FALSE)
  }
  key <- match(data$variant_id, ref$variant_id)
  fallback <- is.na(key)
  if (any(fallback)) {
    data_key <- paste(data$chromosome, data$base_pair_location, sep = ":")
    ref_key <- reference_key %||% paste(ref$chromosome, ref$base_pair_location, sep = ":")
    key[fallback] <- match(data_key[fallback], ref_key)
  }
  unmatched <- is.na(key)
  ref_effect <- ref$effect_allele[key]
  ref_other <- ref$other_allele[key]
  input_effect <- data$effect_allele
  input_other <- data$other_allele
  input_effect_complement <- complement_allele(input_effect)
  input_other_complement <- complement_allele(input_other)
  direct <- !unmatched & !is.na(input_effect) & !is.na(input_other) &
    input_effect == ref_effect & input_other == ref_other
  reverse <- !unmatched & !is.na(input_effect) & !is.na(input_other) &
    input_effect == ref_other & input_other == ref_effect
  complement_direct <- !unmatched & !is.na(input_effect_complement) & !is.na(input_other_complement) &
    input_effect_complement == ref_effect & input_other_complement == ref_other
  complement_reverse <- !unmatched & !is.na(input_effect_complement) & !is.na(input_other_complement) &
    input_effect_complement == ref_other & input_other_complement == ref_effect
  candidates <- cbind(direct, reverse, complement_direct, complement_reverse)
  candidate_count <- rowSums(candidates)
  incompatible <- !unmatched & candidate_count == 0L
  multi <- !unmatched & candidate_count > 1L
  input_eaf <- data$effect_allele_frequency
  reference_eaf <- ref$effect_allele_frequency[key]
  same_error <- abs(input_eaf - reference_eaf)
  flip_error <- abs((1 - input_eaf) - reference_eaf)
  informative <- multi & is.finite(same_error) & is.finite(flip_error)
  enough_separation <- informative & abs(same_error - flip_error) > 0.1
  ambiguous <- multi & !enough_separation
  bad <- unmatched | incompatible | ambiguous
  if (isTRUE(strict) && any(bad)) {
    stop(
      "reference alignment failed: unmatched=", sum(unmatched),
      ", incompatible=", sum(incompatible), ", ambiguous=", sum(ambiguous),
      call. = FALSE
    )
  }
  status <- rep("aligned", input_rows)
  status[unmatched] <- "unmatched"
  status[incompatible] <- "incompatible"
  status[ambiguous] <- "ambiguous"
  flip_flag <- reverse | complement_reverse
  if (any(multi & !ambiguous)) flip_flag[multi & !ambiguous] <- flip_error[multi & !ambiguous] < same_error[multi & !ambiguous]
  flip_flag[bad] <- FALSE

  aligned <- !bad
  if (any(flip_flag)) {
    data$beta[flip_flag] <- -data$beta[flip_flag]
    data$effect_allele_frequency[flip_flag] <- 1 - data$effect_allele_frequency[flip_flag]
  }
  data$chromosome[aligned] <- ref$chromosome[key[aligned]]
  data$base_pair_location[aligned] <- ref$base_pair_location[key[aligned]]
  data$effect_allele[aligned] <- ref$effect_allele[key[aligned]]
  data$other_allele[aligned] <- ref$other_allele[key[aligned]]
  data$variant_id[aligned] <- ref$variant_id[key[aligned]]
  data$harmonisation_status <- status
  data$harmonisation_flip <- flip_flag
  data$input_duplicate <- input_duplicate

  dropped_unmatched <- if (isTRUE(preserve)) 0L else sum(unmatched)
  dropped_incompatible <- if (isTRUE(preserve)) 0L else sum(incompatible)
  dropped_ambiguous <- if (isTRUE(preserve)) 0L else sum(ambiguous)
  if (!isTRUE(preserve)) data <- data[aligned, , drop = FALSE]
  list(data = data, reference_hash = reference_hash(ref), reference_rows = nrow(ref),
       reference_metadata = attr(ref, "reference_metadata"),
       alignment_stats = list(
         input_rows = input_rows,
         output_rows = nrow(data),
         aligned_rows = sum(aligned),
         duplicate_input_rows = sum(input_duplicate),
         unmatched_rows = sum(unmatched),
         incompatible_rows = sum(incompatible),
         ambiguous_rows = sum(ambiguous),
         dropped_duplicates = 0L,
         dropped_unmatched = dropped_unmatched,
         dropped_incompatible = dropped_incompatible,
         dropped_ambiguous = dropped_ambiguous
       ))
}

harmonise_by_chromosome <- function(data, reference, strict = FALSE, preserve = TRUE,
                                    chrom_threads = 1L) {
  ref <- reference_table(reference)
  if (is.null(ref) || nrow(data) < 2L || as.integer(chrom_threads) <= 1L) {
    return(harmonise_to_reference(data, reference, strict = strict, preserve = preserve,
                                  reference_data = ref))
  }
  chrom_threads <- max(1L, as.integer(chrom_threads))
  order_column <- ".compressor_input_order"
  if (order_column %in% names(data)) stop("reserved column name is present: ", order_column, call. = FALSE)
  data[[order_column]] <- seq_len(nrow(data))
  groups <- split(seq_len(nrow(data)), as.character(data$chromosome), drop = TRUE)
  ref_key <- paste(ref$chromosome, ref$base_pair_location, sep = ":")
  worker <- function(idx) harmonise_to_reference(
    data[idx, , drop = FALSE], NULL, strict = strict, preserve = preserve,
    reference_data = ref, reference_key = ref_key
  )
  if (.Platform$OS.type == "windows") {
    cluster <- parallel::makeCluster(min(chrom_threads, length(groups)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parts <- parallel::parLapply(cluster, groups, worker)
  } else {
    parts <- parallel::mclapply(groups, worker, mc.cores = min(chrom_threads, length(groups)),
                                mc.preschedule = FALSE)
  }
  combined <- do.call(rbind, lapply(parts, `[[`, "data"))
  combined <- combined[order(combined[[order_column]]), , drop = FALSE]
  combined[[order_column]] <- NULL
  row.names(combined) <- NULL
  stats_names <- names(parts[[1L]]$alignment_stats)
  stats <- lapply(stats_names, function(nm) sum(vapply(parts, function(x) as.numeric(x$alignment_stats[[nm]]), numeric(1))))
  names(stats) <- stats_names
  stats$input_rows <- nrow(data)
  stats$output_rows <- nrow(combined)
  list(data = combined, reference_hash = parts[[1L]]$reference_hash,
       reference_rows = parts[[1L]]$reference_rows,
       reference_metadata = parts[[1L]]$reference_metadata,
       alignment_stats = stats)
}
