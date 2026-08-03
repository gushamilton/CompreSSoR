normalise_variant_set_columns <- function(data) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  alias_map <- list(
    variant_id = c("variant_id", "variant_id_grch38", "SNPID", "SNP", "snp", "ID", "id"),
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
  data$variant_id <- as.character(data$variant_id)
  data$rsid <- as.character(data$rsid)
  data$chromosome <- sub("^chr", "", as.character(data$chromosome), ignore.case = TRUE)
  data$base_pair_location <- suppressWarnings(as.integer(data$base_pair_location))
  data
}

read_variant_set <- function(variant_set) {
  if (is.data.frame(variant_set)) {
    out <- normalise_variant_set_columns(variant_set)
    metadata <- list(id = "in_memory", rows = nrow(out))
  } else {
    if (length(variant_set) != 1L || !is.character(variant_set) || !file.exists(variant_set)) {
      stop("variant_set must be a data.frame or an existing panel file", call. = FALSE)
    }
    if (grepl("[.]bim(?:[.]gz)?$", variant_set, ignore.case = TRUE)) {
      bim_args <- list(data.table = FALSE, header = FALSE,
                       col.names = c("chromosome", "variant_id", "genetic_distance",
                                     "base_pair_location", "other_allele", "effect_allele"),
                       showProgress = FALSE)
      if (grepl("[.]gz$", variant_set, ignore.case = TRUE)) {
        bim_args$cmd <- paste("gzip -dc", shQuote(normalizePath(variant_set)))
        out <- do.call(data.table::fread, bim_args)
      } else {
        bim_args$input <- variant_set
        out <- do.call(data.table::fread, bim_args)
      }
    } else if (grepl("[.]parquet$", variant_set, ignore.case = TRUE)) {
      out <- arrow::read_parquet(variant_set)
    } else {
      out <- data.table::fread(variant_set, data.table = FALSE, showProgress = FALSE)
    }
    out <- normalise_variant_set_columns(out)
    metadata <- list(
      id = "local",
      local_path = normalizePath(variant_set, mustWork = FALSE),
      sha256 = reference_file_digest(variant_set, "sha256"),
      rows = nrow(out)
    )
  }
  attr(out, "variant_set_metadata") <- metadata
  out
}

normalise_variant_key <- function(x) {
  x <- as.character(x)
  sub("^chr", "", x, ignore.case = TRUE)
}

variant_set_membership <- function(data, panel) {
  panel_id <- unique(normalise_variant_key(panel$variant_id[!is.na(panel$variant_id) & nzchar(panel$variant_id)]))
  panel_rsid <- unique(panel$rsid[!is.na(panel$rsid) & nzchar(panel$rsid)])
  panel_rsid <- unique(c(panel_rsid, panel_id[grepl("^rs", panel_id, ignore.case = TRUE)]))
  keep <- normalise_variant_key(data$variant_id) %in% panel_id
  if ("rsid" %in% names(data)) keep <- keep | as.character(data$rsid) %in% panel_rsid
  panel_coord <- paste(panel$chromosome, panel$base_pair_location, sep = ":")
  data_coord <- paste(data$chromosome, data$base_pair_location, sep = ":")
  keep | (!is.na(data$base_pair_location) & data_coord %in% panel_coord)
}

prepare_sumstats_data <- function(raw, reference, mode = c("qc", "convert", "all", "core", "hm3"),
                                  variant_set = NULL, strict = FALSE, chrom_threads = 1L,
                                  drop_unresolved = TRUE, input_build = "GRCh38", chain = NULL) {
  requested_mode <- match.arg(mode)
  effective_mode <- if (requested_mode == "all") "qc" else requested_mode
  if (effective_mode %in% c("core", "hm3") && is.null(variant_set)) {
    stop("mode='", requested_mode, "' requires variant_set = a core/HM3 panel file or data.frame", call. = FALSE)
  }
  alignment_reference <- if (effective_mode == "convert") NULL else reference
  if (effective_mode != "convert" && is.null(alignment_reference)) {
    stop("reference is required for harmonised CompreSSoR conversion; use mode='convert' explicitly to bypass it",
         call. = FALSE)
  }
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
  data <- alignment$data
  alias_status <- raw$.compressor_reference_alias_status %||% rep("not_attempted", nrow(raw))
  liftover_status <- raw$.compressor_liftover_status %||% rep("not_needed", nrow(raw))
  alignment$alignment_stats$alias_resolution <- as.list(table(alias_status, useNA = "ifany"))
  alignment$alignment_stats$liftover <- as.list(table(liftover_status, useNA = "ifany"))
  filter_stats <- list(name = "all", input_rows = nrow(data), kept_rows = nrow(data), dropped_rows = 0L)
  if (effective_mode %in% c("core", "hm3")) {
    panel <- read_variant_set(variant_set)
    keep <- variant_set_membership(data, panel)
    filter_stats <- c(
      list(name = effective_mode),
      attr(panel, "variant_set_metadata") %||% list(),
      list(input_rows = nrow(data), kept_rows = sum(keep), dropped_rows = sum(!keep))
    )
    data <- data[keep, , drop = FALSE]
    if (!nrow(data)) stop("variant_set filter retained no variants", call. = FALSE)
  }
  alignment$alignment_stats$variant_set <- filter_stats
  list(data = data, alignment = alignment, requested_mode = requested_mode,
       effective_mode = effective_mode,
       genome_build = if (effective_mode == "convert" || is.null(reference)) {
         "unknown"
       } else {
         alignment$reference_metadata$build %||% "GRCh38"
       })
}
