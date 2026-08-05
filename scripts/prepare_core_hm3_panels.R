#!/usr/bin/env Rscript

# Stage the frozen GRCh38 core and HM3 panel artifacts used by the benchmarks.
# The source files remain outside the synced project; this script records their
# hashes and writes compact, canonical identity dictionaries to an external
# output directory.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("--help", "-h"))) {
  cat(paste(
    "usage: prepare_core_hm3_panels.R [--core-source PATH] [--hm3-source PATH]",
    "[--output-dir PATH]", "\n",
    "Defaults use the frozen MR-atlas panel files on the Mac mini and write",
    "prepared panels outside the synced project.\n", sep = " "))
  quit(save = "no", status = 0L)
}

parse_args <- function(values) {
  result <- list(
    core_source = Sys.getenv(
      "COMPRESSOR_CORE_SOURCE",
      "/Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz"
    ),
    hm3_source = Sys.getenv(
      "COMPRESSOR_HM3_SOURCE",
      "/Volumes/crucial_x9/mr_atlas/data/panels/hm3/hm3_grch38_canonical_from_bfile.tsv.gz"
    ),
    output_dir = Sys.getenv(
      "COMPRESSOR_PANEL_OUTPUT",
      "/Volumes/crucial_x9/CompreSSoR-benchmarks/prepared-core-hm3-panels"
    )
  )
  if (!length(values)) return(result)
  if (length(values) %% 2L) stop("options must be supplied as --name value pairs", call. = FALSE)
  for (i in seq(1L, length(values), by = 2L)) {
    key <- sub("^--", "", values[[i]])
    key <- sub("-", "_", key, fixed = TRUE)
    if (!key %in% names(result)) stop("unknown option: --", key, call. = FALSE)
    result[[key]] <- values[[i + 1L]]
  }
  result
}

options <- parse_args(args)
`%||%` <- function(x, y) if (is.null(x)) y else x
for (path in options[c("core_source", "hm3_source")]) {
  if (!file.exists(path)) stop("panel source not found: ", path, call. = FALSE)
}
dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)

shell_quote <- function(x) paste0("'", gsub("'", "'\\''", x, fixed = TRUE), "'")

read_panel_columns <- function(path, aliases) {
  command <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    paste("gzip -dc", shell_quote(normalizePath(path)))
  } else {
    NULL
  }
  header_args <- list(data.table = FALSE, nrows = 0L, showProgress = FALSE,
                      check.names = FALSE)
  if (is.null(command)) header_args$input <- path else header_args$cmd <- command
  header <- do.call(fread, header_args)
  names(header) <- sub("^\\ufeff", "", trimws(names(header)))
  selected <- unname(vapply(aliases, function(options) {
    hit <- intersect(options, names(header))
    if (length(hit)) hit[[1L]] else NA_character_
  }, character(1L)))
  selected <- selected[!is.na(selected)]
  if (!length(selected)) {
    stop("panel has no recognizable identity columns: ", path, call. = FALSE)
  }
  read_args <- list(data.table = FALSE, showProgress = FALSE,
                    check.names = FALSE, select = selected)
  if (is.null(command)) read_args$input <- path else read_args$cmd <- command
  out <- do.call(fread, read_args)
  names(out) <- sub("^\\ufeff", "", trimws(names(out)))
  out
}

rename_first <- function(data, target, aliases) {
  hit <- intersect(aliases, names(data))
  if (length(hit) && !identical(hit[[1L]], target)) names(data)[match(hit[[1L]], names(data))] <- target
  data
}

chromosome_lengths <- c(
  `1` = 248956422, `2` = 242193529, `3` = 198295559,
  `4` = 190214555, `5` = 181538259, `6` = 170805979,
  `7` = 159345973, `8` = 145138636, `9` = 138394717,
  `10` = 133797422, `11` = 135086622, `12` = 133275309,
  `13` = 114364328, `14` = 107043718, `15` = 101991189,
  `16` = 90338345, `17` = 83257441, `18` = 80373285,
  `19` = 58617616, `20` = 64444167, `21` = 46709983,
  `22` = 50818468, X = 156040895, Y = 57227415
)

canonicalise_panel <- function(data, source, kind) {
  alias_map <- list(
    variant_id = c("variant_id", "variant_id_grch38", "vid", "SNPID", "SNP", "snp", "ID", "id"),
    rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID"),
    chromosome = c("chromosome", "chr", "CHR", "#CHROM", "CHROM", "chrom"),
    base_pair_location = c("base_pair_location", "position", "pos", "POS", "BP", "bp"),
    other_allele = c("other_allele", "NEA", "REF", "ref", "A2"),
    effect_allele = c("effect_allele", "EA", "ALT", "alt", "A1")
  )
  for (target in names(alias_map)) data <- rename_first(data, target, alias_map[[target]])
  if (!"variant_id" %in% names(data)) data$variant_id <- NA_character_
  if (!"rsid" %in% names(data)) data$rsid <- NA_character_
  if (!"chromosome" %in% names(data)) data$chromosome <- NA_character_
  if (!"base_pair_location" %in% names(data)) data$base_pair_location <- NA_character_
  if (!"other_allele" %in% names(data)) data$other_allele <- NA_character_
  if (!"effect_allele" %in% names(data)) data$effect_allele <- NA_character_

  chromosome <- toupper(sub("^chr", "", trimws(as.character(data$chromosome)), ignore.case = TRUE))
  chromosome[chromosome == "23"] <- "X"
  chromosome[chromosome == "24"] <- "Y"
  position <- suppressWarnings(as.numeric(data$base_pair_location))
  ref <- toupper(trimws(as.character(data$other_allele)))
  alt <- toupper(trimws(as.character(data$effect_allele)))
  supplied_id <- toupper(sub("^chr", "", trimws(as.character(data$variant_id)), ignore.case = TRUE))
  supplied_id[supplied_id %in% c("", ".", "NA")] <- NA_character_
  can_build_id <- is.na(supplied_id) & !is.na(chromosome) & is.finite(position) &
    !is.na(ref) & !is.na(alt)
  supplied_id[can_build_id] <- paste(chromosome[can_build_id], position[can_build_id],
                                     ref[can_build_id], alt[can_build_id], sep = ":")
  valid <- chromosome %in% names(chromosome_lengths) & is.finite(position) &
    position == floor(position) & position >= 1 &
    position <= unname(chromosome_lengths[chromosome]) &
    grepl("^[ACGT]$", ref) & grepl("^[ACGT]$", alt) & ref != alt &
    grepl("^(?:[1-9]|1[0-9]|2[0-2]|X|Y):[0-9]+:[ACGT]:[ACGT]$", supplied_id)
  invalid_rows <- which(!valid)
  if (length(invalid_rows)) {
    message(kind, ": excluding ", length(invalid_rows),
            " source row(s) outside the canonical biallelic-SNV contract")
  }
  out <- data.frame(
    variant_id = supplied_id[valid],
    chromosome = chromosome[valid],
    base_pair_location = as.integer(position[valid]),
    other_allele = ref[valid],
    effect_allele = alt[valid],
    rsid = as.character(data$rsid)[valid],
    stringsAsFactors = FALSE
  )
  out$rsid[trimws(out$rsid) %in% c("", ".", "NA")] <- NA_character_
  out <- out[!duplicated(out$variant_id), , drop = FALSE]
  row.names(out) <- NULL
  attr(out, "invalid_rows") <- as.integer(length(invalid_rows))
  out
}

write_panel <- function(data, path) {
  temporary <- tempfile("compressor-panel-", fileext = ".tsv")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  # Some data.table builds lack zlib headers even though the gzip executable
  # is present. Keep compression deterministic without making that optional
  # build detail a prerequisite for panel reproduction.
  fwrite(data, temporary, sep = "\t", quote = FALSE, na = "")
  status <- system2("gzip", c("-n", "-c", temporary), stdout = path)
  if (!identical(status, 0L)) stop("could not gzip prepared panel: ", path, call. = FALSE)
  path
}

write_chromosome_shards <- function(data, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  chromosome_order <- c(as.character(1:22), "X", "Y")
  files <- list()
  for (chromosome in chromosome_order) {
    keep <- !is.na(data$chromosome) & data$chromosome == chromosome
    if (!any(keep)) next
    path <- file.path(output_dir, paste0("chr", chromosome, ".tsv.gz"))
    shard <- data[keep, , drop = FALSE]
    write_panel(shard, path)
    files[[chromosome]] <- list(
      file = basename(path),
      rows = nrow(shard),
      sha256 = digest(path, algo = "sha256", file = TRUE)
    )
  }
  files
}

panel_specs <- list(
  core = list(
    source = options$core_source,
    aliases = list(
      variant_id = c("variant_id", "variant_id_grch38", "vid", "SNPID", "SNP", "snp", "ID", "id"),
      rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID"),
      chromosome = c("chromosome", "chr", "CHR", "#CHROM", "CHROM", "chrom"),
      base_pair_location = c("base_pair_location", "position", "pos", "POS", "BP", "bp"),
      other_allele = c("other_allele", "NEA", "REF", "ref", "A2"),
      effect_allele = c("effect_allele", "EA", "ALT", "alt", "A1")
    )
  ),
  hm3 = list(
    source = options$hm3_source,
    aliases = list(
      variant_id = c("variant_id", "variant_id_grch38", "vid", "SNPID", "SNP", "snp", "ID", "id"),
      rsid = c("rsid", "rsID", "RSID", "rs_id", "RS_ID"),
      chromosome = c("chromosome", "chr", "CHR", "#CHROM", "CHROM", "chrom"),
      base_pair_location = c("base_pair_location", "position", "pos", "POS", "BP", "bp"),
      other_allele = c("other_allele", "NEA", "REF", "ref", "A2"),
      effect_allele = c("effect_allele", "EA", "ALT", "alt", "A1")
    )
  )
)

manifest <- list(
  format = "CompreSSoR prepared variant panels",
  format_version = "1",
  build = "GRCh38",
  identity = "chromosome:position:REF:ALT; primary chromosomes 1-22, X and Y; biallelic SNVs",
  transform = c(
    "read only recognizable identity columns from the source TSV/TSV.GZ",
    "normalise chromosome labels and REF/ALT case",
    "validate canonical GRCh38 SNV identity",
    "deduplicate canonical variant_id, retaining first source occurrence"
  ),
  panels = list()
)

for (panel_name in names(panel_specs)) {
  spec <- panel_specs[[panel_name]]
  source <- normalizePath(spec$source, mustWork = TRUE)
  source_data <- read_panel_columns(source, spec$aliases)
  panel <- canonicalise_panel(source_data, source, panel_name)
  output <- file.path(options$output_dir, paste0(panel_name, ".tsv.gz"))
  write_panel(panel, output)
  shard_dir <- file.path(options$output_dir, paste0(panel_name, "_by_chrom"))
  shards <- write_chromosome_shards(panel, shard_dir)
  manifest$panels[[panel_name]] <- list(
    source = source,
    source_sha256 = digest(source, algo = "sha256", file = TRUE),
    source_rows = nrow(source_data),
    excluded_invalid_rows = attr(panel, "invalid_rows") %||% 0L,
    prepared = normalizePath(output, mustWork = FALSE),
    prepared_sha256 = digest(output, algo = "sha256", file = TRUE),
    prepared_rows = nrow(panel),
    columns_read = names(source_data),
    chromosome_shard_directory = normalizePath(shard_dir, mustWork = FALSE),
    chromosome_shards = shards
  )
}

manifest$created_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
manifest_path <- file.path(options$output_dir, "panel_manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE)
cat("Prepared panels in ", options$output_dir, "\n", sep = "")
cat("Manifest: ", manifest_path, "\n", sep = "")
