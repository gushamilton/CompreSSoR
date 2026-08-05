## Build-specific primary chromosome tables used by the self-contained identity
## key. Coordinates are one-based; offsets are zero-based global-position
## offsets. Keep the chromosome order stable because it is part of the file
## identity contract.

compressor_grch37_chromosome_lengths <- c(
  `1` = 249250621, `2` = 243199373, `3` = 198022430,
  `4` = 191154276, `5` = 180915260, `6` = 171115067,
  `7` = 159138663, `8` = 146364022, `9` = 141213431,
  `10` = 135534747, `11` = 135006516, `12` = 133851895,
  `13` = 115169878, `14` = 107349540, `15` = 102531392,
  `16` = 90354753, `17` = 81195210, `18` = 78077248,
  `19` = 59128983, `20` = 63025520, `21` = 48129895,
  `22` = 51304566, X = 155270560, Y = 59373566
)

compressor_grch38_chromosome_lengths <- c(
  `1` = 248956422, `2` = 242193529, `3` = 198295559,
  `4` = 190214555, `5` = 181538259, `6` = 170805979,
  `7` = 159345973, `8` = 145138636, `9` = 138394717,
  `10` = 133797422, `11` = 135086622, `12` = 133275309,
  `13` = 114364328, `14` = 107043718, `15` = 101991189,
  `16` = 90338345, `17` = 83257441, `18` = 80373285,
  `19` = 58617616, `20` = 64444167, `21` = 46709983,
  `22` = 50818468, X = 156040895, Y = 57227415
)

compressor_supported_builds <- c("GRCh37", "GRCh38")
compressor_identity_table_version <- 1L
compressor_identity_schema <- "compressor_variant_identity_v1"
compressor_identity_base_codes <- c(A = 0L, C = 1L, G = 2L, T = 3L)

compressor_normalize_build <- function(build) {
  if (length(build) != 1L || is.na(build) || !nzchar(trimws(as.character(build)))) {
    stop("build must be one of GRCh37/hg19 or GRCh38/hg38", call. = FALSE)
  }
  key <- gsub("[^A-Z0-9]", "", toupper(trimws(as.character(build))))
  if (key %in% c("GRCH37", "HG19", "37")) return("GRCh37")
  if (key %in% c("GRCH38", "HG38", "38")) return("GRCh38")
  stop("build must be one of GRCh37/hg19 or GRCh38/hg38", call. = FALSE)
}

compressor_chromosome_lengths <- function(build = "GRCh38") {
  build <- compressor_normalize_build(build)
  if (identical(build, "GRCh37")) {
    return(compressor_grch37_chromosome_lengths)
  }
  compressor_grch38_chromosome_lengths
}

compressor_chromosome_offsets <- function(build = "GRCh38") {
  lengths <- compressor_chromosome_lengths(build)
  offsets <- c(0, utils::head(cumsum(as.numeric(lengths)), -1L))
  names(offsets) <- names(lengths)
  offsets
}

## Compatibility constants for code that used the original GRCh38-only table
## directly. New code should call compressor_chromosome_lengths() and
## compressor_chromosome_offsets() with an explicit build.
compressor_grch37_chromosome_offsets <- compressor_chromosome_offsets("GRCh37")
compressor_grch38_chromosome_offsets <- compressor_chromosome_offsets("GRCh38")

parse_canonical_variant_keys <- function(x) {
  x <- as.character(x)
  pieces <- strsplit(x, ":", fixed = TRUE)
  out <- data.frame(
    chromosome = rep(NA_character_, length(x)),
    base_pair_location = rep(NA_integer_, length(x)),
    reference_allele = rep(NA_character_, length(x)),
    alternate_allele = rep(NA_character_, length(x)),
    stringsAsFactors = FALSE
  )
  valid <- lengths(pieces) == 4L
  if (any(valid)) {
    valid_pieces <- pieces[valid]
    out$chromosome[valid] <- vapply(valid_pieces, `[[`, character(1), 1L)
    out$base_pair_location[valid] <- suppressWarnings(as.integer(
      vapply(valid_pieces, `[[`, character(1), 2L)
    ))
    out$reference_allele[valid] <- vapply(valid_pieces, `[[`, character(1), 3L)
    out$alternate_allele[valid] <- vapply(valid_pieces, `[[`, character(1), 4L)
  }
  out
}

compressor_normalize_chromosome <- function(chromosome) {
  chromosome <- toupper(trimws(as.character(chromosome)))
  chromosome <- sub("^CHR", "", chromosome, ignore.case = TRUE)
  chromosome[chromosome == "23"] <- "X"
  chromosome[chromosome == "24"] <- "Y"
  valid <- c(as.character(1:22), "X", "Y")
  if (any(is.na(chromosome) | !chromosome %in% valid)) {
    stop("chromosome must be one of 1-22, X, Y, or their chr/23/24 aliases",
         call. = FALSE)
  }
  chromosome
}

compressor_recycle_identity_fields <- function(chromosome, position,
                                                reference_allele,
                                                alternate_allele) {
  lengths <- c(length(chromosome), length(position), length(reference_allele),
               length(alternate_allele))
  n <- max(lengths)
  if (!n && all(lengths == 0L)) {
    return(list(chromosome = character(), position = numeric(),
                reference_allele = character(), alternate_allele = character()))
  }
  if (!n || any(lengths != 1L & lengths != n)) {
    stop("key fields must have length one or a common positive length", call. = FALSE)
  }
  list(
    chromosome = compressor_normalize_chromosome(rep_len(chromosome, n)),
    position = suppressWarnings(as.numeric(as.character(rep_len(position, n)))),
    reference_allele = toupper(trimws(as.character(rep_len(reference_allele, n)))),
    alternate_allele = toupper(trimws(as.character(rep_len(alternate_allele, n))))
  )
}

compressor_validate_identity_alleles <- function(reference_allele, alternate_allele) {
  bases <- names(compressor_identity_base_codes)
  if (any(is.na(reference_allele) | is.na(alternate_allele) |
          !reference_allele %in% bases | !alternate_allele %in% bases |
          reference_allele == alternate_allele)) {
    stop("REF and ALT must be distinct single A/C/G/T alleles for a biallelic SNV",
         call. = FALSE)
  }
  invisible(TRUE)
}

compressor_validate_identity_positions <- function(chromosome, position, build) {
  if (any(!is.finite(position) | position < 1 | position != floor(position))) {
    stop("position must contain positive whole-number coordinates", call. = FALSE)
  }
  lengths <- compressor_chromosome_lengths(build)
  chromosome_limit <- unname(lengths[chromosome])
  if (any(is.na(chromosome_limit) | position > chromosome_limit)) {
    stop("position lies outside its primary chromosome for the selected build",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Encode build-aware CompreSSoR variant identity
#'
#' The returned global position is zero-based and the substitution code is
#' `4 * REF_code + ALT_code`, with A/C/G/T coded as 0/1/2/3. Both values are
#' lossless only together with the selected genome build.
#'
#' @param chromosome Chromosome 1--22, X, or Y, with optional `chr` prefix or
#'   23/24 aliases for X/Y.
#' @param position One-based coordinate on `build`.
#' @param reference_allele Single-base A/C/G/T REF allele.
#' @param alternate_allele Single-base A/C/G/T ALT allele.
#' @param build Genome build, `GRCh37`/`hg19` or `GRCh38`/`hg38`.
#' @return A list containing canonical fields, `global_position`, and
#'   `substitution`.
compressor_encode_variant_identity <- function(chromosome, position,
                                               reference_allele,
                                               alternate_allele,
                                               build = "GRCh38") {
  build <- compressor_normalize_build(build)
  fields <- compressor_recycle_identity_fields(
    chromosome, position, reference_allele, alternate_allele
  )
  if (!length(fields$chromosome)) {
    return(list(build = build, chromosome = character(), position = numeric(),
                reference_allele = character(), alternate_allele = character(),
                global_position = numeric(), substitution = integer()))
  }
  compressor_validate_identity_positions(fields$chromosome, fields$position, build)
  compressor_validate_identity_alleles(fields$reference_allele,
                                       fields$alternate_allele)
  offsets <- compressor_chromosome_offsets(build)
  global_position <- unname(offsets[fields$chromosome]) + fields$position - 1
  substitution <- unname(compressor_identity_base_codes[fields$reference_allele]) * 4L +
    unname(compressor_identity_base_codes[fields$alternate_allele])
  list(
    build = build,
    chromosome = fields$chromosome,
    position = fields$position,
    reference_allele = fields$reference_allele,
    alternate_allele = fields$alternate_allele,
    global_position = as.numeric(global_position),
    substitution = as.integer(substitution)
  )
}

#' Decode build-aware CompreSSoR variant identity
#'
#' @param global_position Zero-based global position produced by
#'   [compressor_encode_variant_identity()].
#' @param substitution Directed four-bit REF-to-ALT code produced by
#'   [compressor_encode_variant_identity()].
#' @param build Genome build used during encoding.
#' @return A list containing canonical chromosome, position, REF, and ALT
#'   fields, plus the supplied encoded values.
compressor_decode_variant_identity <- function(global_position, substitution,
                                               build = "GRCh38") {
  build <- compressor_normalize_build(build)
  lengths <- c(length(global_position), length(substitution))
  n <- max(lengths)
  if (!n && all(lengths == 0L)) {
    return(list(build = build, chromosome = character(), position = numeric(),
                reference_allele = character(), alternate_allele = character(),
                global_position = numeric(), substitution = integer()))
  }
  if (!n || any(lengths != 1L & lengths != n)) {
    stop("identity fields must have length one or a common positive length",
         call. = FALSE)
  }
  global_position <- suppressWarnings(as.numeric(as.character(
    rep_len(global_position, n)
  )))
  substitution <- suppressWarnings(as.numeric(as.character(
    rep_len(substitution, n)
  )))
  if (any(!is.finite(global_position) | global_position < 0 |
          global_position != floor(global_position))) {
    stop("global_position must contain non-negative whole-number values",
         call. = FALSE)
  }
  if (any(!is.finite(substitution) | substitution < 0 | substitution > 15 |
          substitution != floor(substitution))) {
    stop("substitution must contain integer four-bit REF-to-ALT codes", call. = FALSE)
  }
  substitution <- as.integer(substitution)
  ref_code <- substitution %/% 4L
  alt_code <- substitution %% 4L
  if (any(ref_code == alt_code)) {
    stop("substitution contains a REF=ALT code, not a supported biallelic SNV",
         call. = FALSE)
  }
  lengths_by_chromosome <- compressor_chromosome_lengths(build)
  offsets <- compressor_chromosome_offsets(build)
  genome_length <- sum(as.numeric(lengths_by_chromosome))
  if (any(global_position >= genome_length)) {
    stop("global_position lies outside the selected build's primary chromosomes",
         call. = FALSE)
  }
  chromosome_index <- findInterval(global_position, offsets)
  chromosome <- names(lengths_by_chromosome)[chromosome_index]
  position <- global_position - unname(offsets[chromosome]) + 1
  list(
    build = build,
    chromosome = chromosome,
    position = as.numeric(position),
    reference_allele = names(compressor_identity_base_codes)[ref_code + 1L],
    alternate_allele = names(compressor_identity_base_codes)[alt_code + 1L],
    global_position = global_position,
    substitution = substitution
  )
}

## Short internal aliases make the encoding contract convenient for callers
## that refer to the key as identity rather than variant identity.
compressor_encode_identity <- compressor_encode_variant_identity
compressor_decode_identity <- compressor_decode_variant_identity

#' Define manifest metadata for a build-aware variant identity
#'
#' This helper describes the build-specific identity contract. The identity is
#' self-contained; no reference or chain file is a read-time dependency.
#'
#' @param build Default stored build when `stored_build` is omitted.
#' @param input_build Build of the input rows.
#' @param stored_build Build represented by the encoded identity.
#' @return A manifest-ready list containing the build-specific identity table
#'   and provenance fields.
compressor_identity_manifest <- function(build = "GRCh38", input_build = NULL,
                                         stored_build = NULL) {
  default_build <- compressor_normalize_build(build)
  if (is.null(input_build)) input_build <- default_build
  if (is.null(stored_build)) stored_build <- default_build
  input_build <- compressor_normalize_build(input_build)
  stored_build <- compressor_normalize_build(stored_build)
  if (!identical(input_build, stored_build)) {
    stop("input_build and stored_build must be identical in the compression-core " ,
         "identity manifest; liftover is outside the package", call. = FALSE)
  }
  lengths <- compressor_chromosome_lengths(stored_build)
  offsets <- compressor_chromosome_offsets(stored_build)
  table_id <- paste0(tolower(stored_build), "_primary_1_22_X_Y")
  reference_metadata <- list(id = "none", build = stored_build, status = "not_used",
                             external_reference_required = FALSE)
  chain_metadata <- list(required = FALSE, status = "not_used")
  list(
    schema = compressor_identity_schema,
    encoding = "global_position_plus_directed_ref_alt_substitution",
    position_encoding = "zero_based_global_position",
    substitution_encoding = "uint8_4_times_ref_plus_alt",
    external_reference_required = FALSE,
    input_build = input_build,
    stored_build = stored_build,
    genome_build = stored_build,
    assembly = stored_build,
    assembly_identifier = stored_build,
    chromosome_table = list(
      id = table_id,
      version = compressor_identity_table_version,
      chromosomes = names(lengths),
      lengths = as.list(as.numeric(lengths)),
      offsets = as.list(as.numeric(offsets))
    ),
    chromosome_table_id = table_id,
    chromosome_table_version = compressor_identity_table_version,
    chromosome_lengths = as.list(as.numeric(lengths)),
    chromosome_offsets = as.list(as.numeric(offsets)),
    supported_chromosomes = names(lengths),
    supported_variant_class = "biallelic A/C/G/T SNV",
    reference = reference_metadata,
    chain = chain_metadata,
    effect_allele_is_alt = TRUE,
    other_allele_is_ref = TRUE
  )
}

#' Construct canonical CompreSSoR variant keys
#'
#' Canonical keys identify a biallelic SNV within the selected genome build
#' without an rsID or external variant dictionary. The allele order is always
#' `REF:ALT`; CompreSSoR stores effects and frequencies for ALT. The build is
#' metadata for coordinate validation and is not embedded in the text key.
#'
#' @param chromosome Chromosome 1--22, X, or Y, with optional `chr` prefix or
#'   23/24 aliases for X/Y.
#' @param position One-based position on `build`.
#' @param reference_allele Single-base A/C/G/T reference allele.
#' @param alternate_allele Single-base A/C/G/T alternate allele.
#' @param build Genome build, `GRCh37`/`hg19` or `GRCh38`/`hg38`.
#' @return Character vector in `chromosome:position:REF:ALT` form.
#' @examples
#' compressor_variant_key("1", 12345, "A", "G")
#' compressor_variant_key("chr23", 100, "A", "G", build = "GRCh37")
#' @export
compressor_variant_key <- function(chromosome, position,
                                   reference_allele, alternate_allele,
                                   build = "GRCh38") {
  encoded <- compressor_encode_variant_identity(
    chromosome, position, reference_allele, alternate_allele, build = build
  )
  if (!length(encoded$chromosome)) return(character())
  paste(encoded$chromosome, format(encoded$position, scientific = FALSE,
                                   trim = TRUE), encoded$reference_allele,
        encoded$alternate_allele, sep = ":")
}
