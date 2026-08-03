#' Construct canonical CompreSSoR variant keys
#'
#' Canonical keys identify a GRCh38 biallelic SNV without an rsID or external
#' variant dictionary. The allele order is always `REF:ALT`; CompreSSoR stores
#' effects and frequencies for ALT.
#'
#' @param chromosome GRCh38 chromosome (`1`--`22`, `X`, or `Y`).
#' @param position One-based GRCh38 position.
#' @param reference_allele GRCh38 reference allele.
#' @param alternate_allele Alternate allele.
#' @return Character vector in `chromosome:position:REF:ALT` form.
#' @examples
#' compressor_variant_key("1", 12345, "A", "G")
#' @export
compressor_variant_key <- function(chromosome, position,
                                   reference_allele, alternate_allele) {
  lengths <- c(length(chromosome), length(position), length(reference_allele),
               length(alternate_allele))
  n <- max(lengths)
  if (!n || any(lengths != 1L & lengths != n)) {
    stop("key fields must have length one or a common positive length", call. = FALSE)
  }
  chromosome <- toupper(sub("^CHR", "", rep_len(as.character(chromosome), n),
                            ignore.case = TRUE))
  position <- suppressWarnings(as.numeric(rep_len(position, n)))
  reference_allele <- toupper(rep_len(as.character(reference_allele), n))
  alternate_allele <- toupper(rep_len(as.character(alternate_allele), n))
  if (any(is.na(chromosome) | !chromosome %in% c(as.character(1:22), "X", "Y"))) {
    stop("chromosome must be one of 1-22, X, or Y", call. = FALSE)
  }
  if (any(!is.finite(position) | position < 1 | position != floor(position))) {
    stop("position must contain positive whole-number GRCh38 coordinates", call. = FALSE)
  }
  bases <- c("A", "C", "G", "T")
  if (any(is.na(reference_allele) | is.na(alternate_allele) |
          !reference_allele %in% bases | !alternate_allele %in% bases |
          reference_allele == alternate_allele)) {
    stop("REF and ALT must be distinct single A/C/G/T alleles", call. = FALSE)
  }
  paste(chromosome, format(position, scientific = FALSE, trim = TRUE),
        reference_allele, alternate_allele, sep = ":")
}
