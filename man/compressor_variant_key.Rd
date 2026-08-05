\name{compressor_variant_key}
\alias{compressor_variant_key}
\title{Construct Build-Aware Canonical CompreSSoR Variant Keys}
\description{Constructs self-contained biallelic-SNV keys in
\code{chromosome:position:REF:ALT} form for GRCh37 or GRCh38. Chromosome
aliases including a \code{chr} prefix and 23/24 for X/Y are normalized;
coordinates are validated against the selected primary chromosome table.
CompreSSoR stores effects and frequencies for ALT. The build is metadata for
validation and is not embedded in the text key.}
\usage{
compressor_variant_key(
  chromosome, position, reference_allele, alternate_allele,
  build = "GRCh38"
)
}
\arguments{
  \item{chromosome}{Chromosome 1--22, X, or Y, with an optional \code{chr}
    prefix or 23/24 aliases for X/Y.}
  \item{position}{One-based position on \code{build}.}
  \item{reference_allele}{Single-base A/C/G/T reference allele.}
  \item{alternate_allele}{Single-base A/C/G/T alternate allele, distinct from
    \code{reference_allele}.}
  \item{build}{Genome build: \code{GRCh37}/\code{hg19} or
    \code{GRCh38}/\code{hg38}.}
}
\value{A character vector of canonical identity keys.}
\examples{
compressor_variant_key("1", 12345, "A", "G")
compressor_variant_key("chr23", 100, "A", "G", build = "GRCh37")
}
