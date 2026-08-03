\name{compressor_variant_key}
\alias{compressor_variant_key}
\title{Construct Canonical CompreSSoR Variant Keys}
\description{Constructs self-contained GRCh38 biallelic-SNV keys in
\code{chromosome:position:REF:ALT} form. CompreSSoR stores effects and
frequencies for ALT.}
\usage{
compressor_variant_key(chromosome, position, reference_allele, alternate_allele)
}
\arguments{
  \item{chromosome}{GRCh38 chromosome (1--22, X, or Y).}
  \item{position}{One-based GRCh38 position.}
  \item{reference_allele}{GRCh38 reference allele.}
  \item{alternate_allele}{Alternate allele.}
}
\value{A character vector of canonical identity keys.}
\examples{
compressor_variant_key("1", 12345, "A", "G")
}
