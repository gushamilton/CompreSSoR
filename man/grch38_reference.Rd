\name{grch38_reference}
\alias{grch38_reference}
\title{Describe the default GRCh38 reference}
\usage{
grch38_reference(cache_dir = NULL)
}
\arguments{
  \item{cache_dir}{Optional cache directory to use when resolving the descriptor.}
}
\description{
Returns the external EBI/Ensembl release-95 GRCh38 reference descriptor.
Build the complete partitioned dictionary with
\code{\link{build_ebi_reference}} and set
\code{COMPRESSOR_CANONICAL_REFERENCE} to its output directory.
}
\value{A reference descriptor.}
