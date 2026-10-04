\name{validate_compressor}
\alias{validate_compressor}
\title{Validate a CompreSSoR Store}
\description{Checks the manifest, every Pcodec file checksum, indexes and,
optionally, every decoded frame.}
\usage{validate_compressor(store, full = FALSE)}
\arguments{
\item{store}{A \code{compressor_store} object or store directory path.}
\item{full}{For Pcodec stores, verify every payload file against the byte
count and SHA-256 recorded in the manifest, check the aggregate
\code{payload_sha256} against that record, and decode and semantically check
every frame. Without \code{full}, only the manifest checksum (on open), the
native index checksum (whenever the index is parsed) and the structural
metadata are checked; a corrupted byte inside a value or identity stream is
then not detected.}
}
\value{A list with \code{valid}, \code{errors}, \code{rows}, and \code{profile}.}
\examples{
\dontrun{
validate_compressor("gwas.cpr", full = TRUE)
}
}
