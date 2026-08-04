\name{validate_compressor}
\alias{validate_compressor}
\title{Validate a CompreSSoR Store}
\description{Checks the manifest, every Pcodec file checksum, indexes and,
optionally, every decoded frame.}
\usage{validate_compressor(store, full = FALSE)}
\arguments{
\item{store}{A \code{compressor_store} object or store directory path.}
\item{full}{For Pcodec stores, decode and semantically check every frame in
addition to verifying every file, frame, index, and manifest checksum.}
}
\value{A list with \code{valid}, \code{errors}, \code{rows}, and \code{profile}.}
\examples{
\dontrun{
validate_compressor("gwas.cpr", full = TRUE)
}
}
