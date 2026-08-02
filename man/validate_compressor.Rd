\name{validate_compressor}
\alias{validate_compressor}
\title{Validate a CompreSSoR Store}
\description{Checks the manifest, referenced files, and decoded row count.}
\usage{validate_compressor(store)}
\arguments{\item{store}{A \code{compressor_store} object or store directory path.}}
\value{A list with \code{valid}, \code{errors}, \code{rows}, and \code{profile}.}
