\name{download_reference}
\alias{download_reference}
\title{Download a CompreSSoR reference}
\usage{
download_reference(reference = "GRCh38", cache_dir = NULL, overwrite = FALSE)
}
\arguments{
  \item{reference}{\code{"GRCh38"}, a reference descriptor, or a local reference path.}
  \item{cache_dir}{Optional cache directory.}
  \item{overwrite}{Re-download assets that are already cached.}
}
\description{
Downloads and verifies a reference asset, installing it in the local
CompreSSoR cache. Existing verified files are reused.
}
\value{A resolved reference descriptor containing local asset paths and checksums.}
