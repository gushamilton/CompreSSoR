\name{download_reference}
\alias{download_reference}
\title{Download a CompreSSoR reference}
\usage{
download_reference(reference = "GRCh38", cache_dir = NULL, overwrite = FALSE)
}
\arguments{
  \item{reference}{\code{"GRCh38"}, a reference descriptor, or a local
  reference path. The built-in descriptor requires
  \code{COMPRESSOR_CANONICAL_REFERENCE} to point at a dictionary built by
  \code{build_ebi_reference()}; it has no implicit remote URL.}
  \item{cache_dir}{Optional cache directory.}
  \item{overwrite}{Re-download assets that are already cached.}
}
\description{
Resolves and verifies a configured reference asset, downloading only when an
explicit descriptor supplies a URL. Existing verified files are reused.
}
\value{A resolved reference descriptor containing local asset paths and checksums.}
