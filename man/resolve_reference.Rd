\name{resolve_reference}
\alias{resolve_reference}
\title{Resolve a CompreSSoR reference}
\usage{
resolve_reference(reference = "GRCh38", cache_dir = NULL)
}
\arguments{
  \item{reference}{\code{"GRCh38"}, a reference descriptor, or a local reference path.}
  \item{cache_dir}{Optional cache directory.}
}
\description{Resolves a reference descriptor and downloads missing assets.}
\value{A resolved reference descriptor.}
