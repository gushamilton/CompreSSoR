\name{resolve_reference}
\alias{resolve_reference}
\title{Resolve a CompreSSoR reference}
\usage{
resolve_reference(reference = "GRCh38", cache_dir = NULL)
}
\arguments{
  \item{reference}{\code{"GRCh38"}, a reference descriptor, or a local
  reference path. \code{"GRCh38"} requires a preconfigured local dictionary;
  only descriptors with an explicit URL are downloaded.}
  \item{cache_dir}{Optional cache directory.}
}
\description{Resolves a configured reference descriptor and downloads an
explicit remote asset when necessary.}
\value{A resolved reference descriptor.}
