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
Returns the external GRCh38 reference descriptor used by the automatic
reference cache. The reference is downloaded on demand and is not bundled in
the package.
}
\value{A reference descriptor.}
