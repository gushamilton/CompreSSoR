\name{build_cache}
\alias{build_cache}
\title{Build the Optional Regional Cache}
\description{
  Builds an independent q8 framed cache beside a CompreSSoR store. Frames can
  be selected and decompressed for regional reads without scanning the full
  Parquet values file.
}
\usage{
build_cache(store, output = NULL, block_rows = 65536L, overwrite = FALSE)
}
\arguments{
  \item{store}{A \code{compressor_store} object or store directory path.}
  \item{output}{Optional cache directory; defaults to \code{cache.q8} inside the store.}
  \item{block_rows}{Rows per independent cache frame.}
  \item{overwrite}{Whether an existing cache may be replaced.}
}
\value{The cache directory path, invisibly.}
\seealso{\link{read_sumstats}}
