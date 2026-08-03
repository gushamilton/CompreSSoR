\name{read_sumstats}
\alias{read_sumstats}
\title{Read Decoded Summary Statistics}
\description{Reads decoded summary statistics from a CompreSSoR store.}
\usage{
read_sumstats(store, region = NULL, variants = NULL, columns = NULL,
              use_cache = FALSE)
}
\arguments{
  \item{store}{A \code{compressor_store} object or store directory path.}
  \item{region}{Optional region such as \code{"chr1:100-200"}.}
  \item{variants}{Optional variant IDs or zero-based row IDs.}
  \item{columns}{Optional output columns to return.}
  \item{use_cache}{Use an existing q8 cache for a regional read when available.}
}
\value{A data frame with ordinary summary-statistics columns.}
\examples{
\dontrun{
store <- open_compressor("gwas.cpr")
read_sumstats(store, region = "chr1:1000000-2000000",
              columns = c("beta", "standard_error", "p_value"))
}
}
\seealso{\link{compress_sumstats}, \link{decompress_sumstats}}
