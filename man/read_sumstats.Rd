\name{read_sumstats}
\alias{read_sumstats}
\title{Read Decoded Summary Statistics}
\description{Reads decoded summary statistics from a CompreSSoR store.}
\usage{
read_sumstats(store, region = NULL, variants = NULL, columns = NULL,
              use_cache = FALSE, threads = NULL)
}
\arguments{
  \item{store}{A \code{compressor_store} object or store directory path.}
  \item{region}{Optional region such as \code{"chr1:100-200"}.}
  \item{variants}{Optional variant IDs, zero-based row IDs, or, for Pcodec
  stores, canonical \code{chromosome:position:REF:ALT} keys.}
  \item{columns}{Optional output columns to return. Pcodec stores also expose
  \code{global_position} and \code{substitution}, the compact identity fields
  used by native analytical readers.}
  \item{use_cache}{Use an existing q8 cache for a regional read when available.}
  \item{threads}{Number of workers for Pcodec reads. The native full-store
  reader performs indexed stream scans and block decompression in C++; the
  regional, canonical-key and fallback paths use Unix worker processes. Set
  \code{threads} explicitly, or use \code{options(CompreSSoR.pcodec.threads = n)}
  for the non-native paths. Windows remains serial.}
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
