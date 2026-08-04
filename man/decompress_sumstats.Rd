\name{decompress_sumstats}
\alias{decompress_sumstats}
\title{Decompress a CompreSSoR Store}
\description{Alias for \code{read_sumstats}; decodes stored values on read.}
\usage{
decompress_sumstats(store, region = NULL, variants = NULL, columns = NULL,
                    use_cache = FALSE, threads = NULL)
}
\arguments{
  \item{store}{A \code{compressor_store} object or store directory path.}
  \item{region}{Optional region such as \code{"chr1:100-200"}.}
  \item{variants}{Optional variant IDs, zero-based row IDs, or, for Pcodec
  stores, canonical \code{chromosome:position:REF:ALT} keys.}
  \item{columns}{Optional output columns to return.}
  \item{use_cache}{Use an existing q8 cache for a regional read when available.}
  \item{threads}{Number of Unix worker processes for independent Pcodec
  stream/block reads. \code{NULL} selects four workers for whole-store reads
  and one worker for regional or canonical-key reads. Set \code{threads}
  explicitly, or use \code{options(CompreSSoR.pcodec.threads = n)}. Windows
  remains serial.}
}
\value{A data frame with ordinary summary-statistics columns.}
\examples{
\dontrun{
decoded <- decompress_sumstats("gwas.cpr")
}
}
