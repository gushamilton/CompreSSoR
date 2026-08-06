\name{read_pvalue_flag}
\alias{read_pvalue_flag}
\title{Read an aligned p-value flag domain}
\description{
  Reads the optional aligned p-value flag domain from a native Pcodec store.
  The standard native domain is written at p <= 5e-8 and uses the immutable
  zero-based row order of the store.
}
\usage{
read_pvalue_flag(store, name = "pvalue_flag",
                 as = c("row_ids", "logical"), threads = NULL)
}
\arguments{
  \item{store}{A native Pcodec store object or path.}
  \item{name}{Domain name. The current supported name is
    \code{"pvalue_flag"}.}
  \item{as}{Return zero-based row IDs, or a logical vector aligned to all
    store rows.}
  \item{threads}{Number of independent flag frames to decode.}
}
\value{An integer vector of zero-based row IDs, or a logical vector with one
value per store row.}
\seealso{\link{compress_sumstats}, \link{read_sumstats}}
