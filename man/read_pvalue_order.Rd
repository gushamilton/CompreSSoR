\name{read_pvalue_order}
\alias{read_pvalue_order}
\title{Read exact p-value ordering metadata}
\description{
  Reads an optional exact p-value ordering domain from a native Pcodec store.
  Exact ranks are created before lossy encoding from finite supplied p-values
  or exact prepared Z, with canonical variant identity as the deterministic
  tie-break. Stores without the domain fail safely unless reconstructed-p
  fallback is requested explicitly.
}
\usage{
read_pvalue_order(store, threshold = NULL,
                  as = c("row_ids", "ranks"),
                  fallback = c("error", "reconstructed"),
                  name = "pvalue_order", threads = NULL)
}
\arguments{
  \item{store}{A native Pcodec store object or path.}
  \item{threshold}{Requested inclusive candidate threshold. \code{NULL} uses
    the exact domain threshold, or the flag threshold for explicit
    reconstructed fallback when no exact domain is present.}
  \item{as}{Return zero-based native row IDs in deterministic p-value order,
    or a full row-aligned integer rank vector where zero means
    non-candidate.}
  \item{fallback}{Safe default \code{"error"}, or explicit approximate
    ordering from p-values reconstructed from the lossy stored Z stream.}
  \item{name}{Domain name. The current supported name is
    \code{"pvalue_order"}.}
  \item{threads}{Number of independent ordering/value frames to decode.}
}
\value{An integer vector of zero-based row IDs in deterministic p-value
order, or a row-aligned integer rank vector.}
\seealso{\link{compress_sumstats}, \link{read_pvalue_flag},
  \link{read_sumstats}}
