\name{compressor_capabilities}
\alias{compressor_capabilities}
\title{Feature flags of this CompreSSoR build}
\usage{
compressor_capabilities()
}
\value{
A character vector of capability names.
}
\description{
Lets dependent packages feature-detect capabilities without version
parsing. \code{"candidates_one_pass"} means \code{\link[=read_candidates]{read_candidates()}} returns values,
\code{key}, \code{p_value} (bit-identical to \code{\link[=read_sumstats]{read_sumstats()}}) and exact ranks for
the candidate rows in a single pass, and \code{\link[=read_candidates_batch]{read_candidates_batch()}} reuses
same-panel identity. \code{"candidates_batch_rows_checked"} means
\code{\link[=read_candidates_batch]{read_candidates_batch()}} returns exactly the per-store \code{\link[=read_candidates]{read_candidates()}}
result for every batch composition (the shared-identity row loss in
0.7.0 is fixed) and both readers stop, rather than return partial data,
when a decode loses selected rows or a \code{pvalue_flag} read disagrees with
the manifest's flagged-row count.
}
