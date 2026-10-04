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
the candidate rows in a single pass, and \code{\link[=read_candidates_batch]{read_candidates_batch()}} can
reuse same-panel identity. \code{"candidates_batch_rows_checked"} means
\code{\link[=read_candidates_batch]{read_candidates_batch()}} returns exactly the per-store \code{\link[=read_candidates]{read_candidates()}}
result for every batch composition (the shared-identity row loss in
0.7.0 is fixed) and both readers stop, rather than return partial data,
when a decode loses selected rows or a \code{pvalue_flag} read disagrees with
the manifest's flagged-row count. \code{"reads_bit_identical"} means every
reader (full, row-ID, key, region, candidate and batched reads) decodes a
row with the same native decoder, so values are bit-identical to the full
read. \code{"integrity_verified"} means \code{validate_compressor(full = TRUE)}
checks every payload checksum, the native index is checked whenever it is
parsed, and batched identity sharing verifies each member's identity
streams. \code{"key_reads_identity_columns"} means key, row-ID and region reads
(single and batched) return the numeric identity columns \code{global_position}
and \code{substitution} when requested. \code{"request_index"} means
\code{\link[=read_sumstats_batch]{read_sumstats_batch()}} accepts \code{request_index = TRUE}, and
\code{"identity_code"} that \code{\link[=compressor_identity_code]{compressor_identity_code()}} is exported.
}
