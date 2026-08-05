\name{variant_panel}
\alias{write_variant_panel}
\alias{open_variant_panel}
\alias{read_variant_panel}
\alias{validate_variant_panel}
\title{Read and write native Pcodec variant panels}
\description{
  Variant panels use the exact build-aware CompreSSoR identity streams and a
  lossless binary membership flag. They are not summary-statistics stores and
  do not contain fabricated beta, standard-error, or Z values.
}
\usage{
write_variant_panel(data, output, build = "GRCh38",
                   source_provenance = NULL,
                   block_rows = PCODEC_NATIVE_BLOCK_ROWS,
                   key_block_rows = PCODEC_NATIVE_KEY_BLOCK_ROWS,
                   overwrite = FALSE)
open_variant_panel(path, verify_payload = FALSE)
read_variant_panel(panel, chromosomes = NULL, rows = NULL, keys = NULL,
                   hm3_only = FALSE, verify_payload = FALSE,
                   identity_only = FALSE)
validate_variant_panel(panel, full = FALSE)
}
\arguments{
  \item{data}{A data frame containing canonical variant identity and an
    optional binary \code{hm3} column.}
  \item{output}{Directory to create for the native panel store.}
  \item{build}{Genome build for the identity key.}
  \item{source_provenance}{Optional compact source hash and count metadata.}
  \item{block_rows, key_block_rows}{Power-of-two Pcodec frame sizes.}
  \item{overwrite}{Whether an existing output directory may be replaced.}
  \item{path}{A native variant-panel directory.}
  \item{panel}{A panel path or object returned by
    \code{open_variant_panel()}.}
  \item{verify_payload}{Whether to verify every payload checksum while
    opening or reading.}
  \item{chromosomes}{Optional chromosomes to read.}
  \item{rows}{Optional one-based panel row numbers.}
  \item{keys}{Optional canonical \code{chromosome:position:REF:ALT} keys.}
  \item{hm3_only}{Whether to return only rows whose HM3 flag is one.}
  \item{identity_only}{Whether to return numeric identity columns
    \code{global_position} and \code{substitution} instead of materialising
    character variant IDs. This is used by the compression core for bounded
    memory panel matching.}
  \item{full}{Whether to verify payload checksums and decode all streams.}
}
\value{
  \code{write_variant_panel()} and \code{open_variant_panel()} return a
  \code{compressor_variant_panel}; \code{read_variant_panel()} returns a
  data frame with \code{variant_id} and \code{hm3}; validation returns a
  list invisibly.
}
\seealso{\code{read_variant_set()}, \code{variant_set_membership()}}
