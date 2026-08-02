\name{compress_sumstats}
\alias{compress_sumstats}
\title{Write a CompreSSoR Store}
\description{
  Normalises summary-statistics column names, optionally aligns variants to a
  GRCh38 reference spine, and writes a block-indexed Parquet-backed store.
}
\usage{
compress_sumstats(input, output, reference = "GRCh38",
                  profile = c("standard", "exact"), overwrite = FALSE,
                  cache = FALSE, strict = FALSE, block_rows = 65536L,
                  mode = c("qc", "convert", "all", "core", "hm3"),
                  variant_set = NULL, chrom_threads = 1L)
}
\arguments{
  \item{input}{A data frame or a delimited summary-statistics file.}
  \item{output}{Destination directory for the store.}
  \item{reference}{A GRCh38 reference data frame, \code{list(variants = ...)},
    a local path, \code{"GRCh38"} for the automatic cached reference, or
    \code{NULL} to skip reference alignment.}
  \item{profile}{\code{"standard"} uses q9/SE10 compact Parquet;
    \code{"exact"} stores numeric values without quantisation.}
  \item{overwrite}{Whether an existing destination may be replaced.}
  \item{cache}{Whether to build the optional q8 framed regional cache.}
  \item{strict}{Whether absent, incompatible, ambiguous, or duplicated
    variants should fail instead of being preserved and counted in the
    manifest.}
  \item{block_rows}{Rows per Parquet write chunk.}
  \item{mode}{One of \code{"qc"} (default), \code{"convert"}, \code{"all"},
    \code{"core"}, or \code{"hm3"}.}
  \item{variant_set}{A panel data frame or file for \code{mode = "core"} or
    \code{mode = "hm3"}. PLINK \code{.bim}, Parquet and delimited files are
    supported.}
  \item{chrom_threads}{Number of chromosome workers for harmonisation.}
}
\value{A \code{compressor_store} object.}
\seealso{\link{read_sumstats}, \link{validate_compressor}}
