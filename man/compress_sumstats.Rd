\name{compress_sumstats}
\alias{compress_sumstats}
\title{Write a CompreSSoR Store}
\description{
  Normalises summary-statistics column names, optionally aligns variants to a
  GRCh38 reference, and writes a block-indexed Pcodec-backed store by default.
}
\usage{
compress_sumstats(input, output, reference = "GRCh38",
                  profile = c("standard", "exact"), overwrite = FALSE,
                  cache = FALSE, strict = FALSE, keep_extras = FALSE,
                  block_rows = 65536L,
                  mode = c("qc", "convert", "all", "core", "hm3",
                            "pvalue_regions", "core_plus"),
                  variant_set = NULL, chrom_threads = 1L,
                  drop_unresolved = TRUE, input_build = "GRCh38", chain = NULL,
                  backend = c("pcodec", "parquet"),
                  assume_grch38_ref_alt = FALSE, pvalue_threshold = 1e-05,
                  region_padding = 50000L)
}
\arguments{
  \item{input}{A data frame or a delimited summary-statistics file.}
  \item{output}{Destination directory for the store.}
  \item{reference}{A GRCh38 reference data frame, \code{list(variants = ...)},
    a local path, \code{"GRCh38"} for the immutable reference configured by
    \code{COMPRESSOR_CANONICAL_REFERENCE}, or \code{NULL} to skip reference
    alignment.}
\item{profile}{\code{"standard"} uses semantic Z9/EAF8/SE6 streams with
    sparse float32 exceptions; \code{"exact"} is available with the Parquet
    backend. P and beta are derived rather than stored physically.}
  \item{backend}{Storage backend. The default \code{"pcodec"} is the compact
    block-framed format with a self-contained GRCh38 identity key;
    \code{"parquet"} is the interoperable legacy backend.}
  \item{assume_grch38_ref_alt}{Required for Pcodec \code{mode = "convert"}
    unless the input has explicit REF and ALT columns. Set \code{TRUE} only
    when positions are GRCh38, \code{other_allele} is REF,
    \code{effect_allele} is ALT, and beta/Z/EAF refer to ALT.}
  \item{overwrite}{Whether an existing destination may be replaced.}
  \item{cache}{Whether to build the optional q8 framed regional cache.}
  \item{strict}{Whether absent, incompatible, ambiguous, duplicated, or
    backend-unsupported variants should fail instead of being dropped and
    counted in the manifest.}
  \item{drop_unresolved}{Whether unmatched, incompatible, ambiguous and
    duplicate rows are dropped. The default is \code{TRUE}.}
  \item{input_build}{Build of the incoming summary statistics. Non-GRCh38
    input requires \code{chain}.}
  \item{chain}{Optional GRCh37-to-GRCh38 chain file for liftover.}
  \item{pvalue_threshold}{Strict p-value threshold for p-value region or
    core-plus selection; p-values are derived from canonical Z.}
  \item{region_padding}{Base pairs added on each side of significant SNPs for
    p-value region or core-plus selection.}
  \item{keep_extras}{Whether arbitrary non-core input columns should be kept
    in an extras sidecar. The default is \code{FALSE}; harmonisation/QC
    summaries remain in the manifest, while row-wise QC labels are retained
    only when this is explicitly set to \code{TRUE}.}
  \item{block_rows}{Rows per Parquet write chunk; Pcodec uses its measured
    fixed key and value frame sizes.}
  \item{mode}{One of \code{"qc"} (default), \code{"convert"}, \code{"all"},
    \code{"core"}, \code{"hm3"}, \code{"pvalue_regions"}, or
    \code{"core_plus"}.}
  \item{variant_set}{An optional panel data frame, file, or chromosome-shard
    directory. In \code{mode = "core_plus"}, it is unioned with the padded
    p-value regions after harmonisation.}
  \item{chrom_threads}{Number of chromosome workers for harmonisation.}
}
\value{A \code{compressor_store} object.}
\examples{
\dontrun{
example <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
store <- compress_sumstats(
  example, file.path(tempdir(), "example.cpr"), mode = "convert",
  reference = NULL, assume_grch38_ref_alt = TRUE, overwrite = TRUE
)
read_sumstats(store, columns = c("chromosome", "base_pair_location", "z"))
}
}
\seealso{\link{read_sumstats}, \link{validate_compressor}}
