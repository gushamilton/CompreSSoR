\name{compress_sumstats}
\alias{compress_sumstats}
\title{Write a CompreSSoR Store}
\description{
  Writes a compact CompreSSoR store from already-prepared summary statistics.
  The input must contain chromosome, position, explicit REF and ALT,
  effect/non-effect alleles, beta and standard error, or a positive OR and
  standard error that resolve to log-OR. The core performs
  structural QC and compression only: it does not look up a reference, flip
  alleles, resolve rsIDs, or liftover coordinates. GRCh37 and GRCh38 are
  supported when the input and store builds are the same.
}
\usage{
compress_sumstats(input, output, profile = c("standard", "exact"),
                  overwrite = FALSE, cache = FALSE, strict = NULL,
                  qc = c("compact", "none"),
                  keep_extras = FALSE, block_rows = 65536L,
                  variant_set = NULL, input_build = "GRCh38",
                  backend = c("pcodec", "parquet"),
                  pvalue_threshold = 1e-05, region_padding = 50000L,
                  store_build = "GRCh38",
                  selection = c("full", "core", "hm3", "core_plus"),
                  threads = NULL, row_policy = c("report", "error"))
}
\arguments{
  \item{input}{A data frame or delimited summary-statistics file.}
  \item{output}{Destination directory for the store.}
  \item{profile}{\code{"standard"} uses the semantic bounded-lossy
    Z9/EAF8/SE6 profile. \code{"exact"} is available with the Parquet backend.}
  \item{overwrite}{Whether an existing destination may be replaced.}
  \item{cache}{Whether to build the optional Parquet regional cache. Native
    Pcodec stores are already independently block-framed.}
  \item{strict}{Whether structural QC failures should stop the operation.
    This is an alias for \code{row_policy = "error"}.}
  \item{qc}{QC mode. \code{"compact"} is the default bounded-memory
    structural-QC path. \code{"none"} is an explicit fast path for already
    canonical prepared input; it bypasses structural QC, duplicate scans,
    row-level reports, alias resolution, and temporary variant-ID creation.
    Build-aware identity and native codec safety checks still apply.}
  \item{keep_extras}{Whether arbitrary non-core input columns should be kept
    in a Parquet extras sidecar.}
  \item{block_rows}{Rows per Parquet write chunk; native Pcodec uses its
    measured key and value frame sizes.}
  \item{variant_set}{Optional prevalidated panel used for post-import
    filtering with \code{selection = "core"}, \code{"hm3"}, or
    \code{"core_plus"}. It is not used as an allele reference.}
  \item{input_build}{Build of the input rows: \code{"GRCh37"}/\code{"hg19"}
    or \code{"GRCh38"}/\code{"hg38"}.}
  \item{backend}{\code{"pcodec"} is the compact native format;
    \code{"parquet"} is the interoperability backend.}
  \item{pvalue_threshold}{Threshold used by p-value-region or core-plus
    selection, derived from the prepared Z statistic.}
  \item{region_padding}{Base pairs added around p-value-region seeds.}
  \item{store_build}{Build represented by the store. It must be the same as
    \code{input_build}; cross-build conversion is outside the package.}
  \item{selection}{Storage scope: \code{"full"}, \code{"core"},
    \code{"hm3"}, or \code{"core_plus"}.}
  \item{threads}{Positive worker count requested for compression and recorded
    in the manifest.}
  \item{row_policy}{\code{"report"} drops structurally invalid rows and
    records aggregate QC counts; \code{"error"} rejects them.}
}
\value{A \code{compressor_store} object.}
\examples{
\dontrun{
example <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
store <- compress_sumstats(
  example, file.path(tempdir(), "example.cpr"),
  input_build = "GRCh38", store_build = "GRCh38", overwrite = TRUE
)
read_sumstats(store, columns = c("chromosome", "base_pair_location", "z"))
}
}
\seealso{\link{read_sumstats}, \link{validate_compressor}}
