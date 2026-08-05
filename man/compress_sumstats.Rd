\name{compress_sumstats}
\alias{compress_sumstats}
\title{Write a CompreSSoR Store}
\description{
  Imports summary statistics and writes a CompreSSoR store. The canonical
  public path is a self-contained GRCh38 Pcodec store: pass already verified
  GRCh38 REF/ALT data with \code{reference = NULL}. Supplying a reference
  enables optional reference-backed harmonisation before selection and
  compression. Parquet is the explicit interoperability backend for exact
  doubles and arbitrary extra columns.
}
\usage{
compress_sumstats(input, output, reference = NULL,
                  profile = c("standard", "exact"), overwrite = FALSE,
                  cache = FALSE, strict = NULL, keep_extras = FALSE,
                  block_rows = 65536L, mode = NULL, variant_set = NULL,
                  chrom_threads = NULL, drop_unresolved = NULL,
                  input_build = "GRCh38", chain = NULL,
                  backend = c("pcodec", "parquet"),
                  assume_grch38_ref_alt = FALSE,
                  pvalue_threshold = 1e-05, region_padding = 50000L,
                  store_build = "GRCh38",
                  selection = c("full", "core", "hm3", "core_plus"),
                  threads = NULL, row_policy = c("report", "error"))
}
\arguments{
  \item{input}{A data frame or delimited summary-statistics file.}
  \item{output}{Destination directory for the store.}
  \item{reference}{Optional GRCh38 reference data frame, descriptor, or
    readable reference file. \code{NULL} selects the self-contained path and
    skips reference alignment.}
  \item{store_build}{Build represented by the store. The current public
    identity interface supports \code{"GRCh37"} and \code{"GRCh38"}; the
    native writer records the selected build and rejects unsupported reverse
    liftover combinations.}
  \item{selection}{Storage scope: \code{"full"}, \code{"core"},
    \code{"hm3"}, or \code{"core_plus"}. Selection is applied after import,
    structural preparation, and any requested harmonisation. Named core/HM3
    panels are resolved through the existing panel interface.}
  \item{threads}{Positive scalar integer requested for preparation. It is
    passed to eligible chromosome harmonisation and its requested/effective
    values are recorded in the manifest.}
  \item{row_policy}{Either \code{"report"}, which drops unresolved rows and
    records aggregate counts, or \code{"error"}, which rejects them.}
  \item{profile}{\code{"standard"} uses the semantic bounded-lossy numeric
    profile; \code{"exact"} is available with the Parquet backend. Beta and
    p-values are derived rather than stored in either profile.}
  \item{backend}{Storage backend. \code{"pcodec"} is the compact native
    self-contained format; \code{"parquet"} is the interoperability backend.}
  \item{mode}{Compatibility alias for the historical overloaded argument.
    It accepts \code{"qc"}, \code{"convert"}, \code{"all"}, \code{"core"},
    \code{"hm3"}, \code{"pvalue_regions"}, and \code{"core_plus"}. New code
    should use \code{reference}, \code{selection}, and \code{row_policy}
    independently.}
  \item{variant_set}{Optional panel data frame, file, or chromosome-shard
    directory. If omitted for \code{selection = "core"}, \code{"hm3"}, or
    \code{"core_plus"}, the corresponding existing named-panel resolver is
    used.}
  \item{input_build}{Build of the input rows. Non-matching builds require a
    compatible chain when the selected store build needs liftover.}
  \item{chain}{Optional GRCh37-to-GRCh38 chain file. Its path and checksum are
    recorded in provenance metadata.}
  \item{strict}{Legacy alias for \code{row_policy = "error"}.}
  \item{drop_unresolved}{Legacy compatibility switch controlling whether
    unresolved rows are dropped before writing. Pcodec stores cannot retain
    unresolved identity rows.}
  \item{chrom_threads}{Legacy alias for \code{threads}.}
  \item{assume_grch38_ref_alt}{For reference-free Pcodec conversion, asserts
    that coordinates are GRCh38, \code{other_allele} is REF,
    \code{effect_allele} is ALT, and statistics refer to ALT. Explicit REF
    and ALT columns establish the same orientation contract.}
  \item{pvalue_threshold}{Strict threshold for p-value-region or core-plus
    selection; p-values are derived from the pre-encoding canonical Z score.}
  \item{region_padding}{Base pairs added on each side of significant variants
    for p-value-region or core-plus selection.}
  \item{overwrite}{Whether an existing destination may be replaced.}
  \item{cache}{Whether to build the optional Parquet regional cache. Native
    Pcodec stores are already independently block-framed.}
  \item{keep_extras}{Whether arbitrary non-core input columns should be kept
    in a Parquet extras sidecar.}
  \item{block_rows}{Rows per Parquet write chunk; native Pcodec uses its
    measured key and value frame sizes.}
}
\value{A \code{compressor_store} object.}
\examples{
\dontrun{
example <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
store <- compress_sumstats(
  example, file.path(tempdir(), "example.cpr"),
  store_build = "GRCh38", selection = "full", reference = NULL,
  assume_grch38_ref_alt = TRUE, overwrite = TRUE
)
read_sumstats(store, columns = c("chromosome", "base_pair_location", "z"))
}
}
\seealso{\link{read_sumstats}, \link{validate_compressor},
  \link{harmonise_sumstats}}
