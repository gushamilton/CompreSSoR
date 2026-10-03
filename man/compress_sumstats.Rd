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
                  pvalue_threshold = 1e-05, pvalue_flag = NULL,
                  pvalue_flag_threshold = 5e-08, pvalue_order = NULL,
                  pvalue_order_threshold = NULL, region_padding = 10000L,
                  store_build = "GRCh38",
                  selection = c("full", "core", "hm3", "core_plus"),
                  threads = NULL, row_policy = c("report", "error"),
                  allele_columns = NULL)
}
\arguments{
  \item{input}{A data frame or delimited summary-statistics file.}
  \item{output}{Destination directory for the store.}
  \item{profile}{\code{"standard"} uses the semantic bounded-lossy
    Z10/EAF8/SE8 profile. \code{"exact"} is available with the Parquet backend.}
  \item{overwrite}{Whether an existing destination may be replaced.}
  \item{cache}{Whether to build the optional Parquet regional cache. Native
    Pcodec stores are already independently block-framed.}
  \item{strict}{Whether structural QC failures should stop the operation.
    This is an alias for \code{row_policy = "error"}.}
  \item{qc}{QC mode. \code{"compact"} is the default bounded-memory
    structural-QC path. \code{"none"} is an explicit fast path for already
    canonical prepared input; it bypasses structural QC, duplicate scans,
    row-level reports, alias resolution, and temporary variant-ID creation.
    A pre-canonicalization native-identity safety filter drops unsupported
    primary chromosomes, coordinates, or REF/ALT pairs and records their
    counts in preparation metadata. Malformed canonical statistics remain
    fail-closed; native duplicate checks still apply. Use it only for a
    table that is already validated upstream: canonical column names,
    primary chromosomes, whole positive coordinates, distinct single-base
    A/C/G/T REF/ALT with \code{effect_allele} = ALT, no duplicate variants,
    finite beta, positive finite SE, and EAF and p in [0, 1]; statistic
    values are not row-validated. On such a table the payload is identical
    to the \code{"compact"} result (FinnGen, 10M rows: about 26 s instead of
    41 s at 8 threads).}
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
  \item{pvalue_flag}{For native Pcodec stores, whether to write the aligned
    uint8 p-value flag domain. \code{NULL} uses the standard native default
    of \code{TRUE}; set \code{FALSE} to omit it.}
  \item{pvalue_flag_threshold}{Inclusive threshold for the aligned p-value
    flag domain. The standard convention is \code{5e-8}; this is separate
    from \code{pvalue_threshold}, which controls variant selection.}
  \item{pvalue_order}{For native Pcodec stores, whether to write the exact
    p-value ordering domain. \code{NULL} (default) uses the standard default of
    \code{TRUE}; \code{FALSE} omits it. One lossless uint32 rank per native row
    for candidates at \code{pvalue_order_threshold}; zero marks a non-candidate.
    Ranking and membership use the full-precision supplied p-value when present,
    otherwise p derived from exact prepared Z; ties use canonical variant order.
    Stores without the domain remain readable.}
  \item{pvalue_order_threshold}{Inclusive threshold for the exact p-value
    ordering domain. \code{NULL} uses 0.01.}
  \item{region_padding}{Base pairs added around p-value-region seeds. The
    default is 10,000 bp and endpoints are inclusive.}
  \item{store_build}{Build represented by the store. It must be the same as
    \code{input_build}; cross-build conversion is outside the package.}
  \item{selection}{Storage scope: \code{"full"}, \code{"core"},
    \code{"hm3"}, or \code{"core_plus"}.}
  \item{threads}{Positive worker count requested for compression and recorded
    in the manifest.}
  \item{row_policy}{\code{"report"} drops structurally invalid rows and
    records aggregate QC counts; \code{"error"} rejects them.}
  \item{allele_columns}{Optional mapping for inputs without explicit REF/ALT
    columns: \code{c(ref = "<column>", alt = "<column>")} names the input
    columns that hold REF and ALT. A GWAS-SSF table, whose
    \code{effect_allele} is ALT, is written with
    \code{allele_columns = c(ref = "other_allele", alt = "effect_allele")};
    the payload is identical to writing the table with
    \code{reference_allele}/\code{alternate_allele} columns added, without
    copying or rewriting it. The named columns must exist, and an input that
    already has an explicit REF or ALT column cannot also declare one. The
    core still requires \code{effect_allele} = ALT and \code{other_allele} =
    REF (no flipping). The mapping is recorded in the manifest's source
    provenance.}
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
