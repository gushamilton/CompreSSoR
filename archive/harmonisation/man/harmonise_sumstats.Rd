\name{harmonise_sumstats}
\alias{harmonise_sumstats}
\title{Harmonise Summary Statistics to a Target Build}
\description{
  Imports summary statistics, performs structural preparation and optional
  reference-backed allele harmonisation, and returns a canonical data frame
  suitable for \code{compress_sumstats()}. The returned object carries
  alignment, reference, build, and provenance attributes so it can be passed
  to the self-contained compression path without repeating harmonisation.
}
\usage{
harmonise_sumstats(input, reference = "GRCh38", strict = FALSE,
                   mode = c("qc", "convert", "all", "core", "hm3",
                             "pvalue_regions", "core_plus"),
                   variant_set = NULL, chrom_threads = 1L,
                   drop_unresolved = TRUE, input_build = "GRCh38",
                   chain = NULL, pvalue_threshold = 1e-05,
                   region_padding = 50000L, target_build = "GRCh38",
                   qc = NULL, liftover = NULL, strand = NULL,
                   palindromic = NULL, frequency_qc = NULL,
                   frequency_tolerance = NULL, max_examples = NULL)
}
\arguments{
  \item{input}{A data frame or delimited summary-statistics file.}
  \item{reference}{A target-build reference data frame, descriptor, readable
    reference file, or \code{"GRCh38"} for the configured canonical
    reference. Use \code{NULL} for reference-free canonical conversion.}
  \item{input_build}{Build of the incoming summary statistics. A non-target
    build requires \code{chain} when liftover is needed.}
  \item{target_build}{Build represented by the returned coordinates. The
    current reference-backed helper targets GRCh38; already aligned GRCh37
    input may be passed through with \code{reference = NULL}.}
  \item{mode}{Compatibility preparation alias. \code{"qc"} aligns to the
    reference and reports/drops unresolved rows, while \code{"convert"}
    skips reference alignment. The selection-bearing legacy values
    \code{"core"}, \code{"hm3"}, \code{"pvalue_regions"}, and
    \code{"core_plus"} are retained for compatibility; new compression code
    should express storage scope with \code{selection}.}
  \item{variant_set}{Optional panel data frame, file, or chromosome-shard
    directory passed to the existing panel resolver after preparation.}
  \item{strict}{Whether unresolved, incompatible, ambiguous, or duplicate
    rows should raise an error instead of being dropped and counted.}
  \item{drop_unresolved}{Whether unresolved rows are dropped. The default is
    \code{TRUE}; setting \code{FALSE} retains them for an audit-oriented
    Parquet path.}
  \item{chrom_threads}{Number of chromosome workers used by eligible
    reference harmonisation.}
  \item{chain}{Optional GRCh37-to-GRCh38 chain file. Its path and checksum are
    retained in the returned alignment attributes and downstream manifest.}
  \item{pvalue_threshold}{Strict threshold for p-value-region or core-plus
    selection, using the canonical pre-encoding Z score.}
  \item{region_padding}{Base pairs added on both sides of significant variants
    for p-value-region or core-plus selection.}
  \item{qc}{Optional list controlling strand, palindromic, frequency, and
    liftover quality-control policies.}
  \item{liftover}{Optional scalar override for \code{qc$liftover}.}
  \item{strand}{Optional scalar override for \code{qc$strand}.}
  \item{palindromic}{Optional override for \code{qc$palindromic}.}
  \item{frequency_qc}{Optional override for \code{qc$frequency}.}
  \item{frequency_tolerance}{Optional override for
    \code{qc$frequency_tolerance}.}
  \item{max_examples}{Optional override for \code{qc$example_limit}.}
}
\value{A normalized data frame with reference, alignment, build, and
  identity attributes.}
\examples{
study <- data.frame(
  chromosome = "1", position = 100L, effect_allele = "A",
  other_allele = "C", beta = 0.2, standard_error = 0.1, eaf = 0.3
)
reference <- data.frame(
  chromosome = "1", base_pair_location = 100L,
  reference_allele = "A", alternate_allele = "C"
)
harmonised <- harmonise_sumstats(study, reference = reference,
                                  input_build = "GRCh38",
                                  target_build = "GRCh38")
}
\seealso{\link{compress_sumstats}}
