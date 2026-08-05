\name{harmonise_sumstats}
\alias{harmonise_sumstats}
\title{Harmonise Summary Statistics to GRCh38}
\description{
  Normalises common summary-statistics aliases, matches variants to a GRCh38
  reference spine, aligns alleles and flips beta and effect-allele frequency
  when required. The default drops unresolved variants and records their
  counts.
}
\usage{harmonise_sumstats(input, reference = "GRCh38", strict = FALSE,
                          mode = c("qc", "convert", "all", "core", "hm3",
                                   "pvalue_regions", "core_plus"),
                          variant_set = NULL, chrom_threads = 1L,
                          drop_unresolved = TRUE, input_build = "GRCh38",
                          chain = NULL, pvalue_threshold = 1e-05,
                          region_padding = 50000L)}
\arguments{
  \item{input}{A data frame or a delimited summary-statistics file.}
  \item{reference}{A GRCh38 reference data frame, \code{list(variants = ...)},
    a readable Parquet/TSV file, or \code{"GRCh38"} for the immutable reference
    configured by \code{COMPRESSOR_CANONICAL_REFERENCE}.}
  \item{strict}{Whether absent, incompatible, ambiguous, or duplicated
    variants should fail instead of being dropped and counted.}
  \item{drop_unresolved}{Whether unmatched, incompatible, ambiguous and
    duplicate rows are dropped. The default is \code{TRUE}.}
  \item{input_build}{Build of the incoming summary statistics. Non-GRCh38
    input requires \code{chain}.}
  \item{chain}{Optional GRCh37-to-GRCh38 chain file for liftover.}
  \item{pvalue_threshold}{Strict p-value threshold for p-value region or
    core-plus selection. P-values are derived from the canonical Z score.}
  \item{region_padding}{Base pairs added on each side of each significant SNP
    for p-value region or core-plus selection.}
  \item{mode}{One of \code{"qc"} (default), \code{"convert"}, \code{"all"},
    \code{"core"}, \code{"hm3"}, \code{"pvalue_regions"}, or
    \code{"core_plus"}.}
  \item{variant_set}{An optional panel data frame, file, or chromosome-shard
    directory. Staged shards are filtered to chromosomes present after
    harmonisation.}
  \item{chrom_threads}{Number of chromosome workers for harmonisation.}
}
\value{A normalized data frame with reference metadata attributes.}
\examples{
study <- data.frame(
  chromosome = "1", position = 100L, effect_allele = "A",
  other_allele = "C", beta = 0.2, standard_error = 0.1, eaf = 0.3
)
reference <- data.frame(
  chromosome = "1", base_pair_location = 100L,
  reference_allele = "A", alternate_allele = "C"
)
harmonise_sumstats(study, reference)
}
\seealso{\link{compress_sumstats}}
