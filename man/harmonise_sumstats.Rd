\name{harmonise_sumstats}
\alias{harmonise_sumstats}
\title{Harmonise Summary Statistics to GRCh38}
\description{
  Normalises common summary-statistics aliases, matches variants to a GRCh38
  reference spine, aligns alleles and flips beta and effect-allele frequency
  when required. The default preserves unresolved variants and records their
  harmonisation status.
}
\usage{harmonise_sumstats(input, reference = "GRCh38", strict = FALSE,
                          mode = c("qc", "convert", "all", "core", "hm3"),
                          variant_set = NULL, chrom_threads = 1L)}
\arguments{
  \item{input}{A data frame or a delimited summary-statistics file.}
  \item{reference}{A GRCh38 reference data frame, \code{list(variants = ...)},
    a readable Parquet/TSV file, or \code{"GRCh38"} for the automatic cached
    reference.}
  \item{strict}{Whether absent, incompatible, ambiguous, or duplicated
    variants should fail instead of being preserved and counted.}
  \item{mode}{One of \code{"qc"} (default), \code{"convert"}, \code{"all"},
    \code{"core"}, or \code{"hm3"}.}
  \item{variant_set}{A panel data frame or file for \code{mode = "core"} or
    \code{mode = "hm3"}. PLINK \code{.bim}, Parquet and delimited files are
    supported.}
  \item{chrom_threads}{Number of chromosome workers for harmonisation.}
}
\value{A normalized data frame with reference metadata attributes.}
\seealso{\link{compress_sumstats}}
