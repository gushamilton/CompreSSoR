\name{liftover_sumstats}
\alias{liftover_sumstats}
\title{Lift GWAS summary statistics to GRCh38}
\usage{
liftover_sumstats(
  input,
  input_build = "GRCh37",
  chain = NULL,
  drop_unmapped = TRUE
)
}
\arguments{
\item{input}{A data.frame or summary-statistics file accepted by
\code{import_sumstats()}.}

\item{input_build}{Source genome build, normally \code{"GRCh37"}. GRCh38
inputs are returned with status \code{"not_needed"}.}

\item{chain}{A local GRCh37-to-GRCh38 chain file. Required unless
\code{input_build} already denotes GRCh38.}

\item{drop_unmapped}{Drop rows that do not map exactly once.}
}
\value{
A canonical data.frame with a \code{liftover_status} column and
\code{liftover_stats} and \code{genome_build} attributes.
}
\description{
Uses a caller-supplied UCSC chain through Bioconductor's \code{rtracklayer}
and preserves an explicit row-level mapping status. Reverse-strand mappings
reverse-complement both alleles. Zero- and multi-mapped rows are dropped by
default and always counted in the \code{liftover_stats} attribute.
}
\seealso{
\code{\link{import_sumstats}}, \code{\link{harmonise_sumstats}}
}
\examples{
example <- data.frame(
  chromosome = "1", base_pair_location = 10177L,
  effect_allele = "A", other_allele = "C",
  beta = 0.03, standard_error = 0.04,
  effect_allele_frequency = 0.12
)
liftover_sumstats(example, input_build = "GRCh38")
}
