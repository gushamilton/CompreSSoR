\name{import_sumstats}
\alias{import_sumstats}
\title{Import GWAS summary statistics}
\usage{
import_sumstats(input)
}
\arguments{
\item{input}{A data.frame or an existing summary-statistics file.}
}
\value{
A data.frame containing canonical columns such as \code{chromosome},
\code{base_pair_location}, \code{effect_allele}, \code{other_allele},
\code{beta}, \code{standard_error}, \code{z},
\code{effect_allele_frequency}, and \code{p_value}.
}
\description{
Reads a data frame, delimited text file, gzip-compressed text file, or
single-ALT VCF and maps common GWAS column names onto CompreSSoR's canonical
in-memory schema. Unrecognised columns are retained. This function does not
change genome build or allele orientation; use \code{liftover_sumstats()} and
\code{harmonise_sumstats()} for those operations.
}
\seealso{
\code{\link{liftover_sumstats}}, \code{\link{harmonise_sumstats}}
}
\examples{
path <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
example <- import_sumstats(path)
head(example[c("chromosome", "base_pair_location", "beta", "standard_error")])
}
