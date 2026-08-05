\name{import_sumstats}
\alias{import_sumstats}
\title{Import GWAS summary statistics}
\usage{
import_sumstats(input, strict = FALSE, row_policy = c("report", "error"),
                input_build = "GRCh38", allow_p_to_se = FALSE)
}
\arguments{
\item{input}{A data.frame or an existing summary-statistics file.}
\item{strict}{Whether structural QC failures should stop import. The default
is \code{FALSE}; canonical rows are returned with a bounded structural QC
report attribute.}
\item{row_policy}{Either \code{"report"} to continue with a report or
\code{"error"} to stop on any structural rejection. \code{strict = TRUE} is
an alias for \code{"error"}.}
\item{input_build}{Build used for coordinate-bound preflight: \code{"GRCh37"}
or \code{"GRCh38"}, with \code{"hg19"}/\code{"hg38"} aliases.}
\item{allow_p_to_se}{Explicitly enable conflict-checked p-value-to-standard-error
conversion at the import boundary. This is disabled in the strict compression
core and is recorded in resolution provenance.}
}
\value{
A data.frame containing canonical columns such as \code{chromosome},
\code{base_pair_location}, \code{reference_allele},
\code{alternate_allele}, \code{effect_allele}, \code{other_allele}, \code{beta},
\code{standard_error}, \code{z},
\code{effect_allele_frequency}, and \code{p_value}.
The \code{structural_qc_report} attribute records rejected-row counts and
bounded row-number examples. Use the internal preparation helper
\code{preflight_sumstats()} when the filtered data and report are both needed
before strict compression or an external preparation workflow.
}
\description{
Reads a data frame or delimited text file with automatic comma, tab,
semicolon, and whitespace delimiter handling, including gzip/bgzip, bzip2,
xz, and zstd transports. Simple single-trait VCF inputs are converted from
fixed columns plus INFO/FORMAT values. Common aliases, including UKB-PPP
\code{GENPOS}, \code{ALLELE0}, \code{ALLELE1}, and \code{LOG10P}, are mapped onto
CompreSSoR's canonical in-memory schema and unrecognised columns are retained.
MalariaGEN-style \code{alleleA}/\code{alleleB} are mapped to the other/effect
alleles, while a single-letter \code{B} is retained as \code{expected_count}
and is never treated as a beta alias when an exact \code{beta} column is
present. GWAS Catalog \code{hm_*} fields are detected for a structured
orientation diagnostic but are not inferred as REF/ALT automatically.
Chromosome prefixes and 23/24 are normalised to 1--22/X/Y. This function does
not change genome build or allele orientation. External preparation tools must
provide explicit REF/ALT and effect/non-effect allele columns before
compression.
}
\examples{
path <- system.file("extdata", "example-grch38.tsv", package = "CompreSSoR")
example <- import_sumstats(path)
head(example[c("chromosome", "base_pair_location", "beta", "standard_error")])
}
