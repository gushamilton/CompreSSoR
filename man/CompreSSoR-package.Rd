\name{CompreSSoR-package}
\alias{CompreSSoR}
\docType{package}
\title{Compressed GWAS Summary Statistics}
\description{
  CompreSSoR validates and compresses prepared genome-wide association study
  summary statistics. Its standard block-framed Pcodec store has a
  self-contained lossless GRCh37 or GRCh38 SNV identity key and explicitly
  quantised Z9, EAF8 and SE6 streams. Beta and p-values are derived on read;
  Parquet is available as an interoperability backend. Reference lookup,
  allele harmonisation, and liftover are outside the installed core.
}
\seealso{
  \code{\link{import_sumstats}}, \code{\link{compress_sumstats}},
  \code{\link{read_sumstats}}
}
\keyword{package}
