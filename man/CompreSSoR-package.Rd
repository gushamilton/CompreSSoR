\name{CompreSSoR-package}
\alias{CompreSSoR}
\docType{package}
\title{Compressed GWAS Summary Statistics}
\description{
  CompreSSoR imports, lifts, harmonises, compresses, reads and validates
  genome-wide association study summary statistics. Its standard block-framed
  Pcodec store has a self-contained lossless GRCh38 SNV identity key and
  explicitly quantised Z9, EAF8 and SE6 streams. Beta and p-values are derived
  on read; Parquet is available as an interoperability backend.
}
\seealso{
  \code{\link{import_sumstats}}, \code{\link{liftover_sumstats}},
  \code{\link{harmonise_sumstats}}, \code{\link{compress_sumstats}},
  \code{\link{read_sumstats}}
}
\keyword{package}
