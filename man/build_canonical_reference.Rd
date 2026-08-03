\name{build_canonical_reference}
\alias{build_canonical_reference}
\title{Build the Immutable Canonical GRCh38 Variant Reference}
\description{
  Reads an authoritative reference once, normalises and sorts it, and writes a
  single Parquet file whose primary key is \code{chrom:pos:ref:alt}. The
  resulting reference is external to study stores and is never synthesised by
  an individual study conversion.
}
\usage{
build_canonical_reference(source, output, source_build = "GRCh38", chain = NULL,
                          snps_only = FALSE, overwrite = FALSE)
}
\arguments{
  \item{source}{A reference data frame or VCF, PLINK BIM/PVAR, TSV or Parquet
    file. Compressed \code{.gz} and \code{.zst} PVAR sources are supported.}
  \item{output}{Path of the canonical Parquet file to write.}
  \item{source_build}{Build of \code{source}; non-GRCh38 inputs require a
    local GRCh37-to-GRCh38 chain supplied via \code{chain}.}
  \item{chain}{Optional GRCh37-to-GRCh38 chain file used by rtracklayer.}
  \item{snps_only}{Keep only single-base A/C/G/T variants. The default keeps
    SNVs, MNVs and indels with explicit sequence alleles.}
  \item{overwrite}{Whether an existing output may be replaced.}
}
\value{An immutable canonical reference descriptor.}
\seealso{\link{compress_sumstats}, \link{harmonise_sumstats}}
