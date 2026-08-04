\name{benchmark_table}
\alias{benchmark_table}
\title{Read the Shipped Benchmark Table}
\description{Reads one of the measured benchmark tables shipped with the package.}
\usage{benchmark_table(kind = c("finngen_10m", "pcodec_10m_threads", "pareto", "vcf", "finngen", "finngen_optimization",
  "modes", "release_gate", "release_gate_followup", "pcodec_access",
  "storage_size", "storage_amortization", "native_se8", "pareto_chr1"))}
\arguments{\item{kind}{One of \code{"pareto"}, \code{"vcf"}, \code{"finngen"},
\code{"finngen_optimization"}, \code{"modes"}, \code{"release_gate"},
\code{"release_gate_followup"}, \code{"pcodec_access"}, \code{"storage_size"},
or \code{"storage_amortization"}, \code{"native_se8"}, or
\code{"pareto_chr1"}. Use \code{"finngen_10m"} for the current keyed
10-million-row FinnGen format comparison and \code{"pcodec_10m_threads"}
for its native one-versus-four-thread reader comparison. The
\code{"pcodec_access"} table contains
the current v0.3 full-FinnGen exported-API runs; the others are retained
historical design and release-gate measurements.}}
\value{A data frame containing the selected benchmark.}
\seealso{\link{compress_sumstats}}
