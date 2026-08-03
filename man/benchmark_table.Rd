\name{benchmark_table}
\alias{benchmark_table}
\title{Read the Shipped Benchmark Table}
\description{Reads one of the measured benchmark tables shipped with the package.}
\usage{benchmark_table(kind = c("pareto", "vcf", "finngen", "finngen_optimization",
  "modes", "release_gate", "release_gate_followup", "pcodec_access",
  "storage_size", "storage_amortization"))}
\arguments{\item{kind}{One of \code{"pareto"}, \code{"vcf"}, \code{"finngen"},
\code{"finngen_optimization"}, \code{"modes"}, \code{"release_gate"},
\code{"release_gate_followup"}, \code{"pcodec_access"}, \code{"storage_size"},
or \code{"storage_amortization"}. The \code{"pcodec_access"} table contains
the current v0.3 full-FinnGen exported-API runs; the others are retained
historical design and release-gate measurements.}}
\value{A data frame containing the selected benchmark.}
\seealso{\link{compress_sumstats}}
