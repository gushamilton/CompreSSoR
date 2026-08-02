\name{benchmark_table}
\alias{benchmark_table}
\title{Read the Shipped Benchmark Table}
\description{Reads one of the measured benchmark tables shipped with the package.}
\usage{benchmark_table(kind = c("pareto", "vcf", "finngen", "finngen_optimization",
  "modes", "release_gate", "release_gate_followup"))}
\arguments{\item{kind}{One of \code{"pareto"}, \code{"vcf"}, \code{"finngen"}, \code{"finngen_optimization"}, \code{"modes"}, \code{"release_gate"}, or \code{"release_gate_followup"}.}}
\value{A data frame containing the selected benchmark.}
\seealso{\link{compress_sumstats}}
