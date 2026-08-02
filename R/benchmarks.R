benchmark_metadata <- function() {
  list(
    benchmark_id = "first_pareto_10m_full_read",
    description = "Measured five-run full-read benchmark on a real 10-million-row GWAS",
    rows = 10000000L,
    runs = 5L,
    metric = "median_full_read_seconds",
    compression_metric = "compression_ratio_vs_source_gzip",
    data = "inst/benchmarks/first-pareto.csv",
    plot = "inst/figures/compressor-pareto.svg",
    comparison = "inst/benchmarks/vcf-tabix.csv"
  )
}

#' Read the benchmark table shipped with CompreSSoR
#'
#' @param kind Either "pareto" for the first Pareto benchmark, "vcf"
#'   for the VCF.bgz plus Tabix comparison, "finngen" for the real
#'   non-EBI end-to-end conversion, or "finngen_optimization" for the
#'   reference-cache compression-time benchmark, "modes" for the
#'   chromosome/QC/panel edge benchmark, "release_gate" for the repeated
#'   real-GWAS release-gate run, or "release_gate_followup" for its real
#'   chromosome-parallel slice follow-up, or "storage_size" for the measured
#'   on-disk size comparison on the real release-gate GWAS.
#' @return A data.frame containing the selected measured benchmark.
#' @export
benchmark_table <- function(kind = c("pareto", "vcf", "finngen", "finngen_optimization",
                                    "modes", "release_gate", "release_gate_followup",
                                    "storage_size")) {
  kind <- match.arg(kind)
  filename <- switch(kind,
                     pareto = "first-pareto.csv",
                     vcf = "vcf-tabix.csv",
                     finngen = "finngen-end-to-end.csv",
                     finngen_optimization = "finngen-end-to-end-optimization.csv",
                     modes = "modes-edge.csv",
                     release_gate = "release-gate-benchmark.csv",
                     release_gate_followup = "release-gate-followup.csv",
                     storage_size = "storage-size-benchmark.csv")
  path <- system.file("benchmarks", filename, package = "CompreSSoR")
  if (!nzchar(path)) path <- file.path("inst", "benchmarks", filename)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

vcf_benchmark_metadata <- function() {
  list(
    benchmark_id = "vcf_tabix_11m",
    description = "Five-run VCF.bgz plus Tabix full-read and regional-access comparison",
    rows = 11106737L,
    source_build = "GRCh37",
    data = "inst/benchmarks/vcf-tabix.csv",
    runs = 5L
  )
}

finngen_benchmark_metadata <- function() {
  list(
    benchmark_id = "finngen_r2_antidepressants_end_to_end",
    description = "Real non-EBI FinnGen R2 GWAS converted to BP-spine Parquet",
    source = "https://storage.googleapis.com/finngen-public-data-r2/summary_stats/finngen_r2_ANTIDEPRESSANTS.gz",
    input_rows = 16111549L,
    output_rows = 3294606L,
    elapsed_seconds = 68.679,
    output_bytes = 132687271L,
    data = "inst/benchmarks/finngen-end-to-end.csv"
  )
}

finngen_optimization_metadata <- function() {
  list(
    benchmark_id = "finngen_r2_antidepressants_reference_cache_optimization",
    description = "Five warm-cache and one cold-cache end-to-end compression timings on the real FinnGen R2 GWAS",
    source = "https://storage.googleapis.com/finngen-public-data-r2/summary_stats/finngen_r2_ANTIDEPRESSANTS.gz",
    input_rows = 16111549L,
    baseline_elapsed_seconds = 68.679,
    cold_cache_elapsed_seconds = 77.370,
    warm_cache_runs = 5L,
    warm_cache_median_seconds = 63.820,
    speedup_vs_baseline = 1.076,
    output_rows = 3294606L,
    output_bytes = 132687915L,
    normalized_reference_cache_bytes = 142797914L,
    data = "inst/benchmarks/finngen-end-to-end-optimization.csv"
  )
}

mode_benchmark_metadata <- function() {
  list(
    benchmark_id = "compressor_modes_edge_550k",
    description = "Three-run 550,000-row balanced 22-chromosome QC and panel-mode benchmark",
    input_rows = 550000L,
    chromosomes = 22L,
    runs = 3L,
    data = "inst/benchmarks/modes-edge.csv",
    qc_serial_median_seconds = 0.374,
    qc_chrom_parallel_median_seconds = 0.956,
    core_chrom_parallel_median_seconds = 1.130
  )
}

release_gate_benchmark_metadata <- function() {
  list(
    benchmark_id = "compressor_release_gate_finngen_16111549",
    description = "Repeated real-GWAS conversion, decompression, regional and QC release-gate benchmark on the Mac mini",
    source = "https://storage.googleapis.com/finngen-public-data-r2/summary_stats/finngen_r2_ANTIDEPRESSANTS.gz",
    input_rows = 16111549L,
    convert_compress_median_seconds = 46.006,
    convert_decompress_median_seconds = 35.577,
    qc_compress_median_seconds = 96.003,
    qc_decompress_median_seconds = 31.353,
    convert_output_bytes = 892353845L,
    qc_output_bytes = 899137039L,
    data = "inst/benchmarks/release-gate-benchmark.csv"
  )
}

release_gate_followup_metadata <- function() {
  list(
    benchmark_id = "compressor_release_gate_slice_parallel_100k",
    description = "Repeated real-GWAS slice comparison of serial versus four-worker QC and regional access",
    input_rows = 100000L,
    qc_serial_median_seconds = 21.244,
    qc_chrom4_median_seconds = 18.985,
    direct_region_median_seconds = 0.048,
    q8_region_median_seconds = 0.018,
    data = "inst/benchmarks/release-gate-followup.csv"
  )
}

storage_size_benchmark_metadata <- function() {
  list(
    benchmark_id = "compressor_storage_size_finngen_16111549",
    description = "Measured on-disk storage comparison for a real 16.1-million-row GWAS",
    input_rows = 16111549L,
    standard_output_bytes = 892353845L,
    qc_output_bytes = 899137039L,
    data = "inst/benchmarks/storage-size-benchmark.csv"
  )
}
