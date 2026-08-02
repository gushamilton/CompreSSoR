# CompreSSoR benchmark record

first-pareto.csv is the underlying data for the first Pareto plot shipped
with the package at inst/figures/compressor-pareto.svg.

This is a measured five-run full-read benchmark on a real 10-million-row GWAS.
The timing metric is the median wall-clock time to read the complete dataset;
lower is faster. Compression ratio is relative to the source gzip file; higher
is smaller. The Pareto flag identifies the points retained by the plot's
speed/compression frontier.

These numbers are a design record, not a promise of fixed performance. They
depend on hardware, filesystem, software versions, data layout and the exact
read path. They support the package choice of portable q9/SE10 Parquet as the
standard durable store, with q8/framed representations treated as optional
serving caches.

## Conversion/QC modes

`modes-edge.csv` records three runs on a balanced 550,000-row, 22-chromosome
fixture. It checks that serial and chromosome-parallel QC agree, that the
default QC mode retains every row, and that explicit core-panel filtering is
accounted for. Median times are 0.374 seconds for serial QC, 0.956 seconds
for four-worker QC, and 1.130 seconds for four-worker core filtering. The
small fixture is a correctness/mode record; it is not a claim that the
current in-memory parallel path accelerates a full GWAS.

## Release-gate real-GWAS benchmark

`release-gate-benchmark.csv` records the Mac mini run on the real
16,111,549-row FinnGen R2 ANTIDEPRESSANTS GWAS. Five convert-only compression
runs had a 46.006-second median and five whole-store decompressions had a
35.577-second median. Three GRCh38 QC compression runs had a 96.003-second
median and three QC decompressions had a 31.353-second median. The preserve-all
stores were 892,353,845 and 899,137,039 bytes respectively, and all shape and
validation checks passed.

`release-gate-followup.csv` records five exact/standard round-trip checks,
five regional reads through the durable store and q8 cache, and three real
100,000-row serial versus four-worker QC comparisons. The direct-region and
q8-cache medians were 0.048 and 0.018 seconds; serial and four-worker QC
medians were 21.244 and 18.985 seconds, with real serial/parallel output
parity checked separately.

## Storage-size benchmark

`storage-size-benchmark.csv` records absolute on-disk bytes and bytes per row
for the same 16,111,549-row real GWAS. The source gzip was 391,708,712 bytes,
the raw TSV stream was 1,273,085,087 bytes, and the standard q9 Parquet store
was 892,353,845 bytes for convert-only output and 899,137,039 bytes after
GRCh38 QC. Exact Parquet and q8 cache sizes were measured in temporary stores;
the q8 cache is reported as an optional serving layer, not as a standalone
replacement for the q9 variant spine.

## VCF.bgz + Tabix comparison

vcf-tabix.csv and vcf-tabix-runs.csv record a separate comparison against
an indexed VCF representation of 11,106,737 real sumstats rows. The VCF kept
REF/ALT plus beta, SE and p-value in INFO fields, was bgzip-compressed,
bcftools-sorted and Tabix-indexed. It was measured separately because the
source contains 11.1 million rows and is GRCh37, whereas the first Pareto
benchmark is a 10-million-row design record and CompreSSoR targets GRCh38.
