# P-value side-channel benchmark

The compact result files for the Mac mini benchmark are stored here. Raw
FinnGen input, generated stores, and run logs remain on the Mac mini under the
external benchmark workspace.

The benchmark compares a duplicate native Pcodec store of rows selected by
`p <= 1e-7`, `1e-6`, and `1e-5` with an aligned Pcodec `uint8` 0/1 significance
stream. The 100,000-row fixture is real prepared FinnGen data and includes the
source's low-p rows through `p <= 1e-5`. It also measures reading reconstructed
p from the baseline store and filtering at `p <= 1e-5`, plus a 10,000-row
designated-hit cardinality stress case. The 10,000 rows are not source-level
significance labels.

The full-file comparison is in
`pvalue-sidechannel-full-finngen-summary.md`. It uses 14,923,434 valid SNPs
from the full ANTIDEPRESSANTS source and a second real FinnGen chromosome-1
file from the Mac mini benchmark workspace. Thresholds are `1e-4` through
`1e-15`. The two headline paths identify rows with either
whole-store reconstructed-p filtering or aligned flag inspection, then fetch
the same MR projection: identity, beta, standard error, and effect-allele
frequency. The Neale CRP extension adds sparse tail checkpoints through
`1e-200`; its combined comparison plot is
`pvalue-sidechannel-mr-comparison.png`.
