# Current benchmark records

Only the dated 10-million-row FinnGen BluePebble records in
`finngen-10m-bp-20260804/` are current. Everything else is in
`archive/legacy-20260804/` and is retained for provenance only.

The current records contain:

- `format-screen/whole-file-10m-summary.csv`: five-run storage/read summary
  for the keyed format comparison;
- `format-screen/whole-file-10m-frontier.csv`: its measured frontier;
- `pcodec-thread-compare/`: five-run one-thread versus four-thread native
  Pcodec comparison;
- `finngen_10m_metadata.json`: input source, columns, filters, and checksums.

The current Pcodec result is 38,060,022 bytes and 0.451 seconds for the
four-thread whole-file read, against 135,388,867 bytes and 4.955 seconds for
TSV.gz. The Pcodec store includes the complete position plus REF/ALT identity
key.

The archive contains the earlier 1.1m-row chr1 screens, design sweeps,
reference-spine experiments, MR/access benchmarks, Pan-UKB scaling attempt,
and the numeric-only SE6/native-C++ projection experiments. None of those
should be used as the current headline result.
