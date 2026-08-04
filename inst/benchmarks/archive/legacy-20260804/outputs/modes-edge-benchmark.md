# CompreSSoR conversion/QC mode benchmark

This is a three-run mode and edge-case benchmark on 550,000 balanced
synthetic rows across 22 chromosomes.

| Scenario | Median seconds | Input rows | Output rows |
|---|---:|---:|---:|
| Serial QC | 0.374 | 550,000 | 550,000 |
| Four-worker chromosome-parallel QC | 0.956 | 550,000 | 550,000 |
| Four-worker chromosome-parallel core filtering | 1.130 | 550,000 | 275,000 |

Default QC preserves all rows, while `core` retains only the explicit panel.
The tests also cover direct, reverse, complement, palindromic, incompatible,
unmatched and duplicated variants. Parallel output was checked against the
serial result. The small balanced fixture makes parallel startup visible;
the current in-memory chromosome path should not yet be treated as a speed
optimization for the full 16.1-million-row workload.
