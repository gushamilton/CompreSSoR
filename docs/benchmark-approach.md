# CompreSSoR benchmark approach

This is the benchmark contract for the final format comparison. The three
access contracts below are deliberately reported separately: they answer
different questions and must not be collapsed into one Pareto plot.

## Format matrix

Run the agreed keyed-format matrix, not the obsolete identity-free formats.
Every format must contain the same self-contained variant identity:

```text
global GRCh38 position + substitution code
```

The matrix is the final set of approximately forty candidates already
specified in the keyed benchmark. It includes the tested Parquet, fixed-width
binary, framed binary, semantic/binary, Pcodec, VCF+Tabix and TSV.gz controls,
with the same rows, row order and source fields. Do not add formats during the
run unless a format is needed as a clearly labelled control.

### Fixed benchmark input

The full-file benchmark uses exactly **10,000,000 real biallelic SNP rows**
from one non-EBI FinnGen summary-statistics GWAS. The normalized input is
materialized once on BP, hashed, and reused by every job. Indels and malformed
alleles are excluded during this one-time input preparation; no format is
allowed to apply its own row filter. The exact input path, row count, byte
count, SHA-256 and source-column mapping are recorded beside the results.

This 10-million-row input is used for all full-file and reconstruction
measurements. Regional and sparse workloads query that same store, so their
timings are not based on a smaller synthetic fixture.

### Immediate BP run

The current BP run is restricted to the whole-file compact/native access
contract across the complete keyed format matrix. Reconstruction and
end-to-end MR remain specified below but are deliberately deferred until this
baseline has completed.

For each format record:

- on-disk bytes, excluding any shared reference that is stored once;
- write time, reported separately and not used as the main access Pareto axis;
- five independent read runs;
- numerical and identity checksums/tolerances;
- host, CPU allocation, software versions and cache/process condition.

The plots use speed on the x-axis and compression ratio relative to the keyed
TSV.gz control on the y-axis. The shared reference is excluded from every
format equally.

## Benchmark 1: native analytical access

This is the primary CompreSSoR/FastMR contract. The reader returns only the
primitive values needed for analysis:

```text
compact_position, substitution, z, standard_error,
effect_allele_frequency
```

It must not construct chromosome strings, REF/ALT strings, beta, p-values or a
general-purpose reconstructed data frame. A requested-variant query is given
the compact identity keys and returns rows aligned to those keys.

Measure three workloads:

1. Full-GWAS native scan, with a deterministic reduction over every returned
   value so that all streams are genuinely decoded.
2. Regional access over the same fixed regions for every format.
3. Sparse access for a fixed set of 25 randomly selected variant keys.

Report first-access/fresh-process and repeated warm-process results separately.
For the sparse workload include file opening, block selection, key lookup,
numeric decode and alignment; preloading the entire GWAS is not allowed.

## Benchmark 2: standard reconstruction

This measures interoperability and the cost of rebuilding a conventional R
sum-statistics table. The logical output is:

```text
chromosome, base_pair_location, effect_allele, other_allele,
z, standard_error, effect_allele_frequency, beta, p_value
```

Beta and p-value are derived from the stored semantic values where the format
does not store them. The timing includes identity expansion, semantic decoding,
R object materialisation and column projection. It excludes package loading and
writing the reconstructed table unless a separate export result is reported.

Measure full-file and regional reconstruction. Compare this contract only with
other formats reconstructed to the same logical columns. Do not compare it to
a compact `position/substitution` Parquet read.

## Benchmark 3: end-to-end MR

This is the headline application benchmark. Use identical exposure/outcome
files, variant keys, instruments, MR method, numerical tolerances and output
checksum for every format.

The workloads are:

1. One exposure to one outcome, using 25 instruments.
2. A 5 x 5 exposure/outcome grid.
3. A 25 x 25 exposure/outcome grid.

Run each grid in two explicit modes:

- **independent access:** files are opened and queried for each MR comparison;
- **reused access:** each file is opened once and reused for all comparisons.

The primary timer includes file open, variant lookup, decode, alignment and MR
calculation. Report the component timings as diagnostics, but use total wall
time for the MR Pareto plot. Do not allow an already materialised R table or an
unreported OS-cache warm-up to make one format look artificially better.

## Execution rules

- Five fresh-process repeats per format, workload and access mode.
- Randomise format order between repeats.
- Keep all formats on the same machine, CPU allocation and filesystem.
- Load packages before starting the timer, but include normal file opening and
  decoding in the timer.
- Verify that every result contains the same variant keys and row count.
- Use explicit precision gates for z, SE, EAF and MR estimates.
- Report fresh-process and warm-process results rather than mixing them.
- Do not clear the system page cache from a login node; if a true cold-cache
  comparison is required, use isolated compute jobs and label it as such.
- Keep raw run files and the exact benchmark metadata; the README receives only
  the summary Pareto plots and short interpretation.

## Required outputs

Create one result table and one Pareto plot for each benchmark contract:

- `pareto-native-access.csv` and `pareto-native-access.png`;
- `pareto-reconstruction.csv` and `pareto-reconstruction.png`;
- `pareto-end-to-end-mr.csv` and `pareto-end-to-end-mr.png`.

The plot labels must use descriptive format names. The README should link only
the three summary plots; the benchmark directory should retain all per-run
records, checksums, environment metadata and the complete matrix specification.
