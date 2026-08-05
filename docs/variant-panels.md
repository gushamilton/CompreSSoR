# Core, HM3, and core-plus panel provenance

Panel selection is a post-preparation operation. The package first imports and
validates the prepared GWAS; only then does it apply a panel or build the
core-plus union. The four explicit selection scopes are
`full`, `core`, `hm3`, and `core_plus`.

Every canonical panel row is the directed, allele-aware key
`chromosome:position:REF:ALT`, with chromosome labels and allele case
normalised. Coordinate-only or rsID-only legacy files may still be read for
interoperability, but they are recorded as legacy identity and do not acquire
HM3 provenance merely because they were passed as `variant_set`.

The panel-backed modes use frozen GRCh38 identity dictionaries. The original
artifacts used for the current Mac mini benchmarks are:

| Panel | Source file | Data rows | SHA-256 |
|---|---|---:|---|
| core | `/Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz` | 6,271,258 | `167006fc385289160b3ab5b6dc6fede695d899df1f8decf415744da3644e7b4d` |
| HM3 | `/Volumes/crucial_x9/mr_atlas/data/panels/hm3/hm3_grch38_canonical_from_bfile.tsv.gz` | 1,212,965 | `33d158bd879b19c07564b60d968b792f2988909ca44b1d8c2e3ff4dedd33aea4` |

The committed bundled artifact is one combined native Pcodec store:

```text
inst/extdata/panels/core_hm3.cpr/
├── manifest.json
├── manifest.sha256
├── native.index.json
├── position.pco
├── substitution.pco
└── hm3.pco
```

It contains 6,271,256 rows, with core membership implicit by row membership
and one lossless `hm3` flag. There are 1,199,729 HM3-positive rows in the core
universe; 13,236 valid HM3 source rows are outside the frozen core and are not
added. The committed artifact is 11,081,115 bytes. Identity is exact and
build-aware (`chromosome:position:REF:ALT`, stored as the same global-position
and directed REF→ALT streams used by GWAS stores); the HM3 flag is lossless.
The panel is not a FASTA, reference genome, harmonisation table, or liftover
chain.

The core source contributes 6,271,256 valid biallelic A/C/G/T SNVs after two
source rows are excluded. The HM3 source contributes 1,212,965 canonical
GRCh38 rows. Source and prepared hashes are recorded in the bundled manifest;
the source files remain on the Mac mini and are not committed.

These are treated as frozen upstream panel inputs. The core file is not
reconstructed by taking a simple union of the related base tag panel and HM3:
that does not reproduce the original selection. Reproduction in this
repository therefore means deterministic staging and normalization of the
identified source files, with source and output hashes recorded, rather than
claiming to recreate an upstream panel-selection algorithm that is not present
here.

Run the staging script on the machine holding the source files:

```sh
Rscript scripts/prepare_core_hm3_panels.R \
  --core-source /Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz \
  --hm3-source /Volumes/crucial_x9/mr_atlas/data/panels/hm3/hm3_grch38_canonical_from_bfile.tsv.gz \
  --output-dir /Volumes/crucial_x9/CompreSSoR-benchmarks/prepared-core-hm3-panels
```

The script reads only identity columns from the source TSVs, validates
GRCh38 primary-chromosome biallelic SNVs, deduplicates canonical keys in source
order, and writes `core.tsv.gz`, `hm3.tsv.gz`, chromosome-specific
`core_by_chrom/` and `hm3_by_chrom/` shards, and `panel_manifest.json`. In
the frozen core source, two indel rows are excluded by that SNV contract, so
the staged core dictionary has 6,271,256 rows; the exclusion is recorded in
the manifest rather than silently treated as a complete source reproduction.
Neither the raw sources nor the generated panel dictionaries belong in
`inputs/` or in the synced project.

## Named and hash-pinned resolution

Named `core` and `hm3` panels resolve from `COMPRESSOR_CORE_VARIANTS` and
`COMPRESSOR_HM3_VARIANTS` (or the equivalent
`COMPRESSOR_VARIANT_SET_CORE/HM3` variables). A named core/HM3 panel must have
a SHA-256 pin, supplied by `COMPRESSOR_CORE_VARIANTS_SHA256` or
`COMPRESSOR_HM3_VARIANTS_SHA256`, by the matching `COMPRESSOR_VARIANT_SET_*`
or `COMPRESSOR_PANEL_*` variable, or inline as `core@sha256:<digest>` /
`hm3@sha256:<digest>`. A staged `panel_manifest.json` can provide the same
pin through `prepared_sha256`. Resolution verifies the bytes before the panel
is used and records the requested name, observed/expected hashes, build,
source hash, local path, cache status, and identity class in selection
provenance.

When no environment variable or explicit panel path is supplied, the names
`core` and `hm3` resolve to the bundled `core_hm3.cpr` store. The bundled
manifest is protected by its payload and manifest checksums, so it does not
require an external hash environment variable. User-provided named panels
retain the existing hash-pin requirement.

Remote named assets are cached under `COMPRESSOR_VARIANT_SET_CACHE` (or
`COMPRESSOR_PANEL_CACHE`) using the panel name and pinned digest. The cache is
content-verified on reuse; a failed verification is an error, not a silent
refresh. Ordinary local panel files are never copied into the repository.

Compressed TSV.GZ panels use `data.table::fread()` with a gzip command fallback
and do not require the optional `R.utils` package. The reader probes identity
columns first, so large canonical dictionaries do not materialise unused
annotations.

At storage time, a large canonical panel is read as its single normalized
allele-aware key column. The reader deliberately does not expand those keys
back into redundant chromosome/position/allele columns before matching. The
membership operation then uses R/data.table's native hash lookup; splitting a
6-million-row panel into chromosome worker processes was slower and used much
more memory in the FinnGen timing. When a chromosome-shard directory is
provided, only shards represented in the harmonised input are read; this is
the preferred path for chromosome-local jobs. If `pigz` is available,
compressed readers use it automatically; otherwise they retain the `gzip`
fallback.

`read_variant_panel()` supports full reads plus chromosome, one-based row, and
canonical-key access. Row/key/chromosome requests read only the relevant
Pcodec key/flag blocks; full membership operations naturally read the full
core universe. `validate_variant_panel(..., full = TRUE)` verifies the manifest,
index partitions, every payload checksum, sorted identity, valid substitution
codes, and binary HM3 flags.

## Selection and core-plus threshold semantics

The internal selection helpers `select_full_variants()`,
`select_core_variants()`, `select_hm3_variants()`, and
`select_core_plus_variants()` return the same stable result: selected `data`,
the original-row logical `keep` mask, the resolved `panel`, a `regions` table,
and compact `metadata`. The API layer can use this result without depending on
the writer implementation.

Core-plus is the union of the selected core panel and windows around significant
variants. Significance is evaluated from the exact, post-harmonisation,
pre-encoding Z statistic using the strict two-sided rule
`2 * pnorm(-abs(z)) < pvalue_threshold`. It is never evaluated from a decoded
lossy store value. The selection metadata records the source, derivation,
operator, threshold, and `encoding = "not_encoded"`; the same information is
carried into the store manifest.
