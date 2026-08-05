# Core and HM3 panel provenance

The panel-backed modes use frozen GRCh38 identity dictionaries. The original
artifacts used for the current Mac mini benchmarks are:

| Panel | Source file | Data rows | SHA-256 |
|---|---|---:|---|
| core | `/Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz` | 6,271,258 | `167006fc385289160b3ab5b6dc6fede695d899df1f8decf415744da3644e7b4d` |
| HM3 | `/Volumes/crucial_x9/mr_atlas/data/panels/hm3/hm3_grch38_canonical_from_bfile.tsv.gz` | 1,212,965 | `33d158bd879b19c07564b60d968b792f2988909ca44b1d8c2e3ff4dedd33aea4` |

The core artifact is the ready `1kg_all_tag_r2_095_shared_keep_hm3` panel
used by the existing smoke/release benchmarks. The HM3 artifact is the
canonical GRCh38 panel produced from the HapMap3 source list and 1KG reference
mapping; its source/provenance details are recorded in the sibling
`hm3_panel_manifest.json` on the benchmark volume.

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
