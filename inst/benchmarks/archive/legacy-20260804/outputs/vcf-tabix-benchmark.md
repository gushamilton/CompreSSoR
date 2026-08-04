# VCF.bgz + Tabix comparison

This is a separate measured baseline against an indexed VCF representation of
11,106,737 real sumstats rows from the Mac mini. The source was GRCh37, so this
comparison is about storage and access rather than GRCh38 harmonisation.

The VCF stored REF/ALT plus beta, standard error and p-value in INFO fields,
was sorted with bcftools, compressed with bgzip, and indexed with Tabix. The
whole-file scan used bcftools view -H followed by an awk checksum; the
regional scan used Tabix for 1:100000000-110000000 followed by the same
checksum calculation.

| Representation | Rows | VCF.bgz | TBI index | Total | Ratio vs source gzip | Full median (s) | Region median (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| VCF.bgz + Tabix | 11,106,737 | 161,513,759 B | 1,519,452 B | 163,033,211 B | 0.982x | 13.28 | 0.09 |

Full-read runs were 12.08, 13.28, 14.05, 13.65 and 12.92 seconds. Regional
runs were 0.05, 0.05, 0.09, 0.09 and 0.09 seconds, returning 40,822 rows.

This is a useful contrast with CompreSSoR: VCF.bgz + Tabix has excellent
coordinate-indexed regional access, but the textual VCF representation is
larger than the source gzip once the index is included and is much slower for
a whole-genome scan. The source and first Pareto benchmark use different row
counts and readers, so this result is reported separately rather than silently
added to the original Pareto frontier.
