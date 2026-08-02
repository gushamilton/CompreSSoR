# CompreSSoR release-gate benchmark

Source: /Volumes/crucial_x9/CompreSSoR-benchmarks/raw/finngen_r2_ANTIDEPRESSANTS.gz
Reference: /Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz
External working root: /Volumes/crucial_x9/CompreSSoR-benchmarks

This run covers exact and lossy round trips, direct and q8 regional reads,
repeated conversion/decompression, and real FinnGen conversion plus QC.
Large inputs and stores remain outside the synced project.

```text
                         scenario replicate elapsed_seconds   rows   bytes
1    slice_convert_exact_compress         1           0.507 100000 5774048
2    slice_convert_exact_compress         2           0.194 100000 5774048
3    slice_convert_exact_compress         3           0.189 100000 5774048
4    slice_convert_exact_compress         4           0.189 100000 5774048
5    slice_convert_exact_compress         5           0.187 100000 5774048
6          slice_exact_decompress         1           0.299 100000      NA
7     slice_exact_roundtrip_check         1           0.000 100000      NA
8          slice_exact_decompress         2           0.188 100000      NA
9     slice_exact_roundtrip_check         2           0.000 100000      NA
10         slice_exact_decompress         3           0.187 100000      NA
11    slice_exact_roundtrip_check         3           0.000 100000      NA
12         slice_exact_decompress         4           0.188 100000      NA
13    slice_exact_roundtrip_check         4           0.000 100000      NA
14         slice_exact_decompress         5           0.224 100000      NA
15    slice_exact_roundtrip_check         5           0.000 100000      NA
16        slice_standard_compress         1           0.182 100000 5611039
17        slice_standard_compress         2           0.182 100000 5611039
18        slice_standard_compress         3           0.182 100000 5611039
19        slice_standard_compress         4           0.186 100000 5611039
20        slice_standard_compress         5           0.188 100000 5611039
21      slice_standard_decompress         1           0.196 100000      NA
22 slice_standard_roundtrip_check         1           0.000 100000      NA
23      slice_standard_decompress         2           0.197 100000      NA
24 slice_standard_roundtrip_check         2           0.000 100000      NA
25      slice_standard_decompress         3           0.191 100000      NA
26 slice_standard_roundtrip_check         3           0.000 100000      NA
27      slice_standard_decompress         4           0.197 100000      NA
28 slice_standard_roundtrip_check         4           0.000 100000      NA
29      slice_standard_decompress         5           0.199 100000      NA
30 slice_standard_roundtrip_check         5           0.000 100000      NA
31            slice_region_direct         1           0.053     NA      NA
32          slice_region_q8_cache         1           0.018     NA      NA
33            slice_region_direct         2           0.048     NA      NA
34          slice_region_q8_cache         2           0.018     NA      NA
35            slice_region_direct         3           0.049     NA      NA
36          slice_region_q8_cache         3           0.018     NA      NA
37            slice_region_direct         4           0.046     NA      NA
38          slice_region_q8_cache         4           0.018     NA      NA
39            slice_region_direct         5           0.048     NA      NA
40          slice_region_q8_cache         5           0.017     NA      NA
41       slice_qc_serial_compress         1          24.185 100000 5656310
42       slice_qc_chrom4_compress         1          17.055 100000 5656310
43       slice_qc_serial_compress         2          17.184 100000 5656310
44       slice_qc_chrom4_compress         2          19.004 100000 5656310
45       slice_qc_serial_compress         3          21.244 100000 5656310
46       slice_qc_chrom4_compress         3          18.985 100000 5656310
   valid  max_error note
1   TRUE         NA     
2   TRUE         NA     
3   TRUE         NA     
4   TRUE         NA     
5   TRUE         NA     
6     NA         NA     
7   TRUE         NA     
8     NA         NA     
9   TRUE         NA     
10    NA         NA     
11  TRUE         NA     
12    NA         NA     
13  TRUE         NA     
14    NA         NA     
15  TRUE         NA     
16  TRUE         NA     
17  TRUE         NA     
18  TRUE         NA     
19  TRUE         NA     
20  TRUE         NA     
21    NA         NA     
22  TRUE 0.01137196     
23    NA         NA     
24  TRUE 0.01137196     
25    NA         NA     
26  TRUE 0.01137196     
27    NA         NA     
28  TRUE 0.01137196     
29    NA         NA     
30  TRUE 0.01137196     
31    NA         NA     
32    NA         NA     
33    NA         NA     
34    NA         NA     
35    NA         NA     
36    NA         NA     
37    NA         NA     
38    NA         NA     
39    NA         NA     
40    NA         NA     
41  TRUE         NA     
42  TRUE         NA     
43  TRUE         NA     
44  TRUE         NA     
45  TRUE         NA     
46  TRUE         NA     
```

R version: R version 4.5.2 (2025-10-31)
data.table threads: 5
