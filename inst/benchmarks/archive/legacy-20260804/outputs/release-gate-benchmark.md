# CompreSSoR release-gate benchmark

Source: /Volumes/crucial_x9/CompreSSoR-benchmarks/raw/finngen_r2_ANTIDEPRESSANTS.gz
Reference: /Volumes/crucial_x9/mr_atlas/data/panels/1kg_all_tag_r2_095_shared_keep_hm3/variant_dictionary.shared.tsv.gz
External working root: /Volumes/crucial_x9/CompreSSoR-benchmarks

This run covers exact and lossy round trips, direct and q8 regional reads,
repeated conversion/decompression, and real FinnGen conversion plus QC.
Large inputs and stores remain outside the synced project.

```text
                              scenario replicate elapsed_seconds     rows
1         slice_convert_exact_compress         1           0.367   100000
2         slice_convert_exact_compress         2           0.181   100000
3         slice_convert_exact_compress         3           0.189   100000
4         slice_convert_exact_compress         4           0.189   100000
5         slice_convert_exact_compress         5           0.189   100000
6               slice_exact_decompress         1           0.291   100000
7          slice_exact_roundtrip_check         1           0.000   100000
8               slice_exact_decompress         2           0.192   100000
9          slice_exact_roundtrip_check         2           0.000   100000
10              slice_exact_decompress         3           0.192   100000
11         slice_exact_roundtrip_check         3           0.000   100000
12              slice_exact_decompress         4           0.192   100000
13         slice_exact_roundtrip_check         4           0.000   100000
14              slice_exact_decompress         5           0.228   100000
15         slice_exact_roundtrip_check         5           0.000   100000
16             slice_standard_compress         1           0.186   100000
17             slice_standard_compress         2           0.192   100000
18             slice_standard_compress         3           0.196   100000
19             slice_standard_compress         4           0.191   100000
20             slice_standard_compress         5           0.187   100000
21           slice_standard_decompress         1           0.200   100000
22      slice_standard_roundtrip_check         1           0.000   100000
23           slice_standard_decompress         2           0.201   100000
24      slice_standard_roundtrip_check         2           0.000   100000
25           slice_standard_decompress         3           0.200   100000
26      slice_standard_roundtrip_check         3           0.000   100000
27           slice_standard_decompress         4           0.199   100000
28      slice_standard_roundtrip_check         4           0.000   100000
29           slice_standard_decompress         5           0.198   100000
30      slice_standard_roundtrip_check         5           0.000   100000
31                 slice_region_direct         1           0.049       NA
32               slice_region_q8_cache         1           0.017       NA
33                 slice_region_direct         2           0.046       NA
34               slice_region_q8_cache         2           0.017       NA
35                 slice_region_direct         3           0.046       NA
36               slice_region_q8_cache         3           0.017       NA
37                 slice_region_direct         4           0.047       NA
38               slice_region_q8_cache         4           0.017       NA
39                 slice_region_direct         5           0.046       NA
40               slice_region_q8_cache         5           0.017       NA
41   finngen_convert_standard_compress         1          49.501       NA
42   finngen_convert_standard_compress         2          45.663       NA
43   finngen_convert_standard_compress         3          46.006       NA
44   finngen_convert_standard_compress         4          46.512       NA
45   finngen_convert_standard_compress         5          45.459       NA
46 finngen_convert_standard_decompress         1          40.985       NA
47     finngen_convert_roundtrip_shape         1           0.000 16111549
48 finngen_convert_standard_decompress         2          35.706       NA
49     finngen_convert_roundtrip_shape         2           0.000 16111549
50 finngen_convert_standard_decompress         3          33.434       NA
51     finngen_convert_roundtrip_shape         3           0.000 16111549
52 finngen_convert_standard_decompress         4          35.577       NA
53     finngen_convert_roundtrip_shape         4           0.000 16111549
54 finngen_convert_standard_decompress         5          34.565       NA
55     finngen_convert_roundtrip_shape         5           0.000 16111549
56          finngen_qc_serial_compress         1          99.472       NA
57          finngen_qc_serial_compress         2          90.605       NA
58          finngen_qc_serial_compress         3          96.003       NA
59               finngen_qc_decompress         1          39.401       NA
60          finngen_qc_roundtrip_shape         1           0.000 16111549
61               finngen_qc_decompress         2          31.353       NA
62          finngen_qc_roundtrip_shape         2           0.000 16111549
63               finngen_qc_decompress         3          30.075       NA
64          finngen_qc_roundtrip_shape         3           0.000 16111549
       bytes valid  max_error note
1    5774048  TRUE         NA     
2    5774048  TRUE         NA     
3    5774048  TRUE         NA     
4    5774048  TRUE         NA     
5    5774048  TRUE         NA     
6         NA    NA         NA     
7         NA  TRUE         NA     
8         NA    NA         NA     
9         NA  TRUE         NA     
10        NA    NA         NA     
11        NA  TRUE         NA     
12        NA    NA         NA     
13        NA  TRUE         NA     
14        NA    NA         NA     
15        NA  TRUE         NA     
16   5611039  TRUE         NA     
17   5611039  TRUE         NA     
18   5611039  TRUE         NA     
19   5611039  TRUE         NA     
20   5611039  TRUE         NA     
21        NA    NA         NA     
22        NA  TRUE 0.01137196     
23        NA    NA         NA     
24        NA  TRUE 0.01137196     
25        NA    NA         NA     
26        NA  TRUE 0.01137196     
27        NA    NA         NA     
28        NA  TRUE 0.01137196     
29        NA    NA         NA     
30        NA  TRUE 0.01137196     
31        NA    NA         NA     
32        NA    NA         NA     
33        NA    NA         NA     
34        NA    NA         NA     
35        NA    NA         NA     
36        NA    NA         NA     
37        NA    NA         NA     
38        NA    NA         NA     
39        NA    NA         NA     
40        NA    NA         NA     
41 892353845  TRUE         NA     
42 892353845  TRUE         NA     
43 892353845  TRUE         NA     
44 892353845  TRUE         NA     
45 892353845  TRUE         NA     
46 892353845    NA         NA     
47        NA  TRUE         NA     
48 892353845    NA         NA     
49        NA  TRUE         NA     
50 892353845    NA         NA     
51        NA  TRUE         NA     
52 892353845    NA         NA     
53        NA  TRUE         NA     
54 892353845    NA         NA     
55        NA  TRUE         NA     
56 899137039  TRUE         NA     
57 899137039  TRUE         NA     
58 899137039  TRUE         NA     
59 899137039    NA         NA     
60        NA  TRUE         NA     
61 899137039    NA         NA     
62        NA  TRUE         NA     
63 899137039    NA         NA     
64        NA  TRUE         NA     
```

R version: R version 4.5.2 (2025-10-31)
