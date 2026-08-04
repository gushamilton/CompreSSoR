# Final BP factor screen and lock-in

This record contains the five isolated optimization lanes run on BluePebble
before the selected-only final rerun. All lanes used the real FinnGen chr1
SNV store where applicable and were submitted through Slurm compute jobs. The
shared `opt2-LANE` setup attempt was cancelled and is not evidence.

## Locked standard

- native standalone Pcodec, with the self-contained position plus directed
  REF→ALT identity key;
- portable Rust build by default; `COMPRESSOR_RUSTFLAGS` remains an explicit
  site-specific override;
- 131,072-row key blocks and Pcodec pages, 65,536-row value blocks, and the
  existing stream-major layout;
- four worker processes by default for whole-store reads and one for regional
  or canonical-key reads; explicit `threads=` remains available;
- allocation-based R decoding remains the public path.

## Factor decisions

| Lane | BP evidence | Decision |
|---|---|---|
| Thread policy | Full keyed public probe: 0.207, 0.197, 0.200, 0.208 s for 1, 2, 4, 8 threads; all checksums matched. The independent native screen also found full keyed access essentially flat and selective access improved with more workers only for sparse blocks. | Keep four for whole-file reads as a stable default; keep one for selective reads. |
| CPU/build flags | Portable, native, AVX2/BMI, and Skylake-AVX512 builds all passed checksum checks; target-specific flags were tied or slower. | Lock portable builds; preserve an override rather than targeting one CPU. |
| Direct destination buffers | On actual stored streams, reused direct buffers were only about 0.5% faster than the current path; the larger synthetic re-encoded improvement did not transfer to the real store. | Do not complicate the public reader. |
| Fused native reader | Numeric-only fused probe passed, but the keyed checksum was `NaN`. | Exclude from the production path until complete keyed validation passes. |
| I/O geometry | 65,536, 131,072, and 262,144-row geometries differed by less than 1% in size/access in the valid stream-major screen; 65,536 was marginally best for one region probe. | Keep the already validated 131,072/65,536 compromise. |

## Selected-only rerun

The locked Pcodec candidate was rerun five times on the same 1,124,344-row
FinnGen chr1 source. Every run passed identity and numeric round-trip checks.
The median store size was **4,298,204 bytes** versus **15,186,281 bytes** for
the source TSV.gz, or **3.533× compression**. The strict write-then-read
benchmark median was **11.491 s** to write and **0.507 s** to read.

The five factor plots are retained as evidence, while only
`compressor-pareto-summary.png` is the README headline figure.
