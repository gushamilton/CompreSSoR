# Pcodec backend

The R package invokes `compressor_pcodec.py` when
`compress_sumstats(..., backend = "pcodec")` is used.

Install the small optional runtime in an isolated Python environment:

```bash
python3 -m venv .venv-compressor
.venv-compressor/bin/python -m pip install numpy pcodec zstandard
```

Then point R at that interpreter with either:

```r
options(CompreSSoR.python = "/absolute/path/.venv-compressor/bin/python")
```

or the `COMPRESSOR_PYTHON` environment variable. The output is a directory
containing a manifest, block indexes, Pcodec numerical streams, and
value-block-aligned Zstandard exception frames. The identity key is
self-contained and does not require an rsID table or external reference to
decode. Full reads hand compact codes to the package's compiled bridge reader;
regional and sparse reads decode only intersecting frames.
