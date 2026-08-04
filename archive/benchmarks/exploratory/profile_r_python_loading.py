#!/usr/bin/env python3
"""Profile production Pcodec loading before data cross the R boundary.

This script deliberately reports OS-warm component timings. It is a CPU and
bridge profiler, not the cold-cache release benchmark.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import importlib.util
import json
import math
import statistics
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "inst" / "python" / "compressor_pcodec.py"
SPEC = importlib.util.spec_from_file_location("compressor_pcodec", BACKEND)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def timer() -> float:
    return time.perf_counter()


def direct_numpy_load(manifest: dict, context) -> tuple[dict, dict]:
    """Decode the public full-load columns into NumPy arrays in memory."""
    np, _, _, _ = MODULE._dependencies()
    marks: dict[str, float] = {}
    start = timer()

    positions = context.stream("position").read_rows().astype(np.uint64)
    substitutions = context.stream("substitution").read_rows().astype(np.uint64)
    marks["identity_stream_decode_seconds"] = timer() - start
    transform_start = timer()
    keys = MODULE._combine_wrapped_keys(positions, substitutions)
    chromosome, local_position, effect, other = MODULE._keys_to_columns(keys, manifest)
    marks["identity_transform_seconds"] = timer() - transform_start
    marks["identity_decode_seconds"] = (
        marks["identity_stream_decode_seconds"]
        + marks["identity_transform_seconds"]
    )

    code_start = timer()
    codes = {
        name: context.stream(name).read_rows()
        for name in ("z", "se", "eaf")
    }
    exceptions = context.exceptions()
    marks["numeric_code_decode_seconds"] = timer() - code_start

    semantic_start = timer()
    rows = np.arange(int(manifest["rows"]), dtype=np.int64)
    values = MODULE._decode_wrapped_semantics(codes, rows, manifest, exceptions)
    beta = values["z"] * values["se"]
    marks["semantic_decode_seconds"] = timer() - semantic_start

    touch_start = timer()
    checksum = (
        float(np.sum(chromosome, dtype=np.float64))
        + float(np.sum(local_position, dtype=np.float64))
        + float(np.sum(effect, dtype=np.float64))
        + float(np.sum(other, dtype=np.float64))
        + float(np.nansum(values["z"], dtype=np.float64))
        + float(np.nansum(values["se"], dtype=np.float64))
        + float(np.nansum(values["eaf"], dtype=np.float64))
        + float(np.nansum(beta, dtype=np.float64))
    )
    marks["touch_seconds"] = timer() - touch_start
    marks["load_seconds"] = sum(
        marks[name]
        for name in (
            "identity_decode_seconds",
            "numeric_code_decode_seconds",
            "semantic_decode_seconds",
        )
    )
    marks["load_and_touch_seconds"] = timer() - start
    arrays = {
        "chromosome": chromosome,
        "base_pair_location": local_position,
        "effect_allele": effect,
        "other_allele": other,
        "z": values["z"],
        "beta": beta,
        "standard_error": values["se"],
        "effect_allele_frequency": values["eaf"],
    }
    return marks, {"checksum": checksum, "arrays": arrays}


def summarise(records: list[dict]) -> dict:
    fields = [name for name in records[0] if name.endswith("_seconds")]
    return {
        field: {
            "median": statistics.median(float(row[field]) for row in records),
            "min": min(float(row[field]) for row in records),
            "max": max(float(row[field]) for row in records),
        }
        for field in fields
    }


def profile_stream_decode(store: Path, manifest: dict, reps: int) -> dict:
    """Compare sequential and threaded decode of the five physical streams."""
    np, _, _, _ = MODULE._dependencies()
    names = ("position", "substitution", "z", "se", "eaf")
    records = []
    for rep in range(1, reps + 1):
        for mode in ("sequential", "threaded"):
            context = MODULE.WrappedStoreReader(store, manifest)
            try:
                readers = {name: context.stream(name) for name in names}
                start = timer()
                if mode == "sequential":
                    arrays = {name: readers[name].read_rows() for name in names}
                else:
                    with concurrent.futures.ThreadPoolExecutor(
                        max_workers=len(names)
                    ) as executor:
                        values = list(executor.map(
                            lambda name: readers[name].read_rows(), names
                        ))
                    arrays = dict(zip(names, values))
                elapsed = timer() - start
                checksum = sum(
                    float(np.sum(array, dtype=np.float64))
                    for array in arrays.values()
                )
                records.append({
                    "rep": rep,
                    "mode": mode,
                    "stream_decode_seconds": elapsed,
                    "checksum": checksum,
                })
            finally:
                context.close()
    checksums = {row["checksum"] for row in records}
    if len(checksums) != 1:
        raise RuntimeError("stream decode checksums changed across runs")
    summary = {}
    for mode in ("sequential", "threaded"):
        values = [
            row["stream_decode_seconds"] for row in records
            if row["mode"] == mode
        ]
        summary[mode] = {
            "median_seconds": statistics.median(values),
            "min_seconds": min(values),
            "max_seconds": max(values),
        }
    return {"summary": summary, "runs": records}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.reps < 5:
        raise ValueError("at least five repetitions are required")

    store = args.store.expanduser().resolve(strict=True)
    manifest = MODULE.load_manifest(store)
    if manifest["format_version"] != MODULE.VERSION:
        raise ValueError("this profiler requires the wrapped v0.3 store")
    columns = [
        "chromosome", "base_pair_location", "effect_allele", "other_allele",
        "z", "beta", "standard_error", "effect_allele_frequency",
    ]

    direct_records = []
    bridge_records = []
    checksums = []
    with tempfile.TemporaryDirectory(prefix="compressor-loading-profile-") as temporary:
        temporary_path = Path(temporary)
        for rep in range(1, args.reps + 1):
            open_start = timer()
            context = MODULE.WrappedStoreReader(store, manifest)
            open_seconds = timer() - open_start
            try:
                marks, result = direct_numpy_load(manifest, context)
                marks.update({"rep": rep, "context_open_seconds": open_seconds})
                direct_records.append(marks)
                checksums.append(result["checksum"])
            finally:
                context.close()

            bridge = temporary_path / f"bridge-{rep}"
            open_start = timer()
            context = MODULE.WrappedStoreReader(store, manifest)
            open_seconds = timer() - open_start
            try:
                bridge_start = timer()
                MODULE.read_store(
                    store,
                    bridge,
                    columns=columns,
                    output_format="binary",
                    _manifest=manifest,
                    _key_index=context,
                    _value_index=context,
                )
                bridge_seconds = timer() - bridge_start
                bridge_bytes = sum(
                    path.stat().st_size for path in bridge.iterdir() if path.is_file()
                )
                bridge_records.append({
                    "rep": rep,
                    "context_open_seconds": open_seconds,
                    "decode_and_bridge_write_seconds": bridge_seconds,
                    "bridge_bytes": bridge_bytes,
                })
            finally:
                context.close()

    first_checksum = checksums[0]
    if not all(math.isclose(value, first_checksum, rel_tol=0, abs_tol=1e-6)
               for value in checksums[1:]):
        raise RuntimeError("direct NumPy checksums changed across runs")
    result = {
        "schema_version": "1.0.0",
        "benchmark_scope": "OS-warm component profile; not a cold-cache result",
        "store": str(store),
        "rows": int(manifest["rows"]),
        "columns": columns,
        "repetitions": args.reps,
        "direct_numpy": {
            "summary": summarise(direct_records),
            "runs": direct_records,
            "checksum": first_checksum,
        },
        "python_binary_bridge": {
            "summary": summarise(bridge_records),
            "runs": bridge_records,
        },
        "physical_stream_decode": profile_stream_decode(
            store, manifest, args.reps
        ),
    }
    rendered = json.dumps(result, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
