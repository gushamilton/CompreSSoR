#!/usr/bin/env python3
"""Sweep physical Pcodec frame geometry on a real production store.

The sweep reframes already-quantised production codes, so every candidate has
identical variant identity and numerical semantics. Large stores and temporary
bridges stay under the caller-supplied external output directory.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import random
import shutil
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


def directory_bytes(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def load_prefix(store: Path, rows: int):
    np, _, _, _ = MODULE._dependencies()
    manifest = MODULE.load_manifest(store)
    rows = min(rows, int(manifest["rows"]))
    key_index = json.loads((store / manifest["files"]["key_index"]).read_text())
    key_parts = []
    for frame in key_index:
        if int(frame["row_start"]) >= rows:
            break
        decoded = MODULE.decode_key_frame(store, manifest, frame)
        key_parts.append(decoded[: max(0, rows - int(frame["row_start"]))])
    keys = np.concatenate(key_parts)[:rows]

    value_index = json.loads((store / manifest["files"]["value_index"]).read_text())
    parts = {name: [] for name in ("z", "eaf", "se")}
    exception_parts = []
    for frame in value_index:
        if int(frame["row_start"]) >= rows:
            break
        decoded = MODULE.decode_value_frame(store, manifest, frame)
        take = max(0, rows - int(frame["row_start"]))
        for name in parts:
            parts[name].append(decoded[name][:take])
        exceptions = MODULE.decode_exception_frame(store, manifest, frame)
        exception_parts.append(exceptions[exceptions["row"] < rows])
    exception_dtype = np.dtype([
        ("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"),
        ("eaf", "<f4"), ("flags", "u1"),
    ])
    values = {
        name: np.concatenate(chunks)[:rows] for name, chunks in parts.items()
    }
    values["exceptions"] = (
        np.concatenate(exception_parts)
        if exception_parts else np.empty(0, dtype=exception_dtype)
    )
    center_rows = int(
        manifest["semantic_codec"].get(
            "se_center_block_rows", manifest["value_block_rows"]
        )
    )
    centers_needed = (rows + center_rows - 1) // center_rows
    centers = manifest["semantic_codec"]["block_centers_log2_residual"][:centers_needed]
    return manifest, keys, values, center_rows, centers


def write_reframed_store(
    source_manifest,
    keys,
    values,
    center_rows: int,
    centers,
    output: Path,
    key_rows: int,
    value_rows: int,
):
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    key_index, _ = MODULE.encode_key_frames(keys, output, block_rows=key_rows)
    value_index, _ = MODULE.encode_value_frames(values, output, block_rows=value_rows)
    manifest = copy.deepcopy(source_manifest)
    manifest["rows"] = len(keys)
    manifest["n_rows"] = len(keys)
    manifest["key_block_rows"] = key_rows
    manifest["value_block_rows"] = value_rows
    semantic = manifest["semantic_codec"]
    semantic["se_center_block_rows"] = center_rows
    semantic["block_centers_log2_residual"] = list(centers)
    semantic["exception_rows"] = int(len(values["exceptions"]))
    manifest["metadata"] = {
        **manifest.get("metadata", {}),
        "benchmark_geometry": {"key_block_rows": key_rows, "value_block_rows": value_rows},
    }
    (output / "key.index.json").write_text(json.dumps(key_index, indent=2) + "\n")
    (output / "value.index.json").write_text(json.dumps(value_index, indent=2) + "\n")
    integrity = {}
    for relative in manifest["files"].values():
        blob = (output / relative).read_bytes()
        integrity[relative] = {
            "bytes": len(blob),
            "sha256": hashlib.sha256(blob).hexdigest(),
        }
    manifest["integrity"] = {"algorithm": "sha256", "files": integrity}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    MODULE._seal_manifest(output)
    MODULE.load_manifest(output)
    return manifest


def canonical_keys(keys, manifest, row_ids):
    selected = keys[row_ids]
    chromosomes, positions, effects, refs = MODULE._keys_to_rows(selected, manifest)
    return [
        f"{chromosome}:{int(position)}:{ref}:{effect}"
        for chromosome, position, ref, effect in zip(
            chromosomes, positions, refs, effects
        )
    ]


def timed(reps, operation):
    operation(-1)
    values = []
    checksums = []
    for rep in range(reps):
        started = time.perf_counter()
        checksums.append(operation(rep))
        values.append(time.perf_counter() - started)
    if len(set(checksums)) != 1:
        raise RuntimeError("benchmark repetitions produced different outputs")
    return values, checksums[0]


def benchmark_store(store, manifest, key_sets, temp_root, reps):
    workloads = {
        "canonical_25": (key_sets["canonical_25"], [
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error",
        ]),
        "canonical_1000": (key_sets["canonical_1000"], [
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error",
        ]),
    }
    rows = []
    for workload, (keys, columns) in workloads.items():
        def operation(rep):
            output = temp_root / f"{store.name}-{workload}-{rep}"
            if output.exists():
                shutil.rmtree(output)
            MODULE.read_store(
                store, output, identity_keys=keys, columns=columns,
                output_format="binary",
            )
            bridge = json.loads((output / "bridge.json").read_text())
            result = (int(bridge["rows"]), tuple(sorted(bridge["files"])))
            shutil.rmtree(output)
            return result

        timings, checksum = timed(reps, operation)
        rows.append({
            "workload": workload,
            "median_seconds": statistics.median(timings),
            "min_seconds": min(timings),
            "max_seconds": max(timings),
            "runs": reps,
            "rows_returned": checksum[0],
        })

    for workload, columns in (
        ("region_chr1_1mb", ["chromosome", "base_pair_location", "beta", "standard_error"]),
        ("full_beta_se", ["beta", "standard_error"]),
    ):
        def operation(rep):
            output = temp_root / f"{store.name}-{workload}-{rep}"
            if output.exists():
                shutil.rmtree(output)
            kwargs = {}
            if workload.startswith("region"):
                kwargs = {"chromosome": "1", "start": 100_000_000, "end": 101_000_000}
            MODULE.read_store(
                store, output, columns=columns, output_format="binary", **kwargs
            )
            bridge = json.loads((output / "bridge.json").read_text())
            result = (int(bridge["rows"]), tuple(sorted(bridge["files"])))
            shutil.rmtree(output)
            return result

        timings, checksum = timed(reps, operation)
        rows.append({
            "workload": workload,
            "median_seconds": statistics.median(timings),
            "min_seconds": min(timings),
            "max_seconds": max(timings),
            "runs": reps,
            "rows_returned": checksum[0],
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-store", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--rows", type=int, default=2_000_000)
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument(
        "--key-block-rows", default="4096,8192,16384,32768,65536,131072"
    )
    parser.add_argument(
        "--value-block-rows", default="4096,8192,16384,32768,65536"
    )
    args = parser.parse_args()
    if args.reps < 5:
        raise ValueError("at least five repetitions are required")
    key_sizes = [int(value) for value in args.key_block_rows.split(",")]
    value_sizes = [int(value) for value in args.value_block_rows.split(",")]
    args.output_root.mkdir(parents=True, exist_ok=True)
    stores_root = args.output_root / "stores"
    stores_root.mkdir(exist_ok=True)
    temp_root = args.output_root / "tmp"
    temp_root.mkdir(exist_ok=True)

    source_manifest, keys, values, center_rows, centers = load_prefix(
        args.source_store, args.rows
    )
    n = len(keys)
    np, _, _, _ = MODULE._dependencies()
    rows25 = np.asarray(sorted({round(i * (n - 1) / 24) for i in range(25)}), dtype=np.int64)
    random_rows = sorted(random.Random(20260803).sample(range(n), 1000))
    key_sets = {
        "canonical_25": canonical_keys(keys, source_manifest, rows25),
        "canonical_1000": canonical_keys(
            keys, source_manifest, np.asarray(random_rows, dtype=np.int64)
        ),
    }

    records = []
    for key_rows in key_sizes:
        for value_rows in value_sizes:
            label = f"k{key_rows}-v{value_rows}"
            store = stores_root / label
            started = time.perf_counter()
            manifest = write_reframed_store(
                source_manifest, keys, values, center_rows, centers,
                store, key_rows, value_rows,
            )
            encode_seconds = time.perf_counter() - started
            size = directory_bytes(store)
            for result in benchmark_store(store, manifest, key_sets, temp_root, args.reps):
                records.append({
                    "store": label,
                    "rows": n,
                    "key_block_rows": key_rows,
                    "value_block_rows": value_rows,
                    "se_center_block_rows": center_rows,
                    "bytes": size,
                    "encode_seconds": encode_seconds,
                    **result,
                })
            output = args.output_root / "geometry-runs.csv"
            with output.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
            print(json.dumps(records[-4:], indent=2), flush=True)

    (args.output_root / "geometry-summary.json").write_text(json.dumps({
        "source_store": str(args.source_store),
        "rows": n,
        "reps": args.reps,
        "key_sizes": key_sizes,
        "value_sizes": value_sizes,
        "records": records,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
