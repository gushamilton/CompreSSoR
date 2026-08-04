#!/usr/bin/env python3
"""Profile canonical-key reads through the production Pcodec backend."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BACKEND = ROOT / "inst" / "python" / "compressor_pcodec.py"
SPEC = importlib.util.spec_from_file_location("compressor_pcodec", BACKEND)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def elapsed(callable_):
    start = time.perf_counter()
    callable_()
    return time.perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--store", type=Path, required=True)
    parser.add_argument("--keys", type=int, default=25)
    parser.add_argument("--logical-reads", type=int, default=10)
    parser.add_argument("--reps", type=int, default=5)
    args = parser.parse_args()

    manifest = MODULE.load_manifest(args.store)
    n = int(manifest["rows"])
    row_ids = sorted({round(i * (n - 1) / max(1, args.keys - 1)) for i in range(args.keys)})

    with tempfile.TemporaryDirectory(prefix="compressor-production-profile-") as temp:
        root = Path(temp)
        identity_path = root / "identity.tsv"
        MODULE.read_store(
            args.store,
            identity_path,
            row_indices=row_ids,
            columns=["chromosome", "base_pair_location", "effect_allele", "other_allele"],
        )
        with identity_path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        keys = [
            f'{row["chromosome"]}:{row["base_pair_location"]}:'
            f'{row["other_allele"]}:{row["effect_allele"]}'
            for row in rows
        ]
        key_path = root / "keys.txt"
        key_path.write_text("\n".join(keys) + "\n")

        key_lookup = []
        in_process = []
        cold_cli = []
        for rep in range(args.reps):
            def lookup():
                current = MODULE.load_manifest(args.store)
                index = json.loads(
                    (args.store / current["files"]["key_index"]).read_text()
                )
                targets = MODULE.parse_identity_keys(keys, current)
                found = MODULE._rows_for_identity_keys(
                    args.store, current, index, targets
                )
                if len(found) != len(keys):
                    raise RuntimeError("profile keys did not resolve exactly")

            key_lookup.append(elapsed(lookup))

            bridge = root / f"bridge-{rep}"
            in_process.append(elapsed(lambda: MODULE.read_store(
                args.store,
                bridge,
                identity_keys=keys,
                columns=["chromosome", "base_pair_location", "effect_allele",
                         "other_allele", "beta", "standard_error"],
                output_format="binary",
            )))

            cli_bridge = root / f"cli-bridge-{rep}"
            command = [
                sys.executable, str(BACKEND), "read", "--store", str(args.store),
                "--output", str(cli_bridge), "--keys-file", str(key_path),
                "--columns", "chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error",
                "--output-format", "binary",
            ]
            cold_cli.append(elapsed(lambda: subprocess.run(
                command, check=True, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )))

        batch_times = []
        for rep in range(args.reps):
            batch_root = root / f"batch-{rep}"
            batch_root.mkdir()

            def batch():
                for index in range(args.logical_reads):
                    MODULE.read_store(
                        args.store,
                        batch_root / str(index),
                        identity_keys=keys,
                        columns=["chromosome", "base_pair_location", "effect_allele",
                                 "other_allele", "beta", "standard_error"],
                        output_format="binary",
                    )

            batch_times.append(elapsed(batch))

    def summary(values):
        return {
            "median_seconds": statistics.median(values),
            "min_seconds": min(values),
            "runs": len(values),
        }

    result = {
        "store": str(args.store),
        "rows": n,
        "keys": len(keys),
        "logical_reads": args.logical_reads,
        "key_lookup_only": summary(key_lookup),
        "in_process_binary_read": summary(in_process),
        "cold_cli_binary_read": summary(cold_cli),
        "in_process_logical_batch": summary(batch_times),
        "in_process_batch_per_read_seconds": statistics.median(batch_times) / args.logical_reads,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
