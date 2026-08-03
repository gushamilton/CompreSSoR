#!/usr/bin/env python3
"""Small edge/stress test for the optional Pcodec backend.

Run this on the Mac mini in the configured Pcodec environment.  It exercises
the recursive key escapes, exact key round-trip, semantic exceptions, full
decode, row slicing, and coordinate-region access without touching project
large-data directories.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import math
import shutil
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "compressor_pcodec", ROOT / "inst" / "python" / "compressor_pcodec.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_input(path: Path) -> list[dict[str, str]]:
    rows = [
        ("1", 1, "A", "C", 0.10, 0.05, 0.10),
        ("1", 512, "G", "T", -0.20, 0.10, 0.20),
        ("1", 70_000, "C", "A", 4.0, 0.10, 0.30),
        ("1", 200_000, "T", "G", -0.40, 0.20, float("nan")),
        ("2", 1, "A", "G", 0.15, 0.10, 0.40),
        ("X", 100, "C", "T", -0.05, 0.05, 0.60),
    ]
    expected = []
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency",
        ])
        for chrom, pos, alt, ref, beta, se, eaf in rows:
            writer.writerow([chrom, pos, alt, ref, beta, se, eaf])
            expected.append({"chromosome": chrom, "position": pos, "alt": alt, "ref": ref})
    return expected


def key(row: dict[str, str]) -> tuple[int, int, int, int]:
    offsets = MODULE.chromosome_offsets()
    return (
        (offsets[row["chromosome"]] + int(row["position"]) - 1),
        MODULE.BASE_CODE[row["ref"]],
        MODULE.BASE_CODE[row["alt"]],
        0,
    )


def write_rows(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chromosome", "base_pair_location", "effect_allele", "other_allele",
            "beta", "standard_error", "effect_allele_frequency", "z",
        ])
        writer.writerows(rows)


def read_rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="compressor-pcodec-test-") as temp:
        root = Path(temp)
        input_path = root / "input.tsv.gz"
        output = root / "store"
        decoded_path = root / "decoded.tsv"
        expected = write_input(input_path)
        manifest = MODULE.write_store(input_path, output)
        assert manifest["backend"] == "pcodec"
        assert manifest["identity"]["external_reference_required"] is False

        key_index = json.loads((output / "key.index.json").read_text())
        decoded_keys = []
        for frame in key_index:
            decoded_keys.extend(MODULE.decode_key_frame(output, manifest, frame).tolist())
        source_keys, _ = MODULE.identity_keys(MODULE.read_tsv(input_path))
        assert decoded_keys == source_keys.tolist(), "identity key frame round-trip failed"

        MODULE.read_store(output, decoded_path)
        with decoded_path.open(encoding="utf-8", newline="") as handle:
            decoded = list(csv.DictReader(handle, delimiter="\t"))
        assert len(decoded) == len(expected)
        assert [
            (row["chromosome"], int(row["base_pair_location"]), row["effect_allele"], row["other_allele"])
            for row in decoded
        ] == [
            (row["chromosome"], row["position"], row["alt"], row["ref"])
            for row in sorted(expected, key=key)
        ]
        assert float(decoded[2]["z"]) > 3.5, "out-of-range z was not retained as an exception"
        assert decoded[3]["effect_allele_frequency"] == "nan", "missing EAF was not retained"

        region_path = root / "region.tsv"
        MODULE.read_store(output, region_path, chromosome="1", start=60_000, end=80_000)
        with region_path.open(encoding="utf-8", newline="") as handle:
            region = list(csv.DictReader(handle, delimiter="\t"))
        assert len(region) == 1 and region[0]["base_pair_location"] == "70000", region

        slice_path = root / "slice.tsv"
        MODULE.read_store(output, slice_path, row_start=1, row_stop=3)
        with slice_path.open(encoding="utf-8", newline="") as handle:
            sliced = list(csv.DictReader(handle, delimiter="\t"))
        assert [int(row["row"]) for row in sliced] == [1, 2]

        sparse_path = root / "sparse.tsv"
        MODULE.read_store(output, sparse_path, row_indices=[0, len(expected) - 1])
        with sparse_path.open(encoding="utf-8", newline="") as handle:
            sparse = list(csv.DictReader(handle, delimiter="\t"))
        assert [int(row["row"]) for row in sparse] == [0, len(expected) - 1]

        # Dense frames and exact recursive escape boundaries.
        boundary_input = root / "boundaries.tsv"
        boundary_positions = [1, 511, 1022, 66_556, 132_091]
        write_rows(boundary_input, [
            ("1", pos, "C", "A", 0.1, 0.05, 0.2, 2.0)
            for pos in boundary_positions
        ])
        boundary_store = root / "boundary-store"
        boundary_manifest = MODULE.write_store(boundary_input, boundary_store)
        boundary_out = root / "boundary-out.tsv"
        MODULE.read_store(boundary_store, boundary_out)
        assert [int(row["base_pair_location"]) for row in read_rows(boundary_out)] == boundary_positions
        assert MODULE.validate_store(boundary_store, full=True)["valid"]

        dense_input = root / "dense.tsv"
        write_rows(dense_input, [
            ("1", pos, "C", "A", 0.1, 0.05, 0.2, 2.0)
            for pos in range(100, 106)
        ])
        dense_store = root / "dense-store"
        MODULE.write_store(dense_input, dense_store)
        dense_out = root / "dense-out.tsv"
        MODULE.read_store(dense_store, dense_out)
        assert len(read_rows(dense_out)) == 6

        single_input = root / "single.tsv"
        write_rows(single_input, [("1", 100, "C", "A", 0.1, 0.05, 0.2, 2.0)])
        single_store = root / "single-store"
        MODULE.write_store(single_input, single_store)
        single_out = root / "single-out.tsv"
        MODULE.read_store(single_store, single_out)
        assert len(read_rows(single_out)) == 1

        frame_input = root / "frame-boundaries.tsv"
        write_rows(frame_input, (
            ("1", pos, "C", "A",
             0.4 if pos == MODULE.VALUE_BLOCK_ROWS + 1 else 0.1,
             0.05, 0.2,
             8.0 if pos == MODULE.VALUE_BLOCK_ROWS + 1 else 2.0)
            for pos in range(1, MODULE.KEY_BLOCK_ROWS + 2)
        ))
        frame_store = root / "frame-boundary-store"
        frame_manifest = MODULE.write_store(frame_input, frame_store)
        assert len(json.loads((frame_store / frame_manifest["files"]["key_index"]).read_text())) == 2
        assert len(json.loads((frame_store / frame_manifest["files"]["value_index"]).read_text())) == 3
        frame_sparse = root / "frame-sparse.tsv"
        MODULE.read_store(
            frame_store, frame_sparse,
            row_indices=[MODULE.VALUE_BLOCK_ROWS - 1, MODULE.VALUE_BLOCK_ROWS,
                         MODULE.KEY_BLOCK_ROWS - 1, MODULE.KEY_BLOCK_ROWS],
        )
        assert [int(row["row"]) for row in read_rows(frame_sparse)] == [
            MODULE.VALUE_BLOCK_ROWS - 1, MODULE.VALUE_BLOCK_ROWS,
            MODULE.KEY_BLOCK_ROWS - 1, MODULE.KEY_BLOCK_ROWS,
        ]
        frame_keys = [
            f"1:{MODULE.VALUE_BLOCK_ROWS}:A:C",
            f"1:{MODULE.VALUE_BLOCK_ROWS + 1}:A:C",
            f"1:{MODULE.KEY_BLOCK_ROWS}:A:C",
            f"1:{MODULE.KEY_BLOCK_ROWS + 1}:A:C",
            "2:200000:A:C",
        ]
        frame_key_sparse = root / "frame-key-sparse.tsv"
        MODULE.read_store(frame_store, frame_key_sparse, identity_keys=frame_keys)
        assert [int(row["row"]) for row in read_rows(frame_key_sparse)] == [
            MODULE.VALUE_BLOCK_ROWS - 1, MODULE.VALUE_BLOCK_ROWS,
            MODULE.KEY_BLOCK_ROWS - 1, MODULE.KEY_BLOCK_ROWS,
        ]
        try:
            MODULE.read_store(frame_store, root / "bad-key.tsv", identity_keys=["rs123"])
        except ValueError as exc:
            assert "chromosome:position:REF:ALT" in str(exc)
        else:
            raise AssertionError("malformed canonical identity key was accepted")
        assert MODULE.validate_store(frame_store, full=True)["valid"]

        # Exception records are owned by one value frame. Repointing a frame
        # at another frame's exception payload must fail before values decode.
        value_index = json.loads(
            (frame_store / frame_manifest["files"]["value_index"]).read_text()
        )
        assert value_index[1]["exception_count"] == 1
        wrong_owner = dict(value_index[0])
        wrong_owner["exceptions"] = value_index[1]["exceptions"]
        wrong_owner["exception_count"] = value_index[1]["exception_count"]
        try:
            MODULE.decode_exception_frame(frame_store, frame_manifest, wrong_owner)
        except ValueError as exc:
            assert "owning value frame" in str(exc)
        else:
            raise AssertionError("misassigned exception frame was accepted")

        # Full validation recomputes key-frame genomic bounds rather than
        # trusting the region index.
        bounds_store = root / "bounds-tamper-store"
        shutil.copytree(frame_store, bounds_store)
        bounds_manifest = json.loads((bounds_store / "manifest.json").read_text())
        bounds_index_path = bounds_store / bounds_manifest["files"]["key_index"]
        bounds_index = json.loads(bounds_index_path.read_text())
        bounds_index[0]["last_position"] = bounds_index[0]["base_position"]
        bounds_index_path.write_text(json.dumps(bounds_index, indent=2) + "\n")
        bounds_blob = bounds_index_path.read_bytes()
        bounds_manifest["integrity"]["files"][bounds_index_path.name] = {
            "bytes": len(bounds_blob), "sha256": MODULE.hashlib.sha256(bounds_blob).hexdigest()
        }
        (bounds_store / "manifest.json").write_text(json.dumps(bounds_manifest, indent=2) + "\n")
        MODULE._seal_manifest(bounds_store)
        try:
            MODULE.validate_store(bounds_store, full=True)
        except ValueError as exc:
            assert "genomic bounds" in str(exc)
        else:
            raise AssertionError("altered key-frame bounds passed full validation")

        # Decode-critical manifest values remain fixed by the format even when
        # an attacker recomputes the detached manifest digest.
        manifest_store = root / "manifest-tamper-store"
        shutil.copytree(dense_store, manifest_store)
        altered = json.loads((manifest_store / "manifest.json").read_text())
        altered["semantic_codec"]["z_range"] = [-35, 35]
        (manifest_store / "manifest.json").write_text(json.dumps(altered, indent=2) + "\n")
        MODULE._seal_manifest(manifest_store)
        try:
            MODULE.load_manifest(manifest_store)
        except ValueError as exc:
            assert "semantic codec constants" in str(exc)
        else:
            raise AssertionError("altered semantic constants were trusted")

        # A whole-file Z projection must not open identity, EAF, or SE payloads.
        dense_manifest = json.loads((dense_store / "manifest.json").read_text())
        hidden = []
        for name in ("key_gap", "eaf", "se"):
            original = dense_store / dense_manifest["files"][name]
            moved = original.with_suffix(original.suffix + ".hidden")
            original.rename(moved)
            hidden.append((original, moved))
        projected = root / "projected-z.tsv"
        MODULE.read_store(dense_store, projected, columns=["z"])
        assert list(read_rows(projected)[0]) == ["row", "z"]
        for original, moved in hidden:
            moved.rename(original)

        # Missing values are independent: EAF missing cannot erase finite SE,
        # and missing SE cannot be fabricated from the residual predictor.
        missing_input = root / "missing.tsv"
        write_rows(missing_input, [
            ("1", 200, "C", "A", ".", ".", 0.2, 2.0),
            ("1", 201, "G", "A", 0.2, 0.1, ".", 2.0),
            ("1", 202, "T", "A", ".", 0.2, 0.2, "."),
            ("1", 203, "A", "C", 0.2, 0.1, 0.2, 2.0),
        ])
        missing_store = root / "missing-store"
        MODULE.write_store(missing_input, missing_store)
        missing_out = root / "missing-out.tsv"
        MODULE.read_store(missing_store, missing_out)
        missing = read_rows(missing_out)
        assert math.isnan(float(missing[0]["standard_error"]))
        assert abs(float(missing[1]["standard_error"]) - 0.1) < 1e-6
        assert math.isnan(float(missing[1]["effect_allele_frequency"]))
        assert abs(float(missing[2]["standard_error"]) - 0.2) < 0.005
        assert math.isnan(float(missing[2]["z"])) and math.isnan(float(missing[2]["beta"]))

        invalid_input = root / "invalid-se.tsv"
        write_rows(invalid_input, [("1", 300, "C", "A", 0.1, 0.0, 0.2, 2.0)])
        try:
            MODULE.write_store(invalid_input, root / "invalid-store")
        except ValueError as exc:
            assert "standard_error" in str(exc)
        else:
            raise AssertionError("zero standard error was accepted")

        # Column projection avoids unrelated numeric streams.
        identity_out = root / "identity-only.tsv"
        MODULE.read_store(dense_store, identity_out, columns=["chromosome", "base_pair_location"])
        assert list(read_rows(identity_out)[0]) == ["row", "chromosome", "base_pair_location"]

        # Any byte corruption is caught by file checksums.
        z_path = dense_store / dense_manifest_file(dense_store, "z")
        blob = bytearray(z_path.read_bytes())
        blob[len(blob) // 2] ^= 0x01
        z_path.write_bytes(blob)
        try:
            MODULE.validate_store(dense_store)
        except ValueError as exc:
            assert "checksum" in str(exc)
        else:
            raise AssertionError("corrupted middle payload passed validation")

        print(json.dumps({"valid": True, "rows": len(decoded), "frames": len(key_index),
                          "adversarial_cases": 15}))


def dense_manifest_file(store: Path, name: str) -> str:
    return json.loads((store / "manifest.json").read_text())["files"][name]


if __name__ == "__main__":
    main()
