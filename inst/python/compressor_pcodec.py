#!/usr/bin/env python3
"""Block-indexed Pcodec backend for CompreSSoR.

The file format deliberately stores no per-row textual variant ID.  Each row
has a lossless identity key made from a GRCh38 primary-chromosome position and
the complete REF->ALT substitution code.  The chromosome lengths used to
decode the global position are embedded in the manifest, so decoding does not
require a FASTA, canonical spine, or shared study-side reference.

This module is intentionally dependency-light at the boundary: the optional
backend requires numpy, pcodec, and zstandard in the configured Python
environment.  R calls the write/read subcommands; the implementation is also
usable directly for smoke tests and future Python readers.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import struct
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable


CHROM_LENGTHS = {
    "1": 248_956_422,
    "2": 242_193_529,
    "3": 198_295_559,
    "4": 190_214_555,
    "5": 181_538_259,
    "6": 170_805_979,
    "7": 159_345_973,
    "8": 145_138_636,
    "9": 138_394_717,
    "10": 133_797_422,
    "11": 135_086_622,
    "12": 133_275_309,
    "13": 114_364_328,
    "14": 107_043_718,
    "15": 101_991_189,
    "16": 90_338_345,
    "17": 83_257_441,
    "18": 80_373_285,
    "19": 58_617_616,
    "20": 64_444_167,
    "21": 46_709_983,
    "22": 50_818_468,
    "X": 156_040_895,
    "Y": 57_227_415,
}
BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
KEY_BLOCK_ROWS = 131_072
VALUE_BLOCK_ROWS = 65_536
GAP_ESCAPE = 511
GAP_ESCAPE_THRESHOLD = 510
SECOND_GAP_ESCAPE = 65_535
SECOND_GAP_ESCAPE_THRESHOLD = 65_534
FORMAT = "CompreSSoR"
VERSION = "0.2.0-pcodec"
MANIFEST_CHECKSUM = "manifest.sha256"
EXPECTED_FILES = {
    "key_index": "key.index.json",
    "key_gap": "key.gap.frames",
    "key_gap_exceptions": "key.gap.exceptions.frames",
    "key_gap_exceptions2": "key.gap.exceptions2.frames",
    "key_substitution": "key.substitution.frames",
    "value_index": "value.index.json",
    "z": "z.frames",
    "eaf": "eaf.frames",
    "se": "se.frames",
    "exceptions": "exceptions.frames",
}


def _dependencies():
    try:
        import numpy as np
        import zstandard as zstd
        from pcodec import ChunkConfig, standalone
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "The pcodec backend requires numpy, pcodec, and zstandard in "
            "the Python environment selected by COMPRESSOR_PYTHON"
        ) from exc
    return np, zstd, ChunkConfig, standalone


def chromosome_offsets():
    offsets = {}
    current = 0
    for chromosome, length in CHROM_LENGTHS.items():
        offsets[chromosome] = current
        current += length
    return offsets


def _seal_manifest(store: Path) -> None:
    """Write the detached digest checked before trusting manifest metadata."""
    blob = (store / "manifest.json").read_bytes()
    (store / MANIFEST_CHECKSUM).write_text(hashlib.sha256(blob).hexdigest() + "\n")


def _validate_manifest_contract(manifest: dict[str, Any]) -> None:
    """Reject altered decode-critical constants for the fixed v0.2 format."""
    if manifest.get("format") != FORMAT:
        raise ValueError("not a CompreSSoR store")
    if manifest.get("format_version") != VERSION:
        raise ValueError(
            f"unsupported Pcodec format version: {manifest.get('format_version')!r}; expected {VERSION!r}"
        )
    if manifest.get("backend") != "pcodec" or manifest.get("profile") != "standard":
        raise ValueError("store is not a standard Pcodec store")
    try:
        rows = int(manifest["rows"])
        n_rows = int(manifest["n_rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("manifest row count is missing or invalid") from exc
    if rows < 1 or n_rows != rows:
        raise ValueError("manifest row counts disagree or are empty")
    if int(manifest.get("key_block_rows", -1)) != KEY_BLOCK_ROWS:
        raise ValueError("manifest key frame size is incompatible with this format")
    if int(manifest.get("value_block_rows", -1)) != VALUE_BLOCK_ROWS:
        raise ValueError("manifest value frame size is incompatible with this format")

    identity = manifest.get("identity", {})
    expected_identity = {
        "encoding": "global_grch38_primary_position_plus_full_ref_alt_code",
        "external_reference_required": False,
        "chromosome_lengths": CHROM_LENGTHS,
        "chromosome_offsets": chromosome_offsets(),
        "effect_allele_is_alt": True,
        "other_allele_is_ref": True,
    }
    if identity != expected_identity:
        raise ValueError("manifest identity constants differ from the fixed GRCh38 contract")
    chromosome_order = list(CHROM_LENGTHS)
    if (list(identity.get("chromosome_lengths", {})) != chromosome_order or
            list(identity.get("chromosome_offsets", {})) != chromosome_order):
        raise ValueError("manifest chromosomes are not in canonical GRCh38 order")

    semantic = manifest.get("semantic_codec", {})
    fixed_semantic = {
        "name": "z9/eaf8/se6",
        "z_bits": 9,
        "eaf_bits": 8,
        "se_bits": 6,
        "z_range": [-3.5, 3.5],
        "exception_precision": "float32",
        "beta": "derived as z * standard_error",
        "p_value": "derived as erfc(abs(z) / sqrt(2))",
    }
    if any(semantic.get(name) != value for name, value in fixed_semantic.items()):
        raise ValueError("manifest semantic codec constants differ from the fixed Z9/EAF8/SE6 contract")
    centers = semantic.get("block_centers_log2_residual")
    expected_centers = (rows + VALUE_BLOCK_ROWS - 1) // VALUE_BLOCK_ROWS
    if not isinstance(centers, list) or len(centers) != expected_centers:
        raise ValueError("manifest has the wrong number of SE block centres")
    try:
        finite_centers = all(math.isfinite(float(value)) for value in centers)
        exception_rows = int(semantic["exception_rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("manifest semantic metadata is invalid") from exc
    if not finite_centers or exception_rows < 0 or exception_rows > rows:
        raise ValueError("manifest semantic metadata is outside the format contract")
    if manifest.get("files") != EXPECTED_FILES:
        raise ValueError("manifest file layout differs from the fixed Pcodec layout")
    integrity = manifest.get("integrity", {})
    records = integrity.get("files", {})
    if integrity.get("algorithm") != "sha256" or set(records) != set(EXPECTED_FILES.values()):
        raise ValueError("manifest integrity table is incomplete")
    for relative in EXPECTED_FILES.values():
        record = records.get(relative, {})
        digest = record.get("sha256")
        try:
            size = int(record.get("bytes", -1))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"manifest integrity size is invalid for {relative}") from exc
        if size < 0 or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"manifest integrity record is invalid for {relative}")


def _as_float(value: str) -> float:
    if value is None or value == "" or value == "." or value.strip().lower() in {"na", "nan", "null"}:
        return float("nan")
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(f"cannot parse numeric summary-statistic value: {value!r}") from exc


def read_tsv(path: Path):
    np, _, _, _ = _dependencies()
    opener = gzip.open if path.name.lower().endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("input TSV is empty") from exc
        names = {name.strip().lower(): i for i, name in enumerate(header)}

        def column(*aliases):
            for alias in aliases:
                if alias.lower() in names:
                    return names[alias.lower()]
            return None

        chrom_i = column("chromosome", "chrom", "chr")
        pos_i = column("base_pair_location", "position", "pos", "bp")
        effect_i = column("effect_allele", "ea", "a1", "alt")
        other_i = column("other_allele", "oa", "a2", "ref")
        beta_i = column("beta", "b", "effect")
        se_i = column("standard_error", "se", "stderr", "sebeta")
        eaf_i = column("effect_allele_frequency", "eaf", "af")
        z_i = column("z", "zscore", "z_score")
        required = {
            "chromosome": chrom_i,
            "base_pair_location": pos_i,
            "effect_allele": effect_i,
            "other_allele": other_i,
            "standard_error": se_i,
            "effect_allele_frequency": eaf_i,
        }
        missing = [key for key, index in required.items() if index is None]
        if missing:
            raise ValueError("input TSV is missing: " + ", ".join(missing))

        chromosomes = []
        positions = []
        effects = []
        others = []
        betas = []
        ses = []
        eafs = []
        zs = []
        for line_number, row in enumerate(reader, start=2):
            if len(row) < len(header):
                raise ValueError(f"row {line_number} has fewer fields than the header")
            chrom = row[chrom_i].strip().removeprefix("chr").upper()
            if chrom not in CHROM_LENGTHS:
                raise ValueError(f"unsupported GRCh38 primary chromosome at row {line_number}: {chrom}")
            ref = row[other_i].strip().upper()
            alt = row[effect_i].strip().upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in BASE_CODE or alt not in BASE_CODE:
                raise ValueError(f"only biallelic SNVs are currently supported at row {line_number}")
            position = int(row[pos_i])
            if position < 1 or position > CHROM_LENGTHS[chrom]:
                raise ValueError(f"position outside GRCh38 chromosome {chrom} at row {line_number}")
            beta = _as_float(row[beta_i]) if beta_i is not None else float("nan")
            se = _as_float(row[se_i])
            eaf = _as_float(row[eaf_i])
            z = _as_float(row[z_i]) if z_i is not None else float("nan")
            if not math.isnan(beta) and not math.isfinite(beta):
                raise ValueError(f"beta must be finite when supplied at row {line_number}")
            if not math.isnan(z) and not math.isfinite(z):
                raise ValueError(f"z must be finite when supplied at row {line_number}")
            if not math.isnan(se) and (not math.isfinite(se) or se <= 0):
                raise ValueError(f"standard_error must be finite and positive at row {line_number}")
            if not math.isnan(eaf) and (not math.isfinite(eaf) or eaf < 0 or eaf > 1):
                raise ValueError(f"effect_allele_frequency must be between 0 and 1 at row {line_number}")
            if not math.isfinite(z) and math.isfinite(beta) and math.isfinite(se) and se > 0:
                z = beta / se
            if not math.isfinite(beta) and math.isfinite(z) and math.isfinite(se) and se > 0:
                beta = z * se
            chromosomes.append(chrom)
            positions.append(position)
            effects.append(alt)
            others.append(ref)
            betas.append(beta)
            ses.append(se)
            eafs.append(eaf)
            zs.append(z)
    return {
        "chromosome": np.asarray(chromosomes, dtype=object),
        "position": np.asarray(positions, dtype=np.uint32),
        "effect": np.asarray(effects, dtype=object),
        "other": np.asarray(others, dtype=object),
        "beta": np.asarray(betas, dtype=np.float64),
        "se": np.asarray(ses, dtype=np.float64),
        "eaf": np.asarray(eafs, dtype=np.float64),
        "z": np.asarray(zs, dtype=np.float64),
    }


def identity_keys(data):
    np, _, _, _ = _dependencies()
    offsets = chromosome_offsets()
    global_position = np.asarray(
        [offsets[chrom] + int(position) - 1 for chrom, position in zip(data["chromosome"], data["position"])],
        dtype=np.uint64,
    )
    refs = np.asarray([BASE_CODE[base] for base in data["other"]], dtype=np.uint64)
    alts = np.asarray([BASE_CODE[base] for base in data["effect"]], dtype=np.uint64)
    if np.any(refs == alts):
        raise ValueError("REF and ALT must differ")
    keys = (global_position << np.uint64(4)) | (refs << np.uint64(2)) | alts
    order = np.argsort(keys, kind="stable")
    if len(np.unique(keys)) != len(keys):
        raise ValueError("duplicate full REF/ALT identity keys")
    return keys[order], order


def pcodec_compress(values):
    np, _, ChunkConfig, standalone = _dependencies()
    config = ChunkConfig(enable_8_bit=True)
    return standalone.simple_compress(np.ascontiguousarray(values), config)


def pcodec_decompress(blob):
    _, _, _, standalone = _dependencies()
    return standalone.simple_decompress(blob)


def append_blob(handle, blob: bytes) -> dict[str, int]:
    offset = handle.tell()
    handle.write(blob)
    return {"offset": offset, "length": len(blob), "sha256": hashlib.sha256(blob).hexdigest()}


def encode_key_frames(keys, output: Path, block_rows: int = KEY_BLOCK_ROWS):
    np, _, _, _ = _dependencies()
    offsets = chromosome_offsets()
    paths = {
        "gap": output / "key.gap.frames",
        "gap_exceptions": output / "key.gap.exceptions.frames",
        "gap_exceptions2": output / "key.gap.exceptions2.frames",
        "substitution": output / "key.substitution.frames",
    }
    handles = {name: path.open("wb") for name, path in paths.items()}
    index = []
    try:
        for start in range(0, len(keys), block_rows):
            stop = min(start + block_rows, len(keys))
            frame_keys = keys[start:stop]
            positions = (frame_keys >> np.uint64(4)).astype(np.uint64)
            gaps = np.empty(len(positions), dtype=np.uint32)
            gaps[0] = 0
            if len(gaps) > 1:
                gaps[1:] = positions[1:] - positions[:-1]
            mask1 = gaps > GAP_ESCAPE_THRESHOLD
            core1 = np.where(mask1, GAP_ESCAPE, gaps).astype(np.uint16)
            exceptions1 = gaps[mask1].astype(np.uint32)
            mask2 = exceptions1 > SECOND_GAP_ESCAPE_THRESHOLD
            core2 = np.where(mask2, SECOND_GAP_ESCAPE, exceptions1).astype(np.uint16)
            exceptions2 = exceptions1[mask2].astype(np.uint32)
            sub = (frame_keys & np.uint64(15)).astype(np.uint8)
            blobs = {
                "gap": pcodec_compress(core1),
                "gap_exceptions": pcodec_compress(core2) if len(core2) else b"",
                "gap_exceptions2": pcodec_compress(exceptions2) if len(exceptions2) else b"",
                "substitution": pcodec_compress(sub),
            }
            locations = {name: append_blob(handles[name], blob) for name, blob in blobs.items()}
            first_position = int(positions[0])
            last_position = int(positions[-1])
            index.append({
                "row_start": start,
                "row_stop": stop,
                "base_position": first_position,
                "last_position": last_position,
                "gap": locations["gap"],
                "gap_exceptions": locations["gap_exceptions"],
                "gap_exceptions2": locations["gap_exceptions2"],
                "substitution": locations["substitution"],
                "exception_count": int(np.count_nonzero(mask1)),
                "exception2_count": int(np.count_nonzero(mask2)),
            })
    finally:
        for handle in handles.values():
            handle.close()
    return index, paths


def quantise_values(data, block_rows: int = VALUE_BLOCK_ROWS):
    np, _, _, _ = _dependencies()
    n = len(data["beta"])
    z = data["z"].copy()
    se = data["se"].copy()
    eaf = data["eaf"].copy()
    valid_eaf = np.isfinite(eaf) & (eaf >= 0.0) & (eaf <= 1.0)
    safe_eaf = np.clip(np.where(valid_eaf, eaf, 0.5), 1e-12, 1.0 - 1e-12)
    eaf_max_code = 255
    eaf_codes = np.rint(eaf_max_code * (2.0 / np.pi) * np.arcsin(np.sqrt(safe_eaf))).astype(np.uint8)
    eaf_decoded = np.sin(np.pi * eaf_codes.astype(np.float64) / (2.0 * eaf_max_code)) ** 2
    z_valid = np.isfinite(z)
    z_min, z_max = -3.5, 3.5
    z_count = 510
    z_exception = 511
    z_missing = 510
    z_step = (z_max - z_min) / z_count
    z_central = z_valid & (z >= z_min) & (z < z_max)
    z_codes = np.full(n, z_missing, dtype=np.uint16)
    z_codes[z_central] = np.floor((z[z_central] - z_min) / z_step).astype(np.uint16)
    z_codes[~z_central & z_valid] = z_exception
    safe_eaf_for_se = np.clip(eaf_decoded, 1e-12, 1.0 - 1e-12)
    valid_se = np.isfinite(se) & (se > 0)
    residual = np.full(n, np.nan, dtype=np.float64)
    residual_ready = valid_se & valid_eaf
    residual[residual_ready] = (
        np.log2(se[residual_ready])
        + 0.5 * np.log2(
            2.0 * safe_eaf_for_se[residual_ready] * (1.0 - safe_eaf_for_se[residual_ready])
        )
    )
    se_count = 62
    se_exception = 63
    se_missing = 62
    se_step = 2.0 / se_count
    se_codes = np.full(n, se_missing, dtype=np.uint8)
    centers = []
    se_central = np.isfinite(residual) & valid_se & valid_eaf
    for start in range(0, n, block_rows):
        stop = min(start + block_rows, n)
        block_ok = se_central[start:stop]
        center = float(np.median(residual[start:stop][block_ok])) if np.any(block_ok) else 0.0
        centers.append(center)
        delta = residual[start:stop] - center
        in_range = block_ok & (delta >= -1.0) & (delta < 1.0)
        local = np.full(stop - start, se_missing, dtype=np.uint8)
        local[in_range] = np.floor((delta[in_range] + 1.0) / se_step).astype(np.uint8)
        local[block_ok & ~in_range] = se_exception
        local[valid_se[start:stop] & ~valid_eaf[start:stop]] = se_exception
        se_codes[start:stop] = local
    z_exception_mask = (~z_central & z_valid)
    se_exception_mask = (se_codes == se_exception)
    exception_mask = z_exception_mask | se_exception_mask | (~valid_eaf)
    rows = np.flatnonzero(exception_mask)
    exception_dtype = np.dtype(
        [("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"), ("eaf", "<f4"), ("flags", "u1")]
    )
    exceptions = np.empty(len(rows), dtype=exception_dtype)
    exceptions["row"] = rows.astype(np.uint32)
    exceptions["z"] = z[rows].astype(np.float32)
    exceptions["log2se"] = np.log2(np.where(np.isfinite(se[rows]) & (se[rows] > 0), se[rows], 1.0)).astype(np.float32)
    exceptions["eaf"] = eaf[rows].astype(np.float32)
    flags = np.zeros(len(rows), dtype=np.uint8)
    flags |= z_exception_mask[rows].astype(np.uint8)
    flags |= (se_exception_mask[rows].astype(np.uint8) << 1)
    flags |= ((~valid_eaf[rows]).astype(np.uint8) << 2)
    exceptions["flags"] = flags
    if np.any(z_exception_mask[rows] & ~np.isfinite(exceptions["z"])):
        raise ValueError("a Z exception cannot be represented as finite float32")
    if np.any(se_exception_mask[rows] & ~np.isfinite(exceptions["log2se"])):
        raise ValueError("an SE exception cannot be represented as finite float32 log2(SE)")
    finite_eaf_exception = (~valid_eaf[rows]) & np.isfinite(eaf[rows])
    if np.any(finite_eaf_exception & ~np.isfinite(exceptions["eaf"])):
        raise ValueError("an EAF exception cannot be represented as finite float32")
    return {"z": z_codes, "eaf": eaf_codes, "se": se_codes, "centers": centers, "exceptions": exceptions}


def encode_value_frames(values, output: Path, block_rows: int = VALUE_BLOCK_ROWS):
    _, zstd, _, _ = _dependencies()
    paths = {name: output / f"{name}.frames" for name in ("z", "eaf", "se", "exceptions")}
    handles = {name: path.open("wb") for name, path in paths.items()}
    index = []
    exception_rows = values["exceptions"]["row"]
    try:
        for start in range(0, len(values["z"]), block_rows):
            stop = min(start + block_rows, len(values["z"]))
            locations = {}
            for name in ("z", "eaf", "se"):
                locations[name] = append_blob(handles[name], pcodec_compress(values[name][start:stop]))
            lo = int(exception_rows.searchsorted(start))
            hi = int(exception_rows.searchsorted(stop))
            exception_blob = zstd.ZstdCompressor(level=19).compress(values["exceptions"][lo:hi].tobytes())
            locations["exceptions"] = append_blob(handles["exceptions"], exception_blob)
            locations["exception_count"] = hi - lo
            index.append({"row_start": start, "row_stop": stop, **locations})
    finally:
        for handle in handles.values():
            handle.close()
    return index, paths


def write_store(input_path: Path, output: Path, metadata: dict[str, Any] | None = None):
    np, _, _, _ = _dependencies()
    output.mkdir(parents=True, exist_ok=True)
    data = read_tsv(input_path)
    keys, order = identity_keys(data)
    ordered = {name: values[order] for name, values in data.items()}
    values = quantise_values(ordered)
    key_index, key_paths = encode_key_frames(keys, output)
    value_index, value_paths = encode_value_frames(values, output)
    manifest = {
        "format": FORMAT,
        "format_version": VERSION,
        "backend": "pcodec",
        "profile": "standard",
        "rows": len(keys),
        "n_rows": len(keys),
        "key_block_rows": KEY_BLOCK_ROWS,
        "value_block_rows": VALUE_BLOCK_ROWS,
        "identity": {
            "encoding": "global_grch38_primary_position_plus_full_ref_alt_code",
            "external_reference_required": False,
            "chromosome_lengths": CHROM_LENGTHS,
            "chromosome_offsets": chromosome_offsets(),
            "effect_allele_is_alt": True,
            "other_allele_is_ref": True,
        },
        "semantic_codec": {
            "name": "z9/eaf8/se6",
            "z_bits": 9,
            "eaf_bits": 8,
            "se_bits": 6,
            "z_range": [-3.5, 3.5],
            "block_centers_log2_residual": values["centers"],
            "exception_rows": int(len(values["exceptions"])),
            "exception_precision": "float32",
            "beta": "derived as z * standard_error",
            "p_value": "derived as erfc(abs(z) / sqrt(2))",
        },
        "files": EXPECTED_FILES,
        "input_columns": ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "z"],
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pcodec": importlib.metadata.version("pcodec"),
            "zstandard": importlib.metadata.version("zstandard"),
        },
        "metadata": metadata or {},
    }
    (output / "key.index.json").write_text(json.dumps(key_index, indent=2) + "\n")
    (output / "value.index.json").write_text(json.dumps(value_index, indent=2) + "\n")
    integrity = {}
    for relative in manifest["files"].values():
        path = output / relative
        blob = path.read_bytes()
        integrity[relative] = {"bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest()}
    manifest["integrity"] = {"algorithm": "sha256", "files": integrity}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    _seal_manifest(output)
    return manifest


def _verify_file(store: Path, manifest: dict[str, Any], relative: str):
    record = manifest.get("integrity", {}).get("files", {}).get(relative)
    if not record:
        raise ValueError(f"manifest has no integrity record for {relative}")
    path = store / relative
    if not path.is_file():
        raise ValueError(f"store file is missing: {relative}")
    expected_bytes = int(record["bytes"])
    if path.stat().st_size != expected_bytes:
        raise ValueError(f"file size mismatch: {relative}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != record["sha256"]:
        raise ValueError(f"file checksum mismatch: {relative}")


def load_manifest(store: Path, verify_indexes: bool = True):
    path = store / "manifest.json"
    if not path.is_file():
        raise ValueError("store is missing manifest.json")
    checksum_path = store / MANIFEST_CHECKSUM
    if not checksum_path.is_file():
        raise ValueError(f"store is missing {MANIFEST_CHECKSUM}")
    expected = checksum_path.read_text().strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError("manifest checksum record is malformed")
    blob = path.read_bytes()
    if hashlib.sha256(blob).hexdigest() != expected:
        raise ValueError("manifest checksum mismatch")
    manifest = json.loads(blob)
    _validate_manifest_contract(manifest)
    if verify_indexes:
        _verify_file(store, manifest, manifest["files"]["key_index"])
        _verify_file(store, manifest, manifest["files"]["value_index"])
    return manifest


def _read_blob(path: Path, location: dict[str, int]) -> bytes:
    offset = int(location["offset"])
    length = int(location["length"])
    if offset < 0 or length < 0 or offset + length > path.stat().st_size:
        raise ValueError(f"frame byte range is outside {path.name}")
    with path.open("rb") as handle:
        handle.seek(offset)
        blob = handle.read(length)
    if len(blob) != length:
        raise ValueError(f"short frame read from {path.name}")
    expected = location.get("sha256")
    if expected and hashlib.sha256(blob).hexdigest() != expected:
        raise ValueError(f"frame checksum mismatch in {path.name}")
    return blob


def decode_key_frame(store: Path, manifest: dict[str, Any], frame: dict[str, Any]):
    np, _, _, _ = _dependencies()
    files = manifest["files"]
    core = pcodec_decompress(_read_blob(store / files["key_gap"], frame["gap"])).astype(np.uint32)
    exceptions1 = (
        pcodec_decompress(_read_blob(store / files["key_gap_exceptions"], frame["gap_exceptions"])).astype(np.uint32)
        if frame["gap_exceptions"]["length"] else np.empty(0, dtype=np.uint32)
    )
    exceptions2 = pcodec_decompress(_read_blob(store / files["key_gap_exceptions2"], frame["gap_exceptions2"])).astype(np.uint32) if frame["gap_exceptions2"]["length"] else np.empty(0, dtype=np.uint32)
    first = np.flatnonzero(core == GAP_ESCAPE)
    if len(first) != len(exceptions1):
        raise ValueError("first-level key exception count does not match its sentinels")
    core[first] = exceptions1
    second = np.flatnonzero(core == SECOND_GAP_ESCAPE)
    if len(second) != len(exceptions2):
        raise ValueError("second-level key exception count does not match its sentinels")
    core[second] = exceptions2
    positions = np.empty(frame["row_stop"] - frame["row_start"], dtype=np.uint64)
    positions[0] = frame["base_position"]
    if len(positions) > 1:
        positions[1:] = frame["base_position"] + np.cumsum(core[1:].astype(np.uint64), dtype=np.uint64)
    substitutions = pcodec_decompress(_read_blob(store / files["key_substitution"], frame["substitution"])).astype(np.uint8)
    return (positions << np.uint64(4)) | substitutions


def decode_value_frame(
    store: Path,
    manifest: dict[str, Any],
    frame: dict[str, Any],
    streams=("z", "eaf", "se"),
):
    files = manifest["files"]
    return {
        name: pcodec_decompress(_read_blob(store / files[name], frame[name]))
        for name in streams
    }


def _keys_to_columns(keys, manifest):
    np, _, _, _ = _dependencies()
    offsets = np.asarray(list(chromosome_offsets().values()), dtype=np.uint64)
    positions = keys >> np.uint64(4)
    chrom_ids = np.searchsorted(offsets, positions, side="right") - 1
    local_positions = (positions - offsets[chrom_ids] + 1).astype(np.int32)
    substitutions = keys & np.uint64(15)
    refs = (substitutions >> np.uint64(2)).astype(np.uint8)
    alts = (substitutions & np.uint64(3)).astype(np.uint8)
    return (chrom_ids + 1).astype(np.uint8), local_positions, alts, refs


def _keys_to_rows(keys, manifest):
    np, _, _, _ = _dependencies()
    chromosome_codes, local_positions, alts, refs = _keys_to_columns(keys, manifest)
    chromosome_names = list(CHROM_LENGTHS)
    bases = np.asarray(["A", "C", "G", "T"])
    return ([chromosome_names[index - 1] for index in chromosome_codes],
            local_positions, bases[alts], bases[refs])


def _write_binary_bridge(
    output: Path,
    output_rows,
    requested_columns,
    output_values,
    semantic_codes=None,
    exceptions=None,
    manifest=None,
):
    """Write a temporary columnar bridge that R can load with readBin()."""
    np, _, _, _ = _dependencies()
    output.mkdir(parents=True, exist_ok=False)
    files = {}

    def write_column(name, values, dtype, dtype_name):
        filename = f"{name}.bin"
        array = np.ascontiguousarray(values, dtype=dtype)
        array.tofile(output / filename)
        files[name] = {"file": filename, "dtype": dtype_name, "length": int(len(array))}

    identity_dtypes = {
        "chromosome": ("u1", "uint8"),
        "base_pair_location": ("<i4", "int32"),
        "effect_allele": ("u1", "uint8"),
        "other_allele": ("u1", "uint8"),
    }
    for name, (dtype, dtype_name) in identity_dtypes.items():
        if name in requested_columns:
            write_column(name, output_values[name], dtype, dtype_name)

    codec = None
    if semantic_codes is not None:
        code_specs = {
            "z": ("z_code", "<u2", "uint16"),
            "se": ("se_code", "u1", "uint8"),
            "eaf": ("eaf_code", "u1", "uint8"),
        }
        for stream, (name, dtype, dtype_name) in code_specs.items():
            if stream in semantic_codes:
                write_column(name, semantic_codes[stream], dtype, dtype_name)
        write_column("exception_row", exceptions["row"], "<i4", "int32")
        write_column("exception_z", exceptions["z"], "<f4", "float32")
        write_column("exception_se", np.exp2(exceptions["log2se"].astype(np.float64)),
                     "<f8", "float64")
        write_column("exception_eaf", exceptions["eaf"], "<f4", "float32")
        write_column("exception_flags", exceptions["flags"], "u1", "uint8")
        semantic = manifest["semantic_codec"]
        codec = {
            "encoding": "semantic_codes",
            "z_bits": int(semantic["z_bits"]),
            "se_bits": int(semantic["se_bits"]),
            "eaf_bits": int(semantic["eaf_bits"]),
            "z_min": float(semantic["z_range"][0]),
            "z_max": float(semantic["z_range"][1]),
            "z_count": 2**int(semantic["z_bits"]) - 2,
            "se_count": 2**int(semantic["se_bits"]) - 2,
            "eaf_count": 2**int(semantic["eaf_bits"]) - 1,
            "block_rows": int(manifest["value_block_rows"]),
            "block_centers_log2_residual": semantic["block_centers_log2_residual"],
            "streams": sorted(semantic_codes),
        }
    else:
        physical_numeric = []
        if any(name in requested_columns for name in ("z", "beta", "p_value")):
            physical_numeric.append("z")
        if any(name in requested_columns for name in ("standard_error", "beta")):
            physical_numeric.append("standard_error")
        if "effect_allele_frequency" in requested_columns:
            physical_numeric.append("effect_allele_frequency")
        if "beta" in requested_columns:
            physical_numeric.append("beta")
        if "p_value" in requested_columns:
            physical_numeric.append("p_value")
        for name in physical_numeric:
            write_column(name, output_values[name], "<f8", "float64")

    (output / "bridge.json").write_text(json.dumps({
        "format": "CompreSSoR-binary-bridge",
        "version": 1,
        "rows": int(len(output_rows)),
        "requested_columns": requested_columns,
        "files": files,
        "codec": codec,
    }, indent=2))


def decode_full_semantic_codes(store: Path, manifest: dict[str, Any], streams):
    np, _, _, _ = _dependencies()
    streams = set(streams)
    if "se" in streams:
        streams.add("eaf")
    value_index = json.loads((store / manifest["files"]["value_index"]).read_text())
    parts = {name: [] for name in streams}
    exception_parts = []
    for frame in value_index:
        decoded = decode_value_frame(store, manifest, frame, streams=streams)
        for name in parts:
            parts[name].append(decoded[name])
        exception_parts.append(decode_exception_frame(store, manifest, frame))
    codes = {}
    for name in streams:
        dtype = np.uint16 if name == "z" else np.uint8
        codes[name] = np.concatenate(parts[name]) if parts[name] else np.empty(0, dtype=dtype)
    exception_dtype = np.dtype([
        ("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"),
        ("eaf", "<f4"), ("flags", "u1")
    ])
    exceptions = (np.concatenate(exception_parts) if exception_parts
                  else np.empty(0, dtype=exception_dtype))
    return codes, exceptions


def decode_exception_frame(store: Path, manifest: dict[str, Any], frame: dict[str, Any]):
    np, zstd, _, _ = _dependencies()
    relative = manifest["files"]["exceptions"]
    raw = zstd.ZstdDecompressor().decompress(
        _read_blob(store / relative, frame["exceptions"])
    )
    dtype = np.dtype([("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"), ("eaf", "<f4"), ("flags", "u1")])
    if len(raw) % dtype.itemsize:
        raise ValueError("exception frame length is not a whole number of records")
    decoded = np.frombuffer(raw, dtype=dtype)
    if len(decoded) != int(frame["exception_count"]):
        raise ValueError("decoded exception frame count differs from its index")
    if len(decoded):
        rows = decoded["row"].astype(np.int64)
        start = int(frame["row_start"])
        stop = int(frame["row_stop"])
        if np.any(rows < start) or np.any(rows >= stop):
            raise ValueError("exception row is outside its owning value frame")
        if np.any(rows[1:] <= rows[:-1]):
            raise ValueError("exception rows are not strictly increasing within their frame")
        if np.any(decoded["flags"] == 0) or np.any((decoded["flags"] & ~7) != 0):
            raise ValueError("exception flags are invalid")
    return decoded


def decode_values(
    store: Path,
    manifest: dict[str, Any],
    row_start: int,
    row_stop: int,
    row_indices=None,
    needed=("z", "se", "eaf"),
):
    np, _, _, _ = _dependencies()
    needed = set(needed)
    if not needed:
        return {}
    streams = set(needed)
    if "se" in streams:
        streams.add("eaf")
    value_index = json.loads((store / manifest["files"]["value_index"]).read_text())
    if row_indices is None:
        target_rows = np.arange(row_start, row_stop, dtype=np.int64)
    else:
        target_rows = np.asarray(sorted(set(int(row) for row in row_indices)), dtype=np.int64)
    selected = [frame for frame in value_index if np.any(
        (target_rows >= frame["row_start"]) & (target_rows < frame["row_stop"])
    )]
    parts = {name: [] for name in streams}
    exception_parts = []
    for frame in selected:
        decoded = decode_value_frame(store, manifest, frame, streams=streams)
        frame_rows = target_rows[(target_rows >= frame["row_start"]) & (target_rows < frame["row_stop"])]
        local = frame_rows - frame["row_start"]
        for name in streams:
            parts[name].append(decoded[name][local])
        exception_parts.append(decode_exception_frame(store, manifest, frame))
    codes = {}
    for name in streams:
        dtype = np.uint16 if name == "z" else np.uint8
        codes[name] = np.concatenate(parts[name]) if parts[name] else np.empty(0, dtype=dtype)
    result = {}
    if "z" in streams:
        z_count = 2**manifest["semantic_codec"]["z_bits"] - 2
        z = np.full(len(codes["z"]), np.nan, dtype=np.float64)
        z_ok = codes["z"] < z_count
        z_min, z_max = manifest["semantic_codec"]["z_range"]
        z_step = (float(z_max) - float(z_min)) / z_count
        z[z_ok] = float(z_min) + (codes["z"][z_ok].astype(np.float64) + 0.5) * z_step
        result["z"] = z
    if "eaf" in streams:
        eaf_max = float(2**manifest["semantic_codec"]["eaf_bits"] - 1)
        result["eaf"] = np.sin(np.pi * codes["eaf"].astype(np.float64) / (2.0 * eaf_max)) ** 2
    if "se" in streams:
        se_count = 2**manifest["semantic_codec"]["se_bits"] - 2
        se = np.full(len(codes["se"]), np.nan, dtype=np.float64)
        centers = np.asarray(manifest["semantic_codec"]["block_centers_log2_residual"], dtype=np.float64)
        if len(codes["se"]):
            block_rows = int(manifest["value_block_rows"])
            block_centers = centers[target_rows // block_rows]
            ok = codes["se"] < se_count
            safe = np.clip(result["eaf"], 1e-12, 1.0 - 1e-12)
            residual = -1.0 + (codes["se"].astype(np.float64) + 0.5) * (2.0 / se_count)
            se[ok] = 2.0 ** (residual[ok] + block_centers[ok] -
                              0.5 * np.log2(2.0 * safe[ok] * (1.0 - safe[ok])))
        result["se"] = se
    exception_dtype = np.dtype([("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"), ("eaf", "<f4"), ("flags", "u1")])
    exceptions = np.concatenate(exception_parts) if exception_parts else np.empty(0, dtype=exception_dtype)
    rows = exceptions[np.isin(exceptions["row"], target_rows)]
    if len(rows):
        local = np.searchsorted(target_rows, rows["row"].astype(np.int64))
        z_mask = (rows["flags"] & 1) != 0
        se_mask = (rows["flags"] & 2) != 0
        eaf_mask = (rows["flags"] & 4) != 0
        if "z" in result:
            result["z"][local[z_mask]] = rows["z"][z_mask]
        if "se" in result:
            result["se"][local[se_mask]] = np.exp2(
                rows["log2se"][se_mask].astype(np.float64)
            )
        if "eaf" in result:
            result["eaf"][local[eaf_mask]] = rows["eaf"][eaf_mask]
    return {name: result[name] for name in needed}


def _region_global_bounds(manifest: dict[str, Any], chromosome: str, start: int, end: int):
    chromosome = str(chromosome).removeprefix("chr").upper()
    lengths = CHROM_LENGTHS
    offsets = chromosome_offsets()
    if chromosome not in lengths:
        raise ValueError(f"unsupported GRCh38 primary chromosome: {chromosome}")
    if start < 1 or end < start or end > int(lengths[chromosome]):
        raise ValueError(f"invalid region for chromosome {chromosome}: {start}-{end}")
    return int(offsets[chromosome]) + start - 1, int(offsets[chromosome]) + end - 1


def _rows_for_global_region(key_index, global_start: int, global_end: int):
    selected = [frame for frame in key_index
                if frame["last_position"] >= global_start and frame["base_position"] <= global_end]
    if not selected:
        return 0, 0
    return selected[0]["row_start"], selected[-1]["row_stop"]


def parse_identity_keys(values, manifest):
    """Parse canonical ``chromosome:position:REF:ALT`` keys."""
    np, _, _, _ = _dependencies()
    offsets = chromosome_offsets()
    lengths = CHROM_LENGTHS
    parsed = []
    for raw in values:
        value = str(raw).strip()
        fields = value.split(":")
        if len(fields) != 4:
            raise ValueError(
                f"invalid canonical identity key {value!r}; expected chromosome:position:REF:ALT"
            )
        chromosome = fields[0]
        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]
        chromosome = chromosome.upper()
        ref = fields[2].upper()
        alt = fields[3].upper()
        try:
            position = int(fields[1])
        except ValueError as exc:
            raise ValueError(f"invalid position in canonical identity key {value!r}") from exc
        if chromosome not in offsets:
            raise ValueError(f"unsupported chromosome in canonical identity key {value!r}")
        if position < 1 or position > int(lengths[chromosome]):
            raise ValueError(f"position outside GRCh38 in canonical identity key {value!r}")
        if ref not in BASE_CODE or alt not in BASE_CODE or ref == alt:
            raise ValueError(f"invalid REF/ALT substitution in canonical identity key {value!r}")
        global_position = int(offsets[chromosome]) + position - 1
        parsed.append((global_position << 4) | (BASE_CODE[ref] << 2) | BASE_CODE[alt])
    return np.asarray(sorted(set(parsed)), dtype=np.uint64)


def _rows_for_identity_keys(store: Path, manifest: dict[str, Any], key_index, target_keys):
    """Resolve sorted canonical keys to zero-based rows without a whole-key scan."""
    np, _, _, _ = _dependencies()
    if not len(target_keys):
        return np.empty(0, dtype=np.int64)
    target_positions = target_keys >> np.uint64(4)
    hits = []
    for frame in key_index:
        base_position = int(frame["base_position"])
        last_position = int(frame["last_position"])
        lo = int(np.searchsorted(target_positions, base_position, side="left"))
        hi = int(np.searchsorted(target_positions, last_position, side="right"))
        if lo == hi:
            continue
        candidates = target_keys[lo:hi]
        decoded = decode_key_frame(store, manifest, frame)
        locations = np.searchsorted(decoded, candidates)
        valid = locations < len(decoded)
        if np.any(valid):
            valid_locations = locations[valid]
            matched = decoded[valid_locations] == candidates[valid]
            if np.any(matched):
                hits.append(
                    np.asarray(valid_locations[matched] + int(frame["row_start"]), dtype=np.int64)
                )
    if not hits:
        return np.empty(0, dtype=np.int64)
    return np.asarray(sorted(set(np.concatenate(hits).tolist())), dtype=np.int64)


def read_store(
    store: Path,
    output: Path,
    row_start: int | None = None,
    row_stop: int | None = None,
    row_indices=None,
    identity_keys=None,
    chromosome: str | None = None,
    start: int | None = None,
    end: int | None = None,
    include_p: bool = True,
    columns=None,
    output_format: str = "tsv",
):
    np, _, _, _ = _dependencies()
    manifest = load_manifest(store)
    default_columns = [
        "chromosome", "base_pair_location", "effect_allele", "other_allele",
        "z", "beta", "standard_error", "effect_allele_frequency", "p_value",
    ]
    requested_columns = default_columns if columns is None else list(dict.fromkeys(columns))
    if not include_p:
        requested_columns = [name for name in requested_columns if name != "p_value"]
    unknown = set(requested_columns) - set(default_columns)
    if unknown:
        raise ValueError("unknown output column(s): " + ", ".join(sorted(unknown)))
    n = int(manifest["rows"])
    if row_indices is not None and identity_keys is not None:
        raise ValueError("row indices and canonical identity keys are mutually exclusive")
    requested_rows = None if row_indices is None else np.asarray(sorted(set(int(row) for row in row_indices)), dtype=np.int64)
    if requested_rows is not None and np.any((requested_rows < 0) | (requested_rows >= n)):
        raise ValueError("row indices are outside the store")
    key_index = json.loads((store / manifest["files"]["key_index"]).read_text())
    if identity_keys is not None:
        requested_keys = parse_identity_keys(identity_keys, manifest)
        requested_rows = _rows_for_identity_keys(store, manifest, key_index, requested_keys)
    global_region = None
    if chromosome is not None or start is not None or end is not None:
        if chromosome is None or start is None or end is None:
            raise ValueError("chromosome, start, and end must be supplied together")
        global_region = _region_global_bounds(manifest, chromosome, int(start), int(end))
        region_start, region_stop = _rows_for_global_region(key_index, *global_region)
        row_start = region_start if row_start is None else max(row_start, region_start)
        row_stop = region_stop if row_stop is None else min(row_stop, region_stop)
    row_start = 0 if row_start is None else max(0, int(row_start))
    row_stop = n if row_stop is None else min(n, int(row_stop))
    identity_columns = {
        "chromosome", "base_pair_location", "effect_allele", "other_allele"
    }
    need_identity = bool(identity_columns.intersection(requested_columns)) or global_region is not None
    if not need_identity and global_region is None:
        if requested_rows is None:
            output_rows = np.arange(row_start, row_stop, dtype=np.int64)
        else:
            output_rows = requested_rows
        keys = np.empty(0, dtype=np.uint64)
        frames = []
    elif requested_rows is None:
        frames = [frame for frame in key_index if frame["row_stop"] > row_start and frame["row_start"] < row_stop]
    else:
        frames = [frame for frame in key_index if np.any(
            (requested_rows >= frame["row_start"]) & (requested_rows < frame["row_stop"])
        )]
    if need_identity:
        key_parts = []
        row_parts = []
        for frame in frames:
            decoded = decode_key_frame(store, manifest, frame)
            if requested_rows is None:
                slice_start = max(row_start, frame["row_start"]) - frame["row_start"]
                slice_stop = min(row_stop, frame["row_stop"]) - frame["row_start"]
                key_parts.append(decoded[slice_start:slice_stop])
                row_parts.append(np.arange(max(row_start, frame["row_start"]),
                                           min(row_stop, frame["row_stop"]), dtype=np.int64))
            else:
                frame_rows = requested_rows[(requested_rows >= frame["row_start"]) &
                                             (requested_rows < frame["row_stop"])]
                key_parts.append(decoded[frame_rows - frame["row_start"]])
                row_parts.append(frame_rows)
        keys = np.concatenate(key_parts) if key_parts else np.empty(0, dtype=np.uint64)
        output_rows = np.concatenate(row_parts) if row_parts else np.empty(0, dtype=np.int64)
    numeric_needed = set()
    if any(name in requested_columns for name in ("z", "beta", "p_value")):
        numeric_needed.add("z")
    if any(name in requested_columns for name in ("standard_error", "beta")):
        numeric_needed.add("se")
    if "effect_allele_frequency" in requested_columns:
        numeric_needed.add("eaf")
    full_code_bridge = (output_format == "binary" and requested_rows is None and
                        global_region is None and row_start == 0 and row_stop == n and
                        bool(numeric_needed))
    semantic_codes = exceptions = None
    if full_code_bridge:
        semantic_codes, exceptions = decode_full_semantic_codes(
            store, manifest, streams=numeric_needed
        )
        values = {}
    else:
        values = decode_values(
            store, manifest, row_start, row_stop,
            row_indices=output_rows if requested_rows is not None else None,
            needed=numeric_needed,
        )
    if need_identity:
        chromosome_codes, positions, effect_codes, other_codes = _keys_to_columns(keys, manifest)
    else:
        chromosome_codes = np.empty(0, dtype=np.uint8)
        positions = np.empty(0, dtype=np.int32)
        effect_codes = np.empty(0, dtype=np.uint8)
        other_codes = np.empty(0, dtype=np.uint8)
    if global_region is not None:
        chromosome_names = list(CHROM_LENGTHS)
        target_chromosome = str(chromosome).removeprefix("chr").upper()
        target_code = chromosome_names.index(target_chromosome) + 1
        keep = ((chromosome_codes == target_code) &
                (positions >= int(start)) & (positions <= int(end)))
        keys = keys[keep]
        output_rows = output_rows[keep]
        for name in values:
            values[name] = values[name][keep]
        chromosome_codes = chromosome_codes[keep]
        positions, effect_codes, other_codes = (positions[keep], effect_codes[keep],
                                                 other_codes[keep])
    output_values = {
        "chromosome": chromosome_codes,
        "base_pair_location": positions,
        "effect_allele": effect_codes,
        "other_allele": other_codes,
    }
    if "z" in values:
        output_values["z"] = values["z"]
    if "se" in values:
        output_values["standard_error"] = values["se"]
    if "eaf" in values:
        output_values["effect_allele_frequency"] = values["eaf"]
    if semantic_codes is None and "beta" in requested_columns:
        output_values["beta"] = values["z"] * values["se"]
    if semantic_codes is None and "p_value" in requested_columns:
        output_values["p_value"] = np.asarray([
            math.erfc(abs(float(value)) / math.sqrt(2.0))
            if math.isfinite(float(value)) else float("nan")
            for value in values["z"]
        ])
    if output_format == "binary":
        _write_binary_bridge(
            output, output_rows, requested_columns, output_values,
            semantic_codes=semantic_codes, exceptions=exceptions, manifest=manifest,
        )
        return
    if output_format != "tsv":
        raise ValueError(f"unsupported output format: {output_format}")
    chromosome_names = list(CHROM_LENGTHS)
    bases = np.asarray(["A", "C", "G", "T"])
    output_values["chromosome"] = [chromosome_names[index - 1] for index in chromosome_codes]
    output_values["effect_allele"] = bases[effect_codes]
    output_values["other_allele"] = bases[other_codes]
    if "beta" in requested_columns:
        output_values["beta"] = values["z"] * values["se"]
    if "p_value" in requested_columns:
        output_values["p_value"] = np.asarray([
            math.erfc(abs(float(value)) / math.sqrt(2.0)) if math.isfinite(float(value)) else float("nan")
            for value in values["z"]
        ])
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        header = ["row", *requested_columns]
        writer.writerow(header)
        for i in range(len(output_rows)):
            row = [int(output_rows[i])]
            for name in requested_columns:
                value = output_values[name][i]
                row.append(int(value) if name == "base_pair_location" else value)
            writer.writerow(row)


def _validate_index(index, rows: int, label: str):
    if rows and not index:
        raise ValueError(f"{label} index is empty")
    expected = 0
    for frame in index:
        start = int(frame["row_start"])
        stop = int(frame["row_stop"])
        if start != expected or stop <= start or stop > rows:
            raise ValueError(f"{label} index does not cover contiguous valid rows")
        expected = stop
    if expected != rows:
        raise ValueError(f"{label} index ends at {expected}, expected {rows}")


def validate_store(store: Path, full: bool = False):
    np, _, _, _ = _dependencies()
    manifest = load_manifest(store)
    rows = int(manifest["rows"])
    for relative in set(manifest["files"].values()):
        _verify_file(store, manifest, relative)
    key_index = json.loads((store / manifest["files"]["key_index"]).read_text())
    value_index = json.loads((store / manifest["files"]["value_index"]).read_text())
    _validate_index(key_index, rows, "key")
    _validate_index(value_index, rows, "value")
    exception_dtype = np.dtype([("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"), ("eaf", "<f4"), ("flags", "u1")])
    exception_parts = [decode_exception_frame(store, manifest, frame) for frame in value_index]
    exceptions = np.concatenate(exception_parts) if exception_parts else np.empty(0, dtype=exception_dtype)
    if len(exceptions) != int(manifest["semantic_codec"]["exception_rows"]):
        raise ValueError("exception row count differs from the manifest")
    if len(exceptions):
        exception_rows = exceptions["row"].astype(np.int64)
        if np.any(exception_rows < 0) or np.any(exception_rows >= rows):
            raise ValueError("exception row is outside the store")
        if np.any(exception_rows[1:] <= exception_rows[:-1]):
            raise ValueError("exception rows are not strictly increasing")
        if np.any(exceptions["flags"] == 0) or np.any((exceptions["flags"] & ~7) != 0):
            raise ValueError("exception flags are invalid")
    if full:
        previous_key = None
        genome_length = sum(int(value) for value in manifest["identity"]["chromosome_lengths"].values())
        for frame in key_index:
            keys = decode_key_frame(store, manifest, frame)
            expected = int(frame["row_stop"]) - int(frame["row_start"])
            if len(keys) != expected:
                raise ValueError("decoded key frame has the wrong row count")
            if len(keys) > 1 and np.any(keys[1:] <= keys[:-1]):
                raise ValueError("identity keys are not strictly increasing")
            if previous_key is not None and len(keys) and int(keys[0]) <= previous_key:
                raise ValueError("identity keys are not strictly increasing across frames")
            if len(keys):
                previous_key = int(keys[-1])
                positions = keys >> np.uint64(4)
                if int(frame["base_position"]) != int(positions[0]) or \
                   int(frame["last_position"]) != int(positions[-1]):
                    raise ValueError("key frame genomic bounds disagree with decoded positions")
                substitutions = keys & np.uint64(15)
                refs = substitutions >> np.uint64(2)
                alts = substitutions & np.uint64(3)
                if np.any(positions >= genome_length) or np.any(refs == alts):
                    raise ValueError("decoded identity key is outside the GRCh38 SNV contract")
        for frame, frame_exceptions in zip(value_index, exception_parts):
            decoded = decode_value_frame(store, manifest, frame)
            frame_start = int(frame["row_start"])
            frame_stop = int(frame["row_stop"])
            expected = frame_stop - frame_start
            if any(len(decoded[name]) != expected for name in ("z", "eaf", "se")):
                raise ValueError("decoded value frame has the wrong row count")
            if np.any(decoded["z"] > 511) or np.any(decoded["se"] > 63):
                raise ValueError("numeric code exceeds its declared bit width")
            z_flags = np.zeros(expected, dtype=bool)
            se_flags = np.zeros(expected, dtype=bool)
            if len(frame_exceptions):
                local = frame_exceptions["row"].astype(np.int64) - frame_start
                z_flags[local] = (frame_exceptions["flags"] & 1) != 0
                se_flags[local] = (frame_exceptions["flags"] & 2) != 0
            if not np.array_equal(decoded["z"] == 511, z_flags):
                raise ValueError("Z exception flags and sentinel codes disagree")
            if not np.array_equal(decoded["se"] == 63, se_flags):
                raise ValueError("SE exception flags and sentinel codes disagree")
    return {
        "valid": True,
        "errors": [],
        "rows": rows,
        "profile": manifest["profile"],
        "full": bool(full),
        "files_checked": len(set(manifest["files"].values())),
        "frames_checked": len(key_index) + len(value_index) if full else 0,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    write = subparsers.add_parser("write")
    write.add_argument("--input", type=Path, required=True)
    write.add_argument("--output", type=Path, required=True)
    write.add_argument("--metadata", type=Path)
    read = subparsers.add_parser("read")
    read.add_argument("--store", type=Path, required=True)
    read.add_argument("--output", type=Path, required=True)
    read.add_argument("--row-start", type=int)
    read.add_argument("--row-stop", type=int)
    read.add_argument("--rows")
    read.add_argument("--keys-file", type=Path)
    read.add_argument("--chromosome")
    read.add_argument("--start", type=int)
    read.add_argument("--end", type=int)
    read.add_argument("--omit-p", action="store_true")
    read.add_argument("--columns")
    read.add_argument("--output-format", choices=("tsv", "binary"), default="tsv")
    validate = subparsers.add_parser("validate")
    validate.add_argument("--store", type=Path, required=True)
    validate.add_argument("--full", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "write":
        metadata = json.loads(args.metadata.read_text()) if args.metadata else None
        print(json.dumps(write_store(args.input, args.output, metadata=metadata), indent=2))
    elif args.command == "read":
        rows = None if args.rows is None else [int(value) for value in args.rows.split(",") if value]
        keys = None if args.keys_file is None else [
            value.strip() for value in args.keys_file.read_text().splitlines() if value.strip()
        ]
        columns = None if args.columns is None else [value for value in args.columns.split(",") if value]
        read_store(args.store, args.output, row_start=args.row_start, row_stop=args.row_stop,
                   row_indices=rows, identity_keys=keys, chromosome=args.chromosome,
                   start=args.start, end=args.end,
                   include_p=not args.omit_p, columns=columns,
                   output_format=args.output_format)
    else:
        print(json.dumps(validate_store(args.store, full=args.full), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
