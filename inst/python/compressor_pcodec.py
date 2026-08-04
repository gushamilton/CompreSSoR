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
import concurrent.futures
import csv
try:
    import fcntl
except ImportError:  # pragma: no cover - non-POSIX platforms
    fcntl = None
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
import threading
import zlib
from collections import OrderedDict
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
KEY_BLOCK_ROWS = 8_192
VALUE_BLOCK_ROWS = 32_768
SE_CENTER_BLOCK_ROWS = 65_536
MIN_FRAME_ROWS = 4_096
MAX_FRAME_ROWS = 1_048_576
GAP_ESCAPE = 511
GAP_ESCAPE_THRESHOLD = 510
SECOND_GAP_ESCAPE = 65_535
SECOND_GAP_ESCAPE_THRESHOLD = 65_534
FORMAT = "CompreSSoR"
VERSION_V2 = "0.2.0-pcodec"
VERSION = "0.3.0-pcodec"
MANIFEST_CHECKSUM = "manifest.sha256"
MAX_STORE_ROWS = 2**31 - 1
EXPECTED_FILES_V2 = {
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
WRAPPED_CHUNK_ROWS = 262_144
WRAPPED_PAGE_ROWS = 4_096
WRAPPED_PAGE_CACHE_PAGES = 64
WRAPPED_STREAM_DTYPES = {
    "position": ("u32", "<u4", 1),
    "substitution": ("u8", "u1", 2),
    "z": ("u16", "<u2", 3),
    "eaf": ("u8", "u1", 2),
    "se": ("u8", "u1", 2),
}
EXPECTED_FILES_V3 = {
    **{name: f"{name}.pco" for name in WRAPPED_STREAM_DTYPES},
    **{f"{name}_index": f"{name}.index" for name in WRAPPED_STREAM_DTYPES},
    "exceptions": "exceptions.zst",
}
INDEX_MAGIC = b"CPRWIDX1"
INDEX_VERSION = 2
INDEX_FLAG_KEY_BOUNDS = 1
INDEX_HEADER = struct.Struct("<8sBBHQIIIIII")
INDEX_CHUNK = struct.Struct("<QIIII")
INDEX_PAGE = struct.Struct("<QIIQII")
INDEX_KEY_BOUNDS = struct.Struct("<QQ")


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


def _nocache_enabled() -> bool:
    return os.environ.get("COMPRESSOR_NOCACHE", "").strip().lower() in {
        "1", "true", "yes", "on",
    }


def _open_read_fd(path: Path) -> int:
    fd = os.open(path, os.O_RDONLY)
    if _nocache_enabled() and fcntl is not None and hasattr(fcntl, "F_NOCACHE"):
        try:
            fcntl.fcntl(fd, fcntl.F_NOCACHE, 1)
        except Exception:
            os.close(fd)
            raise
    return fd


def _read_path_bytes(path: Path) -> bytes:
    fd = _open_read_fd(path)
    try:
        return _read_fd_bytes(fd, path.stat().st_size, path)
    finally:
        os.close(fd)


def _read_fd_bytes(fd: int, size: int, path: Path) -> bytes:
    """Read one descriptor sequentially in large requests.

    Wrapped pages remain independently decodable, but a complete scan should
    not turn a newly written or fragmented file into hundreds of small random
    reads.  Darwin's F_NOCACHE mode makes that distinction especially large.
    """
    chunks = []
    offset = 0
    while offset < size:
        chunk = _pread(fd, min(8 * 1024 * 1024, size - offset), offset)
        if not chunk:
            raise ValueError(f"could not read complete file: {path}")
        chunks.append(chunk)
        offset += len(chunk)
    return b"".join(chunks)


def _pread(fd: int, size: int, offset: int) -> bytes:
    """Portable positioned read; callers serialize access to each descriptor."""
    if hasattr(os, "pread"):
        return os.pread(fd, size, offset)
    current = os.lseek(fd, 0, os.SEEK_CUR)  # pragma: no cover - Windows
    try:
        os.lseek(fd, offset, os.SEEK_SET)
        return os.read(fd, size)
    finally:
        os.lseek(fd, current, os.SEEK_SET)


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


def _validated_block_rows(value: Any, label: str) -> int:
    try:
        block_rows = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is invalid") from exc
    if (block_rows < MIN_FRAME_ROWS or block_rows > MAX_FRAME_ROWS or
            block_rows & (block_rows - 1)):
        raise ValueError(
            f"{label} must be a power of two from "
            f"{MIN_FRAME_ROWS} through {MAX_FRAME_ROWS}"
        )
    return block_rows


def _validate_manifest_contract_v2(manifest: dict[str, Any]) -> None:
    """Reject altered decode-critical constants for the fixed v0.2 format."""
    if manifest.get("format") != FORMAT:
        raise ValueError("not a CompreSSoR store")
    if manifest.get("format_version") != VERSION_V2:
        raise ValueError(
            f"unsupported Pcodec format version: {manifest.get('format_version')!r}; expected {VERSION_V2!r}"
        )
    if manifest.get("backend") != "pcodec" or manifest.get("profile") != "standard":
        raise ValueError("store is not a standard Pcodec store")
    try:
        rows = int(manifest["rows"])
        n_rows = int(manifest["n_rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("manifest row count is missing or invalid") from exc
    if rows < 1 or rows > MAX_STORE_ROWS or n_rows != rows:
        raise ValueError("manifest row counts disagree or are empty")
    key_block_rows = _validated_block_rows(
        manifest.get("key_block_rows"), "manifest key frame size"
    )
    value_block_rows = _validated_block_rows(
        manifest.get("value_block_rows"), "manifest value frame size"
    )

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
    center_block_rows = _validated_block_rows(
        semantic.get("se_center_block_rows", value_block_rows),
        "manifest SE centre block size",
    )
    centers = semantic.get("block_centers_log2_residual")
    expected_centers = (rows + center_block_rows - 1) // center_block_rows
    if not isinstance(centers, list) or len(centers) != expected_centers:
        raise ValueError("manifest has the wrong number of SE block centres")
    try:
        finite_centers = all(math.isfinite(float(value)) for value in centers)
        exception_rows = int(semantic["exception_rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("manifest semantic metadata is invalid") from exc
    if not finite_centers or exception_rows < 0 or exception_rows > rows:
        raise ValueError("manifest semantic metadata is outside the format contract")
    if manifest.get("files") != EXPECTED_FILES_V2:
        raise ValueError("manifest file layout differs from the fixed Pcodec layout")
    integrity = manifest.get("integrity", {})
    records = integrity.get("files", {})
    if integrity.get("algorithm") != "sha256" or set(records) != set(EXPECTED_FILES_V2.values()):
        raise ValueError("manifest integrity table is incomplete")
    for relative in EXPECTED_FILES_V2.values():
        record = records.get(relative, {})
        digest = record.get("sha256")
        try:
            size = int(record.get("bytes", -1))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"manifest integrity size is invalid for {relative}") from exc
        if size < 0 or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"manifest integrity record is invalid for {relative}")


def _validate_manifest_contract_v3(manifest: dict[str, Any]) -> None:
    if manifest.get("format") != FORMAT or manifest.get("format_version") != VERSION:
        raise ValueError("store is not a supported wrapped Pcodec format")
    if manifest.get("backend") != "pcodec" or manifest.get("profile") != "standard":
        raise ValueError("store is not a standard Pcodec store")
    try:
        rows = int(manifest["rows"])
        n_rows = int(manifest["n_rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("manifest row count is missing or invalid") from exc
    if rows < 1 or rows > MAX_STORE_ROWS or n_rows != rows:
        raise ValueError("manifest row counts disagree or are empty")
    wrapped = manifest.get("wrapped_codec", {})
    if not isinstance(wrapped, dict):
        raise ValueError("manifest wrapped-codec metadata must be an object")
    try:
        wrapped_constants_valid = (
            int(wrapped.get("chunk_rows", 0)) == WRAPPED_CHUNK_ROWS and
            int(wrapped.get("page_rows", 0)) == WRAPPED_PAGE_ROWS and
            wrapped.get("index_version") == INDEX_VERSION and
            wrapped.get("integrity") == "CRC32 per page/meta/header; SHA256 per file"
        )
        page_geometry_valid = (
            int(manifest.get("key_block_rows", 0)) == WRAPPED_PAGE_ROWS and
            int(manifest.get("value_block_rows", 0)) == WRAPPED_PAGE_ROWS
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest wrapped-codec constants are invalid") from exc
    if not wrapped_constants_valid:
        raise ValueError("manifest wrapped-codec constants differ from the format contract")
    if not page_geometry_valid:
        raise ValueError("manifest page geometry differs from the format contract")
    identity = manifest.get("identity", {})
    expected_identity = {
        "encoding": "wrapped_global_position_plus_full_ref_alt_code",
        "external_reference_required": False,
        "chromosome_lengths": CHROM_LENGTHS,
        "chromosome_offsets": chromosome_offsets(),
        "effect_allele_is_alt": True,
        "other_allele_is_ref": True,
    }
    if identity != expected_identity:
        raise ValueError("manifest identity constants differ from the wrapped GRCh38 contract")
    semantic = manifest.get("semantic_codec", {})
    if not isinstance(semantic, dict):
        raise ValueError("manifest semantic metadata must be an object")
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
    center_rows = _validated_block_rows(
        semantic.get("se_center_block_rows"), "manifest SE centre block size"
    )
    centers = semantic.get("block_centers_log2_residual")
    try:
        finite_centers = isinstance(centers, list) and all(
            math.isfinite(float(value)) for value in centers
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("manifest has invalid SE block centres") from exc
    if (not isinstance(centers, list) or
            len(centers) != (rows + center_rows - 1) // center_rows or
            not finite_centers):
        raise ValueError("manifest has invalid SE block centres")
    try:
        exception_rows = int(semantic["exception_rows"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("manifest exception count is invalid") from exc
    if exception_rows < 0 or exception_rows > rows:
        raise ValueError("manifest exception count is outside the format contract")
    if manifest.get("files") != EXPECTED_FILES_V3:
        raise ValueError("manifest file layout differs from the wrapped format")
    integrity = manifest.get("integrity", {})
    if not isinstance(integrity, dict):
        raise ValueError("manifest integrity metadata must be an object")
    records = integrity.get("files", {})
    if (not isinstance(records, dict) or integrity.get("algorithm") != "sha256" or
            set(records) != set(EXPECTED_FILES_V3.values())):
        raise ValueError("manifest integrity table is incomplete")
    for relative in EXPECTED_FILES_V3.values():
        record = records.get(relative, {})
        if not isinstance(record, dict):
            raise ValueError(f"manifest integrity record is invalid for {relative}")
        try:
            size = int(record.get("bytes", -1))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"manifest integrity record is invalid for {relative}") from exc
        if (not isinstance(record.get("sha256"), str) or
                not re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) or size < 0):
            raise ValueError(f"manifest integrity record is invalid for {relative}")


def _validate_manifest_contract(manifest: dict[str, Any]) -> None:
    version = manifest.get("format_version")
    if version == VERSION_V2:
        _validate_manifest_contract_v2(manifest)
    elif version == VERSION:
        _validate_manifest_contract_v3(manifest)
    else:
        raise ValueError(f"unsupported Pcodec format version: {version!r}")


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


def quantise_values(data, block_rows: int = SE_CENTER_BLOCK_ROWS):
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


def write_store_v2(
    input_path: Path,
    output: Path,
    metadata: dict[str, Any] | None = None,
    key_block_rows: int = KEY_BLOCK_ROWS,
    value_block_rows: int = VALUE_BLOCK_ROWS,
    se_center_block_rows: int = SE_CENTER_BLOCK_ROWS,
):
    np, _, _, _ = _dependencies()
    key_block_rows = _validated_block_rows(key_block_rows, "key frame size")
    value_block_rows = _validated_block_rows(value_block_rows, "value frame size")
    se_center_block_rows = _validated_block_rows(
        se_center_block_rows, "SE centre block size"
    )
    output.mkdir(parents=True, exist_ok=True)
    data = read_tsv(input_path)
    keys, order = identity_keys(data)
    if len(keys) > MAX_STORE_ROWS:
        raise ValueError(f"wrapped Pcodec stores support at most {MAX_STORE_ROWS} rows")
    ordered = {name: values[order] for name, values in data.items()}
    values = quantise_values(ordered, block_rows=se_center_block_rows)
    key_index, key_paths = encode_key_frames(keys, output, block_rows=key_block_rows)
    value_index, value_paths = encode_value_frames(
        values, output, block_rows=value_block_rows
    )
    manifest = {
        "format": FORMAT,
        "format_version": VERSION_V2,
        "backend": "pcodec",
        "profile": "standard",
        "rows": len(keys),
        "n_rows": len(keys),
        "key_block_rows": key_block_rows,
        "value_block_rows": value_block_rows,
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
            "se_center_block_rows": se_center_block_rows,
            "block_centers_log2_residual": values["centers"],
            "exception_rows": int(len(values["exceptions"])),
            "exception_precision": "float32",
            "beta": "derived as z * standard_error",
            "p_value": "derived as erfc(abs(z) / sqrt(2))",
        },
        "files": EXPECTED_FILES_V2,
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


def _write_wrapped_stream(output: Path, name: str, values, identity_keys=None):
    np, _, ChunkConfig, _ = _dependencies()
    from pcodec import PagingSpec
    from pcodec.wrapped import FileCompressor

    pco_dtype, numpy_dtype, dtype_code = WRAPPED_STREAM_DTYPES[name]
    values = np.ascontiguousarray(values, dtype=numpy_dtype)
    data_path = output / EXPECTED_FILES_V3[name]
    index_path = output / EXPECTED_FILES_V3[f"{name}_index"]
    compressor = FileCompressor()
    header = compressor.write_header()
    config = ChunkConfig(
        enable_8_bit=True,
        paging_spec=PagingSpec.equal_pages_up_to(WRAPPED_PAGE_ROWS),
    )
    chunks = []
    pages = []
    with data_path.open("wb") as handle:
        handle.write(header)
        for chunk_start in range(0, len(values), WRAPPED_CHUNK_ROWS):
            chunk_stop = min(chunk_start + WRAPPED_CHUNK_ROWS, len(values))
            chunk = compressor.chunk_compressor(values[chunk_start:chunk_stop], config)
            metadata = chunk.write_meta()
            metadata_offset = handle.tell()
            handle.write(metadata)
            first_page = len(pages)
            local_start = 0
            counts = chunk.n_per_page()
            for page_number, count in enumerate(counts):
                blob = chunk.write_page(page_number)
                offset = handle.tell()
                handle.write(blob)
                row_start = chunk_start + local_start
                row_stop = row_start + int(count)
                bounds = None
                if identity_keys is not None:
                    bounds = (
                        int(identity_keys[row_start]),
                        int(identity_keys[row_stop - 1]),
                    )
                pages.append((
                    offset, len(blob), zlib.crc32(blob), row_start,
                    int(count), len(chunks), bounds,
                ))
                local_start += int(count)
            chunks.append((
                metadata_offset, len(metadata), zlib.crc32(metadata),
                first_page, len(counts),
            ))
    flags = INDEX_FLAG_KEY_BOUNDS if identity_keys is not None else 0
    with index_path.open("wb") as handle:
        handle.write(INDEX_HEADER.pack(
            INDEX_MAGIC, INDEX_VERSION, dtype_code, flags, len(values),
            WRAPPED_CHUNK_ROWS, WRAPPED_PAGE_ROWS, len(header),
            zlib.crc32(header),
            len(chunks), len(pages),
        ))
        for record in chunks:
            handle.write(INDEX_CHUNK.pack(*record))
        for offset, length, crc, row_start, count, chunk_id, bounds in pages:
            handle.write(INDEX_PAGE.pack(
                offset, length, crc, row_start, count, chunk_id
            ))
            if bounds is not None:
                handle.write(INDEX_KEY_BOUNDS.pack(*bounds))
    return {
        "dtype": pco_dtype,
        "chunks": len(chunks),
        "pages": len(pages),
        "bytes": data_path.stat().st_size + index_path.stat().st_size,
    }


def write_store(
    input_path: Path,
    output: Path,
    metadata: dict[str, Any] | None = None,
    key_block_rows: int = WRAPPED_PAGE_ROWS,
    value_block_rows: int = WRAPPED_PAGE_ROWS,
    se_center_block_rows: int = SE_CENTER_BLOCK_ROWS,
):
    """Write the wrapped-page v0.3 format; v0.2 remains readable."""
    np, zstd, _, _ = _dependencies()
    if int(key_block_rows) != WRAPPED_PAGE_ROWS or int(value_block_rows) != WRAPPED_PAGE_ROWS:
        raise ValueError(
            f"wrapped Pcodec stores require {WRAPPED_PAGE_ROWS}-row pages"
        )
    se_center_block_rows = _validated_block_rows(
        se_center_block_rows, "SE centre block size"
    )
    output.mkdir(parents=True, exist_ok=True)
    data = read_tsv(input_path)
    keys, order = identity_keys(data)
    ordered = {name: column[order] for name, column in data.items()}
    values = quantise_values(ordered, block_rows=se_center_block_rows)
    positions = (keys >> np.uint64(4)).astype(np.uint32)
    substitutions = (keys & np.uint64(15)).astype(np.uint8)
    stream_values = {
        "position": positions,
        "substitution": substitutions,
        "z": values["z"],
        "eaf": values["eaf"],
        "se": values["se"],
    }
    stream_metadata = {
        name: _write_wrapped_stream(
            output, name, column,
            identity_keys=keys if name == "position" else None,
        )
        for name, column in stream_values.items()
    }
    exception_path = output / EXPECTED_FILES_V3["exceptions"]
    exception_path.write_bytes(
        zstd.ZstdCompressor(level=19).compress(values["exceptions"].tobytes())
    )
    manifest = {
        "format": FORMAT,
        "format_version": VERSION,
        "backend": "pcodec",
        "profile": "standard",
        "rows": len(keys),
        "n_rows": len(keys),
        "key_block_rows": WRAPPED_PAGE_ROWS,
        "value_block_rows": WRAPPED_PAGE_ROWS,
        "wrapped_codec": {
            "chunk_rows": WRAPPED_CHUNK_ROWS,
            "page_rows": WRAPPED_PAGE_ROWS,
            "index_version": INDEX_VERSION,
            "integrity": "CRC32 per page/meta/header; SHA256 per file",
            "streams": stream_metadata,
        },
        "identity": {
            "encoding": "wrapped_global_position_plus_full_ref_alt_code",
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
            "se_center_block_rows": se_center_block_rows,
            "block_centers_log2_residual": values["centers"],
            "exception_rows": int(len(values["exceptions"])),
            "exception_precision": "float32",
            "beta": "derived as z * standard_error",
            "p_value": "derived as erfc(abs(z) / sqrt(2))",
        },
        "files": EXPECTED_FILES_V3,
        "input_columns": [
            "chromosome", "base_pair_location", "effect_allele",
            "other_allele", "beta", "standard_error",
            "effect_allele_frequency", "z",
        ],
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pcodec": importlib.metadata.version("pcodec"),
            "zstandard": importlib.metadata.version("zstandard"),
        },
        "metadata": metadata or {},
    }
    integrity = {}
    for relative in manifest["files"].values():
        path = output / relative
        blob = path.read_bytes()
        integrity[relative] = {
            "bytes": len(blob),
            "sha256": hashlib.sha256(blob).hexdigest(),
        }
    manifest["integrity"] = {"algorithm": "sha256", "files": integrity}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    _seal_manifest(output)
    return manifest


class WrappedStreamReader:
    def __init__(self, store: Path, manifest: dict[str, Any], name: str):
        np, _, _, _ = _dependencies()
        from pcodec.wrapped import FileDecompressor

        self.store = store
        self.name = name
        self.pco_dtype, numpy_dtype, expected_dtype_code = WRAPPED_STREAM_DTYPES[name]
        self.numpy_dtype = np.dtype(numpy_dtype)
        self.data_path = store / manifest["files"][name]
        index_path = store / manifest["files"][f"{name}_index"]
        blob = _read_path_bytes(index_path)
        self.bytes_read = len(blob)
        if len(blob) < INDEX_HEADER.size:
            raise ValueError(f"wrapped {name} index is truncated")
        fields = INDEX_HEADER.unpack_from(blob)
        (magic, version, dtype_code, flags, rows, chunk_rows, page_rows,
         header_length, header_crc, chunk_count, page_count) = fields
        if (magic != INDEX_MAGIC or version != INDEX_VERSION or
                dtype_code != expected_dtype_code or rows != int(manifest["rows"]) or
                chunk_rows != WRAPPED_CHUNK_ROWS or page_rows != WRAPPED_PAGE_ROWS):
            raise ValueError(f"wrapped {name} index header is invalid")
        self.flags = flags
        self.rows = int(rows)
        self.header_length = int(header_length)
        offset = INDEX_HEADER.size
        self.chunks = []
        for _ in range(chunk_count):
            if offset + INDEX_CHUNK.size > len(blob):
                raise ValueError(f"wrapped {name} chunk index is truncated")
            self.chunks.append(INDEX_CHUNK.unpack_from(blob, offset))
            offset += INDEX_CHUNK.size
        self.pages = []
        bounds = []
        for _ in range(page_count):
            if offset + INDEX_PAGE.size > len(blob):
                raise ValueError(f"wrapped {name} page index is truncated")
            self.pages.append(INDEX_PAGE.unpack_from(blob, offset))
            offset += INDEX_PAGE.size
            if flags & INDEX_FLAG_KEY_BOUNDS:
                if offset + INDEX_KEY_BOUNDS.size > len(blob):
                    raise ValueError(f"wrapped {name} key bounds are truncated")
                bounds.append(INDEX_KEY_BOUNDS.unpack_from(blob, offset))
                offset += INDEX_KEY_BOUNDS.size
        if offset != len(blob):
            raise ValueError(f"wrapped {name} index has trailing bytes")
        expected_chunks = (self.rows + WRAPPED_CHUNK_ROWS - 1) // WRAPPED_CHUNK_ROWS
        expected_pages = (self.rows + WRAPPED_PAGE_ROWS - 1) // WRAPPED_PAGE_ROWS
        if len(self.chunks) != expected_chunks or len(self.pages) != expected_pages:
            raise ValueError(f"wrapped {name} index has invalid chunk/page counts")
        self.page_starts = np.asarray([page[3] for page in self.pages], dtype=np.int64)
        self.page_stops = self.page_starts + np.asarray(
            [page[4] for page in self.pages], dtype=np.int64
        )
        if (not len(self.pages) or self.page_starts[0] != 0 or
                self.page_stops[-1] != self.rows or
                np.any(self.page_starts[1:] != self.page_stops[:-1])):
            raise ValueError(f"wrapped {name} pages do not cover every row")
        self.first_keys = np.asarray([item[0] for item in bounds], dtype=np.uint64)
        self.last_keys = np.asarray([item[1] for item in bounds], dtype=np.uint64)
        data_size = self.data_path.stat().st_size
        cursor = int(header_length)
        page_cursor = 0
        for chunk_id, chunk in enumerate(self.chunks):
            metadata_offset, metadata_length, _, first_page, pages_in_chunk = chunk
            if (metadata_length <= 0 or metadata_offset != cursor or
                    first_page != page_cursor or pages_in_chunk <= 0 or
                    page_cursor + pages_in_chunk > len(self.pages)):
                raise ValueError(f"wrapped {name} chunk layout is invalid")
            cursor += int(metadata_length)
            for page_id in range(page_cursor, page_cursor + pages_in_chunk):
                page_offset, page_length, _, _, count, owning_chunk = self.pages[page_id]
                # Constant Pcodec pages legitimately have a zero-byte payload;
                # their values are represented entirely by chunk metadata.
                if (page_offset != cursor or count <= 0 or
                        count > WRAPPED_PAGE_ROWS or owning_chunk != chunk_id):
                    raise ValueError(f"wrapped {name} page layout is invalid")
                cursor += int(page_length)
            page_cursor += int(pages_in_chunk)
        if page_cursor != len(self.pages) or cursor != data_size:
            raise ValueError(f"wrapped {name} payload layout or length is invalid")
        self.fd = None
        try:
            self.fd = _open_read_fd(self.data_path)
            header = _pread(self.fd, self.header_length, 0)
            self.bytes_read += len(header)
            if (len(header) != self.header_length or
                    zlib.crc32(header) != int(header_crc)):
                raise ValueError(f"wrapped {name} header checksum mismatch")
            self.file_decompressor, consumed = FileDecompressor.new(header)
            if consumed != self.header_length:
                raise ValueError(f"wrapped {name} header length is invalid")
        except Exception:
            self.close()
            raise
        self.chunk_decompressors = {}
        try:
            cache_pages = int(os.environ.get(
                "COMPRESSOR_PAGE_CACHE_PAGES", WRAPPED_PAGE_CACHE_PAGES
            ))
        except ValueError:
            cache_pages = WRAPPED_PAGE_CACHE_PAGES
        self.page_cache_limit = max(0, min(4096, cache_pages))
        self.page_cache = OrderedDict()

    def close(self):
        if getattr(self, "fd", None) is not None:
            os.close(self.fd)
            self.fd = None
        if hasattr(self, "page_cache"):
            self.page_cache.clear()

    def _read_blob(self, offset: int, length: int, crc: int) -> bytes:
        blob = _pread(self.fd, int(length), int(offset))
        self.bytes_read += len(blob)
        if len(blob) != int(length) or zlib.crc32(blob) != int(crc):
            raise ValueError(f"wrapped {self.name} page or metadata is corrupt")
        return blob

    def _chunk_decompressor(self, chunk_id: int):
        if chunk_id not in self.chunk_decompressors:
            offset, length, crc, _, _ = self.chunks[chunk_id]
            metadata = self._read_blob(offset, length, crc)
            decompressor, consumed = self.file_decompressor.chunk_decompressor(
                metadata, self.pco_dtype
            )
            if consumed != length:
                raise ValueError(f"wrapped {self.name} chunk metadata length is invalid")
            self.chunk_decompressors[chunk_id] = decompressor
        return self.chunk_decompressors[chunk_id]

    def read_page(self, page_id: int):
        np, _, _, _ = _dependencies()
        page_id = int(page_id)
        cached = self.page_cache.pop(page_id, None)
        if cached is not None:
            self.page_cache[page_id] = cached
            return cached
        offset, length, crc, _, count, chunk_id = self.pages[page_id]
        destination = np.empty(int(count), dtype=self.numpy_dtype)
        progress, consumed = self._chunk_decompressor(int(chunk_id)).read_page_into(
            self._read_blob(offset, length, crc), int(count), destination
        )
        if consumed != length or not progress.finished or progress.n_processed != count:
            raise ValueError(f"wrapped {self.name} page did not decode completely")
        if self.page_cache_limit:
            self.page_cache[page_id] = destination
            while len(self.page_cache) > self.page_cache_limit:
                self.page_cache.popitem(last=False)
        return destination

    def read_pages(self, page_ids):
        """Decode several pages using an adaptive physical access plan."""
        np, _, _, _ = _dependencies()
        page_ids = [int(page_id) for page_id in page_ids]
        if not page_ids:
            return []
        if any(
            page_id < 0 or page_id >= len(self.pages)
            for page_id in page_ids
        ):
            raise ValueError(f"wrapped {self.name} page ID is outside the stream")
        if len(set(page_ids)) != len(page_ids):
            raise ValueError(f"wrapped {self.name} page IDs must be unique")

        first_offset = min(int(self.pages[page_id][0]) for page_id in page_ids)
        last_offset = max(
            int(self.pages[page_id][0] + self.pages[page_id][1])
            for page_id in page_ids
        )
        try:
            threshold = int(os.environ.get(
                "COMPRESSOR_SEQUENTIAL_PAGE_THRESHOLD", "0"
            ))
        except ValueError:
            threshold = 0
        sequential = (
            _nocache_enabled() and threshold >= 2 and
            len(page_ids) >= threshold and
            last_offset - first_offset >= self.data_path.stat().st_size // 2
        )
        if not sequential:
            return [self.read_page(page_id) for page_id in page_ids]

        # Genome-scattered instruments can touch many tiny pages over most of
        # a stream.  On a true cold read, one sequential scan is much cheaper
        # than many seeks, while only the selected Pcodec pages are decoded.
        file_blob = _read_fd_bytes(
            self.fd, self.data_path.stat().st_size, self.data_path
        )
        self.bytes_read += len(file_blob)
        chunk_decompressors = {}
        decoded_pages = []
        for page_id in page_ids:
            cached = self.page_cache.pop(page_id, None)
            if cached is not None:
                self.page_cache[page_id] = cached
                decoded_pages.append(cached)
                continue
            offset, length, crc, _, count, chunk_id = self.pages[page_id]
            chunk_id = int(chunk_id)
            decompressor = chunk_decompressors.get(chunk_id)
            if decompressor is None:
                metadata_offset, metadata_length, metadata_crc, _, _ = (
                    self.chunks[chunk_id]
                )
                metadata = file_blob[
                    int(metadata_offset):int(metadata_offset + metadata_length)
                ]
                if (len(metadata) != int(metadata_length) or
                        zlib.crc32(metadata) != int(metadata_crc)):
                    raise ValueError(
                        f"wrapped {self.name} page or metadata is corrupt"
                    )
                decompressor, consumed = self.file_decompressor.chunk_decompressor(
                    metadata, self.pco_dtype
                )
                if consumed != int(metadata_length):
                    raise ValueError(
                        f"wrapped {self.name} chunk metadata length is invalid"
                    )
                chunk_decompressors[chunk_id] = decompressor
            page_blob = file_blob[int(offset):int(offset + length)]
            if (len(page_blob) != int(length) or
                    zlib.crc32(page_blob) != int(crc)):
                raise ValueError(
                    f"wrapped {self.name} page or metadata is corrupt"
                )
            destination = np.empty(int(count), dtype=self.numpy_dtype)
            progress, consumed = decompressor.read_page_into(
                page_blob, int(count), destination
            )
            if (consumed != int(length) or not progress.finished or
                    progress.n_processed != int(count)):
                raise ValueError(
                    f"wrapped {self.name} page did not decode completely"
                )
            decoded_pages.append(destination)
        return decoded_pages

    def read_rows(self, rows=None):
        np, _, _, _ = _dependencies()
        if rows is None:
            # Full scans use one sequential read of the physical stream.  Page
            # framing is still validated and decoded independently, but this
            # avoids a pread for every metadata/page blob on cold or freshly
            # staged APFS files.
            file_blob = _read_fd_bytes(
                self.fd, self.data_path.stat().st_size, self.data_path
            )
            self.bytes_read += len(file_blob)
            output = np.empty(self.rows, dtype=self.numpy_dtype)
            active_chunk = None
            active_decompressor = None
            for page_id in range(len(self.pages)):
                start = int(self.page_starts[page_id])
                stop = int(self.page_stops[page_id])
                cached = self.page_cache.pop(page_id, None)
                if cached is not None:
                    self.page_cache[page_id] = cached
                    output[start:stop] = cached
                    continue
                offset, length, crc, _, count, chunk_id = self.pages[page_id]
                if active_chunk != int(chunk_id):
                    metadata_offset, metadata_length, metadata_crc, _, _ = (
                        self.chunks[int(chunk_id)]
                    )
                    metadata = file_blob[
                        int(metadata_offset):int(metadata_offset + metadata_length)
                    ]
                    if (len(metadata) != int(metadata_length) or
                            zlib.crc32(metadata) != int(metadata_crc)):
                        raise ValueError(
                            f"wrapped {self.name} page or metadata is corrupt"
                        )
                    active_decompressor, consumed = (
                        self.file_decompressor.chunk_decompressor(
                            metadata, self.pco_dtype
                        )
                    )
                    if consumed != int(metadata_length):
                        raise ValueError(
                            f"wrapped {self.name} chunk metadata length is invalid"
                        )
                    active_chunk = int(chunk_id)
                page_blob = file_blob[int(offset):int(offset + length)]
                if (len(page_blob) != int(length) or
                        zlib.crc32(page_blob) != int(crc)):
                    raise ValueError(
                        f"wrapped {self.name} page or metadata is corrupt"
                    )
                destination = output[start:stop]
                progress, consumed = active_decompressor.read_page_into(
                    page_blob, int(count), destination
                )
                if (consumed != length or not progress.finished or
                        progress.n_processed != count):
                    raise ValueError(
                        f"wrapped {self.name} page did not decode completely"
                    )
            return output
        rows = np.asarray(rows, dtype=np.int64)
        if not len(rows):
            return np.empty(0, dtype=self.numpy_dtype)
        if len(rows) > 1 and np.any(rows[1:] <= rows[:-1]):
            raise ValueError("wrapped stream rows must be sorted and unique")

        # Whole-genome and region reads are contiguous.  Avoid constructing a
        # page id for every row and, especially, avoid repeatedly scanning that
        # vector once per page.  Decode each intersecting page exactly once and
        # trim only the boundary pages.
        first_row = int(rows[0])
        last_row = int(rows[-1])
        if len(rows) == last_row - first_row + 1:
            first_page = int(np.searchsorted(self.page_stops, first_row, side="right"))
            last_page = int(np.searchsorted(self.page_stops, last_row, side="right"))
            decoded = np.concatenate(self.read_pages(
                range(first_page, last_page + 1)
            ))
            offset = first_row - int(self.page_starts[first_page])
            return decoded[offset:offset + len(rows)]

        # Sparse rows are sorted by every caller.  Their page ids are therefore
        # sorted too, so contiguous groups can be sliced directly.  The former
        # flatnonzero(page_ids == page_id) loop was O(rows * pages) and became
        # pathological for widely scattered selections.
        page_ids = np.searchsorted(self.page_stops, rows, side="right")
        output = np.empty(len(rows), dtype=self.numpy_dtype)
        group_starts = np.flatnonzero(np.r_[True, page_ids[1:] != page_ids[:-1]])
        group_stops = np.r_[group_starts[1:], len(rows)]
        unique_page_ids = [int(page_ids[index]) for index in group_starts]
        decoded_pages = self.read_pages(unique_page_ids)
        for group_start, group_stop, page_id, page in zip(
                group_starts, group_stops, unique_page_ids, decoded_pages):
            selected = slice(int(group_start), int(group_stop))
            output[selected] = page[rows[selected] - self.page_starts[page_id]]
        return output


class WrappedStoreReader:
    def __init__(self, store: Path, manifest: dict[str, Any]):
        self.store = store
        self.manifest = manifest
        self.streams = {}
        self._exceptions = None
        self.exception_bytes_read = 0
        # Pcodec chunk-decompressor objects are cached and are not documented
        # as thread-safe.  Batch reads may share one store context, so serialize
        # operations on that context while still allowing different stores to
        # decode concurrently.
        self.lock = threading.RLock()

    def stream(self, name: str) -> WrappedStreamReader:
        if name not in self.streams:
            self.streams[name] = WrappedStreamReader(self.store, self.manifest, name)
        return self.streams[name]

    def exceptions(self):
        if self._exceptions is None:
            np, zstd, _, _ = _dependencies()
            relative = self.manifest["files"]["exceptions"]
            path = self.store / relative
            record = self.manifest["integrity"]["files"][relative]
            blob = _read_path_bytes(path)
            self.exception_bytes_read += len(blob)
            if (len(blob) != int(record["bytes"]) or
                    hashlib.sha256(blob).hexdigest() != record["sha256"]):
                raise ValueError("wrapped exception stream checksum mismatch")
            raw = zstd.ZstdDecompressor().decompress(blob)
            dtype = np.dtype([
                ("row", "<u4"), ("z", "<f4"), ("log2se", "<f4"),
                ("eaf", "<f4"), ("flags", "u1"),
            ])
            if len(raw) % dtype.itemsize:
                raise ValueError("wrapped exception stream has a partial record")
            self._exceptions = np.frombuffer(raw, dtype=dtype)
            if len(self._exceptions) != int(self.manifest["semantic_codec"]["exception_rows"]):
                raise ValueError("wrapped exception count differs from the manifest")
            if len(self._exceptions):
                rows = self._exceptions["row"].astype(np.int64)
                flags = self._exceptions["flags"].astype(np.uint16)
                if (np.any(rows < 0) or np.any(rows >= int(self.manifest["rows"])) or
                        np.any(rows[1:] <= rows[:-1])):
                    raise ValueError("wrapped exception rows are invalid")
                if np.any(flags == 0) or np.any((flags & np.uint16(0xF8)) != 0):
                    raise ValueError("wrapped exception flags are invalid")
        return self._exceptions

    @property
    def bytes_read(self):
        return self.exception_bytes_read + sum(
            stream.bytes_read for stream in self.streams.values()
        )

    def close(self):
        for stream in self.streams.values():
            stream.close()
        self.streams.clear()


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
    digest = hashlib.sha256(_read_path_bytes(path)).hexdigest()
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
        if manifest["format_version"] == VERSION_V2:
            _verify_file(store, manifest, manifest["files"]["key_index"])
            _verify_file(store, manifest, manifest["files"]["value_index"])
        else:
            for name in WRAPPED_STREAM_DTYPES:
                _verify_file(store, manifest, manifest["files"][f"{name}_index"])
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
    identity_codes=None,
    exceptions=None,
    manifest=None,
    source_bytes_read=None,
):
    """Write a temporary columnar bridge that R can load with readBin()."""
    np, _, _, _ = _dependencies()
    output_count = int(output_rows) if isinstance(output_rows, (int, np.integer)) else len(output_rows)
    output.mkdir(parents=True, exist_ok=False)
    files = {}

    def write_column(name, values, dtype, dtype_name):
        filename = f"{name}.bin"
        array = np.ascontiguousarray(values, dtype=dtype)
        array.tofile(output / filename)
        files[name] = {"file": filename, "dtype": dtype_name, "length": int(len(array))}

    identity = None
    if identity_codes is not None:
        write_column(
            "global_position_code", identity_codes["position"], "<u4", "uint32"
        )
        write_column(
            "substitution_code", identity_codes["substitution"], "u1", "uint8"
        )
        lengths = manifest["identity"]["chromosome_lengths"]
        identity = {
            "encoding": "global_position_substitution",
            "chromosome_lengths": list(lengths.values()),
        }
    else:
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
            "block_rows": int(
                semantic.get("se_center_block_rows", manifest["value_block_rows"])
            ),
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
        "rows": output_count,
        "requested_columns": requested_columns,
        "files": files,
        "codec": codec,
        "identity": identity,
        "source_bytes_read": source_bytes_read,
    }, indent=2))


def decode_full_semantic_codes(
    store: Path, manifest: dict[str, Any], streams, value_index=None
):
    np, _, _, _ = _dependencies()
    streams = set(streams)
    if "se" in streams:
        streams.add("eaf")
    if value_index is None:
        value_index = json.loads(
            (store / manifest["files"]["value_index"]).read_text()
        )
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
        if (np.any(decoded["flags"] == 0) or
                np.any((decoded["flags"].astype(np.uint16) & np.uint16(0xF8)) != 0)):
            raise ValueError("exception flags are invalid")
    return decoded


def decode_values(
    store: Path,
    manifest: dict[str, Any],
    row_start: int,
    row_stop: int,
    row_indices=None,
    needed=("z", "se", "eaf"),
    value_index=None,
):
    np, _, _, _ = _dependencies()
    needed = set(needed)
    if not needed:
        return {}
    streams = set(needed)
    if "se" in streams:
        streams.add("eaf")
    if value_index is None:
        value_index = json.loads(
            (store / manifest["files"]["value_index"]).read_text()
        )
    if row_indices is None:
        target_rows = np.arange(row_start, row_stop, dtype=np.int64)
    else:
        target_rows = np.asarray(sorted(set(int(row) for row in row_indices)), dtype=np.int64)
    selected = _frames_containing_rows(value_index, target_rows)
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
            block_rows = int(
                manifest["semantic_codec"].get(
                    "se_center_block_rows", manifest["value_block_rows"]
                )
            )
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


def _rows_for_identity_keys(
    store: Path,
    manifest: dict[str, Any],
    key_index,
    target_keys,
    return_keys: bool = False,
):
    """Resolve sorted canonical keys to zero-based rows without a whole-key scan."""
    np, _, _, _ = _dependencies()
    if not len(target_keys):
        empty_rows = np.empty(0, dtype=np.int64)
        if return_keys:
            return empty_rows, np.empty(0, dtype=np.uint64)
        return empty_rows
    target_positions = target_keys >> np.uint64(4)
    frame_starts = np.fromiter(
        (int(frame["base_position"]) for frame in key_index),
        dtype=np.uint64,
        count=len(key_index),
    )
    frame_stops = np.fromiter(
        (int(frame["last_position"]) for frame in key_index),
        dtype=np.uint64,
        count=len(key_index),
    )
    candidate_ids = set()
    for position in np.unique(target_positions):
        frame_id = int(np.searchsorted(frame_stops, position, side="left"))
        while frame_id < len(key_index) and frame_starts[frame_id] <= position:
            candidate_ids.add(frame_id)
            frame_id += 1
    hits = []
    key_hits = []
    for frame_id in sorted(candidate_ids):
        frame = key_index[frame_id]
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
                hits.append(np.asarray(
                    valid_locations[matched] + int(frame["row_start"]),
                    dtype=np.int64,
                ))
                key_hits.append(np.asarray(candidates[valid][matched], dtype=np.uint64))
    if not hits:
        empty_rows = np.empty(0, dtype=np.int64)
        if return_keys:
            return empty_rows, np.empty(0, dtype=np.uint64)
        return empty_rows
    rows = np.concatenate(hits)
    matched_keys = np.concatenate(key_hits)
    order = np.argsort(rows, kind="stable")
    rows = rows[order]
    matched_keys = matched_keys[order]
    keep = np.concatenate((np.asarray([True]), rows[1:] != rows[:-1]))
    rows = rows[keep]
    matched_keys = matched_keys[keep]
    if return_keys:
        return rows, matched_keys
    return rows


def _frames_containing_rows(index, rows):
    """Return only indexed frames containing at least one sorted row ID."""
    np, _, _, _ = _dependencies()
    if not len(rows) or not index:
        return []
    row_ids = np.asarray(rows, dtype=np.int64)
    frame_stops = np.fromiter(
        (int(frame["row_stop"]) for frame in index),
        dtype=np.int64,
        count=len(index),
    )
    frame_ids = np.unique(np.searchsorted(frame_stops, row_ids, side="right"))
    selected = []
    for frame_id in frame_ids:
        if frame_id >= len(index):
            continue
        frame = index[int(frame_id)]
        if np.any(
            (row_ids >= int(frame["row_start"])) &
            (row_ids < int(frame["row_stop"]))
        ):
            selected.append(frame)
    return selected


def read_store_v2(
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
    _manifest=None,
    _key_index=None,
    _value_index=None,
):
    np, _, _, _ = _dependencies()
    manifest = load_manifest(store) if _manifest is None else _manifest
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
    key_index = _key_index
    if key_index is None:
        key_index = json.loads(
            (store / manifest["files"]["key_index"]).read_text()
        )
    resolved_identity_keys = None
    if identity_keys is not None:
        requested_keys = parse_identity_keys(identity_keys, manifest)
        requested_rows, resolved_identity_keys = _rows_for_identity_keys(
            store, manifest, key_index, requested_keys, return_keys=True
        )
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
        frames = _frames_containing_rows(key_index, requested_rows)
    if need_identity and resolved_identity_keys is not None:
        keys = resolved_identity_keys
        output_rows = requested_rows
    elif need_identity:
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
            store, manifest, streams=numeric_needed, value_index=_value_index
        )
        values = {}
    else:
        values = decode_values(
            store, manifest, row_start, row_stop,
            row_indices=output_rows if requested_rows is not None else None,
            needed=numeric_needed,
            value_index=_value_index,
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


def _combine_wrapped_keys(positions, substitutions):
    np, _, _, _ = _dependencies()
    positions = np.asarray(positions, dtype=np.uint64)
    substitutions = np.asarray(substitutions, dtype=np.uint64)
    if np.any(substitutions > 15):
        raise ValueError("wrapped substitution code exceeds four bits")
    return (positions << np.uint64(4)) | substitutions


def _wrapped_rows_for_keys(context: WrappedStoreReader, target_keys):
    np, _, _, _ = _dependencies()
    position_stream = context.stream("position")
    substitution_stream = context.stream("substitution")
    candidate_pages = set()
    for key in target_keys:
        page_id = int(np.searchsorted(position_stream.last_keys, key, side="left"))
        while (page_id < len(position_stream.pages) and
               position_stream.first_keys[page_id] <= key):
            candidate_pages.add(page_id)
            page_id += 1
    rows = []
    keys = []
    candidate_pages = sorted(candidate_pages)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        position_future = executor.submit(
            position_stream.read_pages, candidate_pages
        )
        substitution_future = executor.submit(
            substitution_stream.read_pages, candidate_pages
        )
        position_pages = position_future.result()
        substitution_pages = substitution_future.result()
    for page_id, position_page, substitution_page in zip(
            candidate_pages, position_pages, substitution_pages):
        positions = position_page.astype(np.uint64)
        substitutions = substitution_page.astype(np.uint64)
        decoded = _combine_wrapped_keys(positions, substitutions)
        lo = int(np.searchsorted(target_keys, position_stream.first_keys[page_id]))
        hi = int(np.searchsorted(
            target_keys, position_stream.last_keys[page_id], side="right"
        ))
        candidates = target_keys[lo:hi]
        locations = np.searchsorted(decoded, candidates)
        valid = locations < len(decoded)
        if np.any(valid):
            valid_locations = locations[valid]
            matched = decoded[valid_locations] == candidates[valid]
            if np.any(matched):
                rows.append(
                    valid_locations[matched].astype(np.int64) +
                    position_stream.page_starts[page_id]
                )
                keys.append(candidates[valid][matched])
    if not rows:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.uint64)
    rows = np.concatenate(rows)
    keys = np.concatenate(keys)
    order = np.argsort(rows, kind="stable")
    return rows[order], keys[order]


def _wrapped_rows_for_region(
    context: WrappedStoreReader, global_start: int, global_end: int
):
    np, _, _, _ = _dependencies()
    position_stream = context.stream("position")
    first_positions = position_stream.first_keys >> np.uint64(4)
    last_positions = position_stream.last_keys >> np.uint64(4)
    page_start = int(np.searchsorted(last_positions, global_start, side="left"))
    pages = []
    while page_start < len(first_positions) and first_positions[page_start] <= global_end:
        pages.append(page_start)
        page_start += 1
    rows = []
    keys = []
    substitution_stream = context.stream("substitution")
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        position_future = executor.submit(position_stream.read_pages, pages)
        substitution_future = executor.submit(
            substitution_stream.read_pages, pages
        )
        position_pages = position_future.result()
        substitution_pages = substitution_future.result()
    for page_id, position_page, substitution_page in zip(
            pages, position_pages, substitution_pages):
        positions = position_page.astype(np.uint64)
        substitutions = substitution_page.astype(np.uint64)
        keep = (positions >= global_start) & (positions <= global_end)
        if np.any(keep):
            local = np.flatnonzero(keep)
            rows.append(local.astype(np.int64) + position_stream.page_starts[page_id])
            keys.append(_combine_wrapped_keys(positions[keep], substitutions[keep]))
    if not rows:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.uint64)
    return np.concatenate(rows), np.concatenate(keys)


def _wrapped_exception_subset(exceptions, target_rows):
    np, _, _, _ = _dependencies()
    if not len(exceptions) or not len(target_rows):
        return exceptions[:0]
    exception_rows = exceptions["row"].astype(np.int64)
    locations = np.searchsorted(exception_rows, target_rows)
    valid = locations < len(exception_rows)
    matched_locations = locations[valid]
    matched = exception_rows[matched_locations] == target_rows[valid]
    return exceptions[matched_locations[matched]]


def _decode_wrapped_semantics(codes, target_rows, manifest, exceptions):
    np, _, _, _ = _dependencies()
    semantic = manifest["semantic_codec"]
    result = {}
    if "z" in codes:
        z_count = 2**int(semantic["z_bits"]) - 2
        z = np.full(len(target_rows), np.nan, dtype=np.float64)
        ok = codes["z"] < z_count
        z_min, z_max = semantic["z_range"]
        z[ok] = float(z_min) + (codes["z"][ok].astype(np.float64) + 0.5) * (
            (float(z_max) - float(z_min)) / z_count
        )
        result["z"] = z
    if "eaf" in codes:
        eaf_count = 2**int(semantic["eaf_bits"]) - 1
        result["eaf"] = np.sin(
            np.pi * codes["eaf"].astype(np.float64) / (2.0 * eaf_count)
        ) ** 2
    if "se" in codes:
        se_count = 2**int(semantic["se_bits"]) - 2
        se = np.full(len(target_rows), np.nan, dtype=np.float64)
        ok = codes["se"] < se_count
        centres = np.asarray(
            semantic["block_centers_log2_residual"], dtype=np.float64
        )
        block_rows = int(semantic["se_center_block_rows"])
        safe = np.clip(result["eaf"], 1e-12, 1.0 - 1e-12)
        residual = -1.0 + (codes["se"].astype(np.float64) + 0.5) * (
            2.0 / se_count
        )
        se[ok] = 2.0 ** (
            residual[ok] + centres[target_rows[ok] // block_rows] -
            0.5 * np.log2(2.0 * safe[ok] * (1.0 - safe[ok]))
        )
        result["se"] = se
    selected = _wrapped_exception_subset(exceptions, target_rows)
    if len(selected):
        local = np.searchsorted(target_rows, selected["row"].astype(np.int64))
        z_mask = (selected["flags"] & 1) != 0
        se_mask = (selected["flags"] & 2) != 0
        eaf_mask = (selected["flags"] & 4) != 0
        if "z" in result:
            result["z"][local[z_mask]] = selected["z"][z_mask]
        if "se" in result:
            result["se"][local[se_mask]] = np.exp2(
                selected["log2se"][se_mask].astype(np.float64)
            )
        if "eaf" in result:
            result["eaf"][local[eaf_mask]] = selected["eaf"][eaf_mask]
    return result


def read_store_v3(
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
    compact_identity_bridge: bool = True,
    manifest=None,
    context=None,
):
    np, _, _, _ = _dependencies()
    manifest = load_manifest(store) if manifest is None else manifest
    owned_context = context is None
    context = WrappedStoreReader(store, manifest) if context is None else context
    if os.environ.get("COMPRESSOR_RELOAD_EXCEPTIONS", "").strip().lower() in {
            "1", "true", "yes", "on"}:
        context._exceptions = None
    bytes_read_at_start = context.bytes_read
    try:
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
        requested_rows = None if row_indices is None else np.asarray(
            sorted(set(int(row) for row in row_indices)), dtype=np.int64
        )
        if requested_rows is not None and np.any(
            (requested_rows < 0) | (requested_rows >= n)
        ):
            raise ValueError("row indices are outside the store")
        resolved_keys = None
        if identity_keys is not None:
            requested_keys = parse_identity_keys(identity_keys, manifest)
            requested_rows, resolved_keys = _wrapped_rows_for_keys(context, requested_keys)
        global_region = None
        if chromosome is not None or start is not None or end is not None:
            if chromosome is None or start is None or end is None:
                raise ValueError("chromosome, start, and end must be supplied together")
            global_region = _region_global_bounds(
                manifest, chromosome, int(start), int(end)
            )
            region_rows, region_keys = _wrapped_rows_for_region(context, *global_region)
            if requested_rows is None:
                requested_rows, resolved_keys = region_rows, region_keys
            else:
                keep = np.isin(requested_rows, region_rows)
                requested_rows = requested_rows[keep]
                resolved_keys = region_keys[np.isin(region_rows, requested_rows)]
        row_start = 0 if row_start is None else max(0, int(row_start))
        row_stop = n if row_stop is None else min(n, int(row_stop))
        full_binary_selection = (
            output_format == "binary" and requested_rows is None and
            row_start == 0 and row_stop == n
        )
        if requested_rows is None:
            output_rows = (
                None if full_binary_selection else
                np.arange(row_start, row_stop, dtype=np.int64)
            )
        else:
            output_rows = requested_rows[
                (requested_rows >= row_start) & (requested_rows < row_stop)
            ]
            if resolved_keys is not None and len(resolved_keys) != len(output_rows):
                resolved_keys = None
        identity_columns = {
            "chromosome", "base_pair_location", "effect_allele", "other_allele"
        }
        need_identity = bool(identity_columns.intersection(requested_columns))
        full_identity_bridge = (
            full_binary_selection and need_identity and resolved_keys is None and
            compact_identity_bridge
        )
        identity_codes = None
        if need_identity:
            if full_identity_bridge:
                identity_codes = {
                    "position": context.stream("position").read_rows(),
                    "substitution": context.stream("substitution").read_rows(),
                }
                chromosome_codes = np.empty(0, dtype=np.uint8)
                positions = np.empty(0, dtype=np.int32)
                effect_codes = np.empty(0, dtype=np.uint8)
                other_codes = np.empty(0, dtype=np.uint8)
            else:
                if resolved_keys is None:
                    positions = context.stream("position").read_rows(output_rows).astype(np.uint64)
                    substitutions = context.stream("substitution").read_rows(output_rows).astype(np.uint64)
                    keys = _combine_wrapped_keys(positions, substitutions)
                else:
                    keys = resolved_keys
                chromosome_codes, positions, effect_codes, other_codes = _keys_to_columns(
                    keys, manifest
                )
        else:
            chromosome_codes = np.empty(0, dtype=np.uint8)
            positions = np.empty(0, dtype=np.int32)
            effect_codes = np.empty(0, dtype=np.uint8)
            other_codes = np.empty(0, dtype=np.uint8)
        numeric_needed = set()
        if any(name in requested_columns for name in ("z", "beta", "p_value")):
            numeric_needed.add("z")
        if any(name in requested_columns for name in ("standard_error", "beta")):
            numeric_needed.add("se")
        if "effect_allele_frequency" in requested_columns:
            numeric_needed.add("eaf")
        if "se" in numeric_needed:
            numeric_needed.add("eaf")
        sparse_parallel = output_rows is not None and len(numeric_needed) > 1
        if sparse_parallel:
            stream_names = sorted(numeric_needed)
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=len(stream_names) + 1) as executor:
                code_futures = {
                    name: executor.submit(
                        context.stream(name).read_rows, output_rows
                    )
                    for name in stream_names
                }
                exception_future = executor.submit(context.exceptions)
                codes = {
                    name: code_futures[name].result() for name in stream_names
                }
                exceptions = exception_future.result()
        else:
            codes = {
                name: context.stream(name).read_rows(output_rows)
                for name in numeric_needed
            }
            exceptions = context.exceptions() if numeric_needed else None
        full_code_bridge = (
            full_binary_selection and bool(numeric_needed)
        )
        values = {} if full_code_bridge else _decode_wrapped_semantics(
            codes, output_rows, manifest, exceptions if exceptions is not None else []
        )
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
        if not full_code_bridge and "beta" in requested_columns:
            output_values["beta"] = values["z"] * values["se"]
        if not full_code_bridge and "p_value" in requested_columns:
            output_values["p_value"] = np.asarray([
                math.erfc(abs(float(value)) / math.sqrt(2.0))
                if math.isfinite(float(value)) else float("nan")
                for value in values["z"]
            ])
        if output_format == "binary":
            _write_binary_bridge(
                output, n if output_rows is None else output_rows,
                requested_columns, output_values,
                semantic_codes=codes if full_code_bridge else None,
                identity_codes=identity_codes,
                exceptions=exceptions if full_code_bridge else None,
                manifest=manifest,
                source_bytes_read=context.bytes_read - bytes_read_at_start,
            )
            return
        if output_format != "tsv":
            raise ValueError(f"unsupported output format: {output_format}")
        chromosome_names = list(CHROM_LENGTHS)
        bases = np.asarray(["A", "C", "G", "T"])
        output_values["chromosome"] = [
            chromosome_names[index - 1] for index in chromosome_codes
        ]
        output_values["effect_allele"] = bases[effect_codes]
        output_values["other_allele"] = bases[other_codes]
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["row", *requested_columns])
            for index, row in enumerate(output_rows):
                writer.writerow([
                    int(row),
                    *[
                        int(output_values[name][index])
                        if name == "base_pair_location" else output_values[name][index]
                        for name in requested_columns
                    ],
                ])
    finally:
        if owned_context:
            context.close()


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
    compact_identity_bridge: bool = True,
    _manifest=None,
    _key_index=None,
    _value_index=None,
):
    manifest = load_manifest(store) if _manifest is None else _manifest
    if manifest["format_version"] == VERSION:
        context = _key_index if isinstance(_key_index, WrappedStoreReader) else None
        return read_store_v3(
            store, output, row_start=row_start, row_stop=row_stop,
            row_indices=row_indices, identity_keys=identity_keys,
            chromosome=chromosome, start=start, end=end, include_p=include_p,
            columns=columns, output_format=output_format,
            compact_identity_bridge=compact_identity_bridge,
            manifest=manifest, context=context,
        )
    return read_store_v2(
        store, output, row_start=row_start, row_stop=row_stop,
        row_indices=row_indices, identity_keys=identity_keys,
        chromosome=chromosome, start=start, end=end, include_p=include_p,
        columns=columns, output_format=output_format,
        _manifest=manifest, _key_index=_key_index, _value_index=_value_index,
    )


def _load_cached_store(store: Path, cache=None):
    """Load validated store metadata, reusing it while the manifest is unchanged."""
    manifest_path = store / "manifest.json"
    checksum_path = store / MANIFEST_CHECKSUM
    fingerprint = tuple(
        (path.stat().st_mtime_ns, path.stat().st_size)
        for path in (manifest_path, checksum_path)
    )
    cache_key = str(store)
    if cache is not None and cache_key in cache:
        cached_fingerprint, manifest, key_index, value_index = cache.pop(cache_key)
        if cached_fingerprint == fingerprint:
            cache[cache_key] = (
                cached_fingerprint, manifest, key_index, value_index
            )
            return manifest, key_index, value_index
        if isinstance(key_index, WrappedStoreReader):
            key_index.close()
    manifest = load_manifest(store)
    if manifest["format_version"] == VERSION:
        key_index = WrappedStoreReader(store, manifest)
        value_index = key_index
    else:
        key_index = json.loads(
            (store / manifest["files"]["key_index"]).read_text()
        )
        value_index = json.loads(
            (store / manifest["files"]["value_index"]).read_text()
        )
    if cache is not None:
        cache[cache_key] = (fingerprint, manifest, key_index, value_index)
        # A wrapped context can hold five stream descriptors.  Keep enough
        # contexts for ordinary reuse without approaching macOS's common
        # 256-descriptor soft limit in large GWAS grids.
        while len(cache) > 8:
            _, _, evicted, _ = cache.pop(next(iter(cache)))
            if isinstance(evicted, WrappedStoreReader):
                evicted.close()
    return manifest, key_index, value_index


def read_stores_batch_payload(
    payload: dict[str, Any],
    output: Path,
    threads: int = 1,
    cache=None,
    write_response: bool = True,
):
    """Execute many sparse canonical reads in one Python runtime."""
    if threads < 1:
        raise ValueError("batch threads must be at least one")
    reads = payload.get("reads") if isinstance(payload, dict) else None
    if not isinstance(reads, list) or not reads:
        raise ValueError("batch request must contain a non-empty reads list")
    coalesce = payload.get("coalesce", True)
    reload_contexts = payload.get("reload_contexts", False)
    if not isinstance(coalesce, bool) or not isinstance(reload_contexts, bool):
        raise ValueError("batch coalesce and reload_contexts must be booleans")
    output.mkdir(parents=True, exist_ok=False)
    if cache is None:
        cache = {}
    prepared = []
    bridge_for = {}
    seen_requests = {}
    for index, item in enumerate(reads):
        if not isinstance(item, dict):
            raise ValueError("every batch read must be an object")
        try:
            store = Path(item["store"]).expanduser().resolve(strict=True)
        except (KeyError, OSError) as exc:
            raise ValueError(f"batch read {index} has an invalid store") from exc
        if not store.is_dir():
            raise ValueError(f"batch read {index} store is not a directory")
        keys = item.get("keys")
        columns = item.get("columns")
        if isinstance(keys, str):
            keys = [keys]
        if isinstance(columns, str):
            columns = [columns]
        if not isinstance(keys, list) or any(not isinstance(key, str) for key in keys):
            raise ValueError(f"batch read {index} keys must be a string list")
        if (not isinstance(columns, list) or not columns or
                any(not isinstance(column, str) or not column for column in columns)):
            raise ValueError(f"batch read {index} columns must be non-empty strings")
        cache_key = str(store)
        manifest, key_index, value_index = _load_cached_store(
            store, None if reload_contexts else cache
        )
        signature = (cache_key, tuple(keys), tuple(columns))
        if coalesce and signature in seen_requests:
            bridge_for[index] = seen_requests[signature]
        else:
            seen_requests[signature] = index
            bridge_for[index] = index
            prepared.append(
                (index, store, keys, columns, manifest, key_index, value_index)
            )

    def execute(item):
        index, store, keys, columns, manifest, key_index, value_index = item
        bridge = output / str(index)
        try:
            if isinstance(key_index, WrappedStoreReader):
                with key_index.lock:
                    read_store(
                        store,
                        bridge,
                        identity_keys=keys,
                        columns=columns,
                        output_format="binary",
                        _manifest=manifest,
                        _key_index=key_index,
                        _value_index=value_index,
                    )
            else:
                read_store(
                    store,
                    bridge,
                    identity_keys=keys,
                    columns=columns,
                    output_format="binary",
                    _manifest=manifest,
                    _key_index=key_index,
                    _value_index=value_index,
                )
        except Exception as exc:
            raise RuntimeError(f"failed to read {store}: {exc}") from exc
        finally:
            retained = (not reload_contexts) and any(
                cached[2] is key_index for cached in cache.values()
            )
            if isinstance(key_index, WrappedStoreReader) and not retained:
                key_index.close()
        bridge_manifest = json.loads((bridge / "bridge.json").read_text())
        return {
            "index": index,
            "rows": int(bridge_manifest["rows"]),
            "source_bytes_read": int(bridge_manifest.get("source_bytes_read", 0)),
        }

    workers = min(int(threads), len(prepared))
    if workers == 1:
        results = [execute(item) for item in prepared]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(execute, prepared))
    unique_results = {item["index"]: item for item in results}
    results = [
        {
            "index": index,
            "bridge": bridge_for[index],
            "rows": unique_results[bridge_for[index]]["rows"],
        }
        for index in range(len(reads))
    ]
    response = {
        "format": "CompreSSoR-batch-bridge",
        "version": 1,
        "reads": results,
        "source_bytes_read": sum(
            int(item.get("source_bytes_read", 0)) for item in unique_results.values()
        ),
    }
    if write_response:
        temporary = output / "batch.json.tmp"
        temporary.write_text(json.dumps(response, indent=2) + "\n")
        os.replace(temporary, output / "batch.json")
    return response


def read_stores_batch(request_path: Path, output: Path, threads: int = 1):
    payload = json.loads(request_path.read_text())
    return read_stores_batch_payload(payload, output, threads=threads)


def _close_wrapped_cache(cache) -> None:
    seen = set()
    for _, _, key_index, _ in cache.values():
        if isinstance(key_index, WrappedStoreReader) and id(key_index) not in seen:
            key_index.close()
            seen.add(id(key_index))


def serve() -> None:
    """Serve newline-delimited read requests for one R session."""
    _dependencies()
    cache = {}
    print(json.dumps({
        "format": "CompreSSoR-worker",
        "version": 1,
        "ready": True,
    }), flush=True)
    for line in sys.stdin:
        if not line.strip():
            continue
        request_id = None
        try:
            request = json.loads(line)
            if not isinstance(request, dict):
                raise ValueError("worker request must be an object")
            request_id = request.get("id")
            command = request.get("command")
            if command == "shutdown":
                response = {"id": request_id, "ok": True}
                print(json.dumps(response, separators=(",", ":")), flush=True)
                _close_wrapped_cache(cache)
                return
            if command == "batch-read":
                output = Path(request["output"]).expanduser()
                try:
                    batch_result = read_stores_batch_payload(
                        {
                            "reads": request.get("reads"),
                            "coalesce": request.get("coalesce", True),
                        },
                        output,
                        threads=int(request.get("threads", 1)),
                        cache=cache,
                        write_response=True,
                    )
                except Exception as exc:
                    output.mkdir(parents=True, exist_ok=True)
                    temporary = output / "batch-error.json.tmp"
                    temporary.write_text(json.dumps({
                        "error": f"{type(exc).__name__}: {exc}",
                    }) + "\n")
                    os.replace(temporary, output / "batch-error.json")
                    raise
                result = {
                    "format": "CompreSSoR-batch-ack",
                    "version": 1,
                    "batch_file": "batch.json",
                    "reads": len(batch_result["reads"]),
                }
            elif command == "read":
                store = Path(request["store"]).expanduser().resolve(strict=True)
                output = Path(request["output"]).expanduser()
                manifest, key_index, value_index = _load_cached_store(store, cache)
                rows = request.get("rows")
                keys = request.get("keys")
                columns = request.get("columns")
                if isinstance(rows, int):
                    rows = [rows]
                if isinstance(keys, str):
                    keys = [keys]
                if isinstance(columns, str):
                    columns = [columns]
                read_store(
                    store,
                    output,
                    row_start=request.get("row_start"),
                    row_stop=request.get("row_stop"),
                    row_indices=rows,
                    identity_keys=keys,
                    chromosome=request.get("chromosome"),
                    start=request.get("start"),
                    end=request.get("end"),
                    columns=columns,
                    output_format="binary",
                    compact_identity_bridge=bool(
                        request.get("compact_identity", True)
                    ),
                    _manifest=manifest,
                    _key_index=key_index,
                    _value_index=value_index,
                )
                bridge_manifest = json.loads((output / "bridge.json").read_text())
                result = {"rows": int(bridge_manifest["rows"])}
            else:
                raise ValueError(f"unsupported worker command: {command}")
            response = {"id": request_id, "ok": True, "result": result}
        except Exception as exc:
            response = {
                "id": request_id,
                "ok": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        print(json.dumps(response, separators=(",", ":")), flush=True)
    _close_wrapped_cache(cache)


def _validate_index(index, rows: int, label: str, block_rows: int):
    if rows and not index:
        raise ValueError(f"{label} index is empty")
    expected = 0
    for frame_number, frame in enumerate(index):
        start = int(frame["row_start"])
        stop = int(frame["row_stop"])
        if start != expected or stop <= start or stop > rows:
            raise ValueError(f"{label} index does not cover contiguous valid rows")
        frame_rows = stop - start
        is_final = frame_number == len(index) - 1
        if (not is_final and frame_rows != block_rows) or frame_rows > block_rows:
            raise ValueError(f"{label} index frame size differs from the manifest")
        expected = stop
    if expected != rows:
        raise ValueError(f"{label} index ends at {expected}, expected {rows}")


def validate_store_v2(store: Path, full: bool = False):
    np, _, _, _ = _dependencies()
    manifest = load_manifest(store)
    rows = int(manifest["rows"])
    for relative in set(manifest["files"].values()):
        _verify_file(store, manifest, relative)
    key_index = json.loads((store / manifest["files"]["key_index"]).read_text())
    value_index = json.loads((store / manifest["files"]["value_index"]).read_text())
    _validate_index(key_index, rows, "key", int(manifest["key_block_rows"]))
    _validate_index(value_index, rows, "value", int(manifest["value_block_rows"]))
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
        if (np.any(exceptions["flags"] == 0) or
                np.any((exceptions["flags"].astype(np.uint16) & np.uint16(0xF8)) != 0)):
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


def validate_store_v3(store: Path, full: bool = False):
    np, _, _, _ = _dependencies()
    manifest = load_manifest(store)
    rows = int(manifest["rows"])
    for relative in set(manifest["files"].values()):
        _verify_file(store, manifest, relative)
    context = WrappedStoreReader(store, manifest)
    try:
        streams = [context.stream(name) for name in WRAPPED_STREAM_DTYPES]
        page_layout = [
            (stream.page_starts.tolist(), stream.page_stops.tolist())
            for stream in streams
        ]
        if any(layout != page_layout[0] for layout in page_layout[1:]):
            raise ValueError("wrapped stream page layouts disagree")
        position_stream = context.stream("position")
        if (len(position_stream.first_keys) != len(position_stream.pages) or
                np.any(position_stream.first_keys > position_stream.last_keys) or
                np.any(position_stream.first_keys[1:] <= position_stream.last_keys[:-1])):
            raise ValueError("wrapped identity page bounds are invalid")
        exceptions = context.exceptions()
        if len(exceptions):
            exception_rows = exceptions["row"].astype(np.int64)
            if (np.any(exception_rows < 0) or np.any(exception_rows >= rows) or
                    np.any(exception_rows[1:] <= exception_rows[:-1]) or
                    np.any(exceptions["flags"] == 0) or
                    np.any((exceptions["flags"].astype(np.uint16) & np.uint16(0xF8)) != 0)):
                raise ValueError("wrapped exception records are invalid")
        frames_checked = 0
        if full:
            decoded = {
                name: context.stream(name).read_rows()
                for name in WRAPPED_STREAM_DTYPES
            }
            frames_checked = sum(
                len(context.stream(name).pages) for name in WRAPPED_STREAM_DTYPES
            )
            keys = _combine_wrapped_keys(
                decoded["position"], decoded["substitution"]
            )
            substitutions = decoded["substitution"].astype(np.uint64)
            refs = substitutions >> np.uint64(2)
            alts = substitutions & np.uint64(3)
            genome_length = sum(CHROM_LENGTHS.values())
            if (len(keys) != rows or np.any(keys[1:] <= keys[:-1]) or
                    np.any(decoded["position"] >= genome_length) or np.any(refs == alts)):
                raise ValueError("wrapped identity stream violates the GRCh38 contract")
            if np.any(decoded["z"] > 511) or np.any(decoded["se"] > 63):
                raise ValueError("wrapped semantic code exceeds its bit width")
            z_flags = np.zeros(rows, dtype=bool)
            se_flags = np.zeros(rows, dtype=bool)
            if len(exceptions):
                exception_rows = exceptions["row"].astype(np.int64)
                z_flags[exception_rows] = (exceptions["flags"] & 1) != 0
                se_flags[exception_rows] = (exceptions["flags"] & 2) != 0
            if (not np.array_equal(decoded["z"] == 511, z_flags) or
                    not np.array_equal(decoded["se"] == 63, se_flags)):
                raise ValueError("wrapped exception flags and sentinels disagree")
        return {
            "valid": True,
            "errors": [],
            "rows": rows,
            "profile": manifest["profile"],
            "full": bool(full),
            "files_checked": len(set(manifest["files"].values())),
            "frames_checked": frames_checked,
        }
    finally:
        context.close()


def validate_store(store: Path, full: bool = False):
    manifest = load_manifest(store)
    if manifest["format_version"] == VERSION:
        return validate_store_v3(store, full=full)
    return validate_store_v2(store, full=full)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    write = subparsers.add_parser("write")
    write.add_argument("--input", type=Path, required=True)
    write.add_argument("--output", type=Path, required=True)
    write.add_argument("--metadata", type=Path)
    write.add_argument("--key-block-rows", type=int, default=WRAPPED_PAGE_ROWS)
    write.add_argument("--value-block-rows", type=int, default=WRAPPED_PAGE_ROWS)
    write.add_argument("--se-center-block-rows", type=int, default=SE_CENTER_BLOCK_ROWS)
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
    read.add_argument("--expanded-identity-bridge", action="store_true")
    batch_read = subparsers.add_parser("batch-read")
    batch_read.add_argument("--request", type=Path, required=True)
    batch_read.add_argument("--output", type=Path, required=True)
    batch_read.add_argument("--threads", type=int, default=1)
    subparsers.add_parser("serve")
    validate = subparsers.add_parser("validate")
    validate.add_argument("--store", type=Path, required=True)
    validate.add_argument("--full", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "write":
        metadata = json.loads(args.metadata.read_text()) if args.metadata else None
        print(json.dumps(write_store(
            args.input,
            args.output,
            metadata=metadata,
            key_block_rows=args.key_block_rows,
            value_block_rows=args.value_block_rows,
            se_center_block_rows=args.se_center_block_rows,
        ), indent=2))
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
                   output_format=args.output_format,
                   compact_identity_bridge=not args.expanded_identity_bridge)
    elif args.command == "batch-read":
        print(json.dumps(
            read_stores_batch(args.request, args.output, threads=args.threads),
            indent=2,
        ))
    elif args.command == "serve":
        serve()
    else:
        print(json.dumps(validate_store(args.store, full=args.full), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
