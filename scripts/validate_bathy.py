#!/usr/bin/env python3
"""Validate CTW/SVB final bathymetry files used by rerun directories."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RERUNS = ROOT / "ctw_reruns"
NX, NY = 512, 612
EXPECTED_SIZE = NX * NY * 8
EXPECTED_WET = {
    "plus": 211096,
    "minus": 204908,
}
README_SHA256 = {
    "bathy_realistic_minusSVB_FINAL.bin": "b9cd1468a972f6ce3b46b970a382bd79b5716acb27cd8279524e3a5a416ad0e1",
    "bathy_realistic_plusSVB_FINAL.bin": "5bd3f87847c6fa179953303795bf395a542dffc45f1e31b8249674e3e9293fa2",
    "bathy_smooth30_minusSVB_FINAL.bin": "38fad368107da97a3cb1d0caca723040fe329cc417b41a1af984908447e5808e",
    "bathy_smooth30_plusSVB_FINAL.bin": "a55f58aea88a393c9bbf8d7c2b62fd923d8511156a853030733e1f938f888fed",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    rows = []
    masks = {}
    for path in sorted((RERUNS / "bathymetry").glob("bathy_*_FINAL.bin")):
        arr = np.fromfile(path, dtype=">f8").reshape(NY, NX)
        wet = arr < 0
        kind = "plus" if "plusSVB" in path.name else "minus"
        west = int(wet[:, 0].sum())
        south = int(wet[0, :].sum())
        north = int(wet[-1, :].sum())
        digest = sha256(path.resolve())
        rows.append({
            "file": path.name,
            "size_ok": path.stat().st_size == EXPECTED_SIZE,
            "sha256": digest,
            "sha256_matches_readme": digest == README_SHA256[path.name],
            "wet_columns": int(wet.sum()),
            "wet_columns_ok": int(wet.sum()) == EXPECTED_WET[kind],
            "min_depth": float(arr.min()),
            "min_wet_depth": float(-arr[wet].max()),
            "land_zero_count": int((arr == 0).sum()),
            "west_wet_rows": west,
            "south_wet_cols": south,
            "north_wet_cols": north,
        })
        masks[path.name] = wet

    out = RERUNS / "validation_bathy.csv"
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    comparisons = {
        "Rplus_vs_Splus_mask_equal": np.array_equal(
            masks["bathy_realistic_plusSVB_FINAL.bin"],
            masks["bathy_smooth30_plusSVB_FINAL.bin"],
        ),
        "Rminus_vs_Sminus_mask_equal": np.array_equal(
            masks["bathy_realistic_minusSVB_FINAL.bin"],
            masks["bathy_smooth30_minusSVB_FINAL.bin"],
        ),
    }
    for key, val in comparisons.items():
        print(f"{key}: {val}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
