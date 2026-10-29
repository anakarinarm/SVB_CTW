#!/usr/bin/env python3
"""Print basic sanity checks for CTW rerun OBC T/S time-series files."""

from __future__ import annotations

from pathlib import Path

import numpy as np


BASE = Path(__file__).resolve().parents[1] / "input" / "common_input_512x612" / "obcs_ts_from_initial_hydrography_16day"
COMMON = Path(__file__).resolve().parents[1] / "input" / "common_input_512x612"
SPECS = {
    "OBWt.bin": (768, 100, 612),
    "OBWs.bin": (768, 100, 612),
    "OBSt.bin": (768, 100, 512),
    "OBSs.bin": (768, 100, 512),
    "OBNt.bin": (768, 100, 512),
    "OBNs.bin": (768, 100, 512),
}


def main() -> None:
    theta = np.fromfile(COMMON / "feb_Tref_100zlev_512x612.bin", ">f8").reshape(100, 612, 512)[:, 0, 0]
    salt = np.fromfile(COMMON / "feb_Sref_100zlev_512x612.bin", ">f8").reshape(100, 612, 512)[:, 0, 0]
    for name, shape in SPECS.items():
        path = BASE / name
        values = np.fromfile(path, ">f8")
        expected = int(np.prod(shape))
        print(f"{name}: bytes={path.stat().st_size} values={values.size} expected={expected}")
        records = values.reshape(shape)
        for idx in (0, 1, 2, shape[0] - 1):
            rec = records[idx]
            print(
                f"  rec {idx + 1:03d}: min={np.nanmin(rec):.10g} "
                f"max={np.nanmax(rec):.10g} mean={np.nanmean(rec):.10g}"
            )
        print(f"  max |rec2-rec1| = {np.max(np.abs(records[1] - records[0])):.3e}")
        print(f"  max |last-rec1| = {np.max(np.abs(records[-1] - records[0])):.3e}")
        profile = theta if name.endswith("t.bin") else salt
        expected = np.repeat(profile[:, None], shape[2], axis=1)
        print(f"  max |rec1-profile| = {np.max(np.abs(records[0] - expected)):.3e}")


if __name__ == "__main__":
    main()
