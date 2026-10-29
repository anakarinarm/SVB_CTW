#!/usr/bin/env python3
"""Create non-wrapping forcing files for one-piece 15-day CTW/SVB runs.

MITgcm's external forcing and OBCS prescribed fields share the
periodicExternalForcing machinery in this setup. To run a 15-day experiment
with time-varying winds and steady boundary T/S, use a 16-day forcing cycle:
the model stops before the periodic interpolation wraps back to record 1.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "input"
COMMON = ROOT / "ctw_reruns" / "common_input"

NX, NY, NR = 512, 612, 100
REC_SPACING_SECONDS = 1800
RUN_DAYS = 15
FORCING_CYCLE_DAYS = 16
FORCING_RECORDS = FORCING_CYCLE_DAYS * 86400 // REC_SPACING_SECONDS

WIND_FILES = {
    "upwelling_merid_halfgauss2hr_512x612.bin": "upwelling_merid_halfgauss2hr_16day_512x612.bin",
    "upwelling_zonal_halfgauss2hr_512x612.bin": "upwelling_zonal_halfgauss2hr_16day_512x612.bin",
}


def replace_with_tmp(path: Path, writer) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    writer(tmp)
    tmp.replace(path)


def extend_wind_file(src: Path, dst: Path) -> None:
    record_bytes = NX * NY * 8
    src_size = src.stat().st_size
    if src_size % record_bytes:
        raise ValueError(f"{src} size is not an integer number of 64-bit XY records")
    src_records = src_size // record_bytes
    if src_records > FORCING_RECORDS:
        raise ValueError(f"{src} has {src_records} records, more than target {FORCING_RECORDS}")

    def write(tmp: Path) -> None:
        data = np.memmap(src, dtype=">f8", mode="r", shape=(src_records, NY, NX))
        with tmp.open("wb") as out_f:
            for i in range(src_records):
                np.asarray(data[i], dtype=">f8").tofile(out_f)
            zero_record = b"\0" * record_bytes
            for _ in range(FORCING_RECORDS - src_records):
                out_f.write(zero_record)

    replace_with_tmp(dst, write)
    print(f"{dst.name}: copied {src_records} stress records, padded to {FORCING_RECORDS}")


def read_profile(name: str) -> np.ndarray:
    field = np.fromfile(COMMON / name, dtype=">f8").reshape(NR, NY, NX)
    profile = field[:, 0, 0]
    max_dev = float(np.max(np.abs(field - profile[:, None, None])))
    if max_dev > 0:
        raise ValueError(f"{name} is not horizontally homogeneous: max deviation {max_dev}")
    return profile.astype(">f8", copy=False)


def write_repeated_boundary(path: Path, record: np.ndarray) -> None:
    record_be = np.ascontiguousarray(record, dtype=">f8")

    def write(tmp: Path) -> None:
        with tmp.open("wb") as out_f:
            payload = record_be.tobytes(order="C")
            for _ in range(FORCING_RECORDS):
                out_f.write(payload)

    replace_with_tmp(path, write)


def write_obcs_timeseries() -> None:
    outdir = COMMON / "obcs_ts_from_initial_hydrography_16day"
    outdir.mkdir(exist_ok=True)
    theta = read_profile("feb_Tref_100zlev_512x612.bin")
    salt = read_profile("feb_Sref_100zlev_512x612.bin")
    fields = {
        "OBWt.bin": np.repeat(theta[:, None], NY, axis=1),
        "OBWs.bin": np.repeat(salt[:, None], NY, axis=1),
        "OBSt.bin": np.repeat(theta[:, None], NX, axis=1),
        "OBSs.bin": np.repeat(salt[:, None], NX, axis=1),
        "OBNt.bin": np.repeat(theta[:, None], NX, axis=1),
        "OBNs.bin": np.repeat(salt[:, None], NX, axis=1),
    }
    for name, record in fields.items():
        write_repeated_boundary(outdir / name, record)
        print(f"{outdir.name}/{name}: wrote {FORCING_RECORDS} steady records")


def main() -> None:
    COMMON.mkdir(parents=True, exist_ok=True)
    for src_name, dst_name in WIND_FILES.items():
        extend_wind_file(INPUT / src_name, COMMON / dst_name)
    write_obcs_timeseries()
    print(
        f"ready for {RUN_DAYS}-day run using {FORCING_CYCLE_DAYS}-day "
        f"cycle ({FORCING_RECORDS} records at {REC_SPACING_SECONDS}s)"
    )


if __name__ == "__main__":
    main()
