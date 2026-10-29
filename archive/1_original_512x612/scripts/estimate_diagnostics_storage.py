#!/usr/bin/env python3
"""Estimate storage for the requested CTW/SVB diagnostics."""

from __future__ import annotations

NX, NY = 512, 612
LEVELS_3D = list(range(1, 66)) + [67, 69, 71, 73, 75, 77, 79, 84, 90, 100]
N_3D_FIELDS = 5
N_2D_FIELDS = 4
SECONDS = 1296000
CADENCE = 1200
BYTES_PER_VALUE = 8


def gib(n: int) -> float:
    return n / 1024**3


def main() -> None:
    records = SECONDS // CADENCE
    grid = NX * NY
    bytes_3d = records * grid * len(LEVELS_3D) * N_3D_FIELDS * BYTES_PER_VALUE
    bytes_2d = records * grid * N_2D_FIELDS * BYTES_PER_VALUE
    per_run = bytes_3d + bytes_2d
    print(f"records_per_run: {records}")
    print(f"3d_levels: {len(LEVELS_3D)}")
    print(f"estimated_binary_gib_per_run: {gib(per_run):.2f}")
    print(f"estimated_binary_gib_four_runs: {gib(per_run * 4):.2f}")
    print("note: filesystem metadata, tiled files, checkpoints, and monitor output are extra")


if __name__ == "__main__":
    main()
