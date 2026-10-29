#!/usr/bin/env python3
"""Prepare external-drive execution directories for CTW/SVB reruns.

Run MITgcm from these external directories so diagnostics, checkpoints, monitor
files, and any binary dumps are written to the external drive by default.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RERUNS = ROOT / "ctw_reruns"
RUNS = ["Rplus", "Rminus", "Splus", "Sminus", "S50plus", "S50minus", "S70plus", "S70minus"]
RUN_FILES = [
    "data",
    "data.obcs",
    "data.diagnostics",
    "data.pkg",
    "eedata",
]
TOP_LEVEL_LINKS = [
    "common_input",
    "bathymetry",
    "telescoped_input",
    "south_extended_input",
    "west_south_input",
]


def link_force(dst: Path, target: Path) -> None:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    rel = os.path.relpath(target, dst.parent)
    dst.symlink_to(rel)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reruns-root",
        type=Path,
        default=DEFAULT_RERUNS,
        help="Local run-family directory containing Rplus/Rminus/Splus/Sminus and optional S50 runs.",
    )
    parser.add_argument(
        "--external-root",
        type=Path,
        default=Path("/Volumes/LaCie/CTW_SVB_reruns"),
        help="External-drive directory that will contain the per-case run folders.",
    )
    parser.add_argument(
        "--executable",
        type=Path,
        default=None,
        help="MITgcm executable to symlink into each external run directory if present.",
    )
    args = parser.parse_args()

    reruns = args.reruns_root.resolve()
    executable = args.executable or (reruns / "build" / "mitgcmuv")
    external_root = args.external_root
    external_root.mkdir(parents=True, exist_ok=True)
    for name in TOP_LEVEL_LINKS:
        target = reruns / name
        if target.exists() or target.is_symlink():
            link_force(external_root / name, target)

    for run in RUNS:
        out = external_root / run
        out.mkdir(parents=True, exist_ok=True)
        for name in RUN_FILES:
            target = reruns / run / name
            if target.exists() or target.is_symlink():
                link_force(out / name, target)
        for name in TOP_LEVEL_LINKS:
            target = reruns / name
            if target.exists() or target.is_symlink():
                link_force(out / name, target)
        if executable.exists():
            link_force(out / executable.name, executable)

    (external_root / "README_EXTERNAL_RUNS.txt").write_text(
        "Run MITgcm from the per-case directories in this folder so output lands on the external drive.\n"
        "Example after building/linking the executable:\n"
        "  cd /Volumes/LaCie/CTW_SVB_reruns/Rplus\n"
        "  mpirun -np 16 ./mitgcmuv\n"
    )
    print(f"prepared external run directories under {external_root}")


if __name__ == "__main__":
    main()
