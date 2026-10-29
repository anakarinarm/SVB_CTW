#!/usr/bin/env python3
"""Check that CTW/SVB rerun configs differ only in approved places."""

from __future__ import annotations

import difflib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "ctw_reruns"
RUNS = ["Rplus", "Rminus", "Splus", "Sminus"]
FILES = ["data", "data.obcs", "data.diagnostics", "data.pkg", "data.mnc", "data.exch2", "eedata"]


def normalized(run: str, file_name: str) -> list[str]:
    text = (ROOT / run / file_name).read_text().splitlines()
    out = []
    for line in text:
        stripped = line.strip()
        if stripped.startswith("bathyFile="):
            out.append(" bathyFile='<BATHY>',")
        elif stripped.startswith("mnc_outdir_str="):
            out.append(" mnc_outdir_str='<MNC_OUTDIR>',")
        else:
            out.append(line)
    return out


def main() -> None:
    report = []
    ok = True
    base = RUNS[0]
    for run in RUNS[1:]:
        for name in FILES:
            diff = list(difflib.unified_diff(
                normalized(base, name),
                normalized(run, name),
                fromfile=f"{base}/{name}",
                tofile=f"{run}/{name}",
                lineterm="",
            ))
            if diff:
                ok = False
                report.extend(diff)
                report.append("")
    out = ROOT / "validation_config_diff.txt"
    out.write_text("\n".join(report) if report else "OK: only normalized bathyFile differs.\n")
    print(f"config_diff_ok: {ok}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
